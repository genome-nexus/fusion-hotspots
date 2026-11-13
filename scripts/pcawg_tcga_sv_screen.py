"""Reproducible PCAWG TCGA DNA screen; run from repository root.

uv run python scripts/pcawg_tcga_sv_screen.py --stage screen|flags|sequence
Downloads are supplied as data in data/raw/pcawg_sv (see run REPORT.md).
Each stage finishes in the foreground and writes only inside this worktree.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import sys
import time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any

import pandas as pd
import requests

from cfh.cohort.artefact_flags import (
    EnsemblSequenceSource,
    ParalogLookup,
    genomic_flags,
    paralog_flag,
    prefetch_sequences,
)
from cfh.cohort.auto_config import batch_fetch_canonical_transcripts
from cfh.cohort.domain_screen import screen_cohort
from cfh.domain_function.pfam_hmm import fetch_pfam_release
from cfh.domain_function.record import TRANSLATED_CLASSES, build_records
from cfh.domain_function.run import (
    assess_all,
    load_hmms,
    load_structures,
    tier_counts,
    write_records,
    write_tsv,
)
from cfh.domain_function.uniprot import UniProtClient
from cfh.ingestion.pcawg_sv import pair_genes, read_archive, read_gencode, sample_index
from cfh.mapping.junction_domain_classifier import RECONSTITUTION_CLASSES
from cfh.mapping.transcript_models import TranscriptModelSource

csv.field_size_limit(sys.maxsize)
RAW = Path("data/raw/pcawg_sv")
WORK = Path("runs/.pcawg_tcga_sv_work")
OUT = Path("runs/pcawg_tcga_sv_screen_2026-11-13")
CACHE = WORK / "transcripts_GRCh37"
COHORT = "pcawg_tcga_sv"
EXTRA = (
    "donor_unique_id",
    "icgc_donor_id",
    "submitter_donor_id",
    "aliquot_id",
    "icgc_specimen_id",
    "submitter_specimen_id",
    "donor_wgs_exclusion_white_gray",
    "dcc_project_code",
    "sv_id",
    "svclass",
    "strand1",
    "strand2",
    "svmethod",
    "caller_support",
    "pe_support",
    "chrom1",
    "start1",
    "end1",
    "chrom2",
    "start2",
    "end2",
    "source_file",
    "source_line",
    "gene_pairing_status",
    "five_gene_id",
    "three_gene_id",
    "five_gene_strand",
    "three_gene_strand",
)
NONTRIVIAL = {
    *TRANSLATED_CLASSES,
    "out_of_frame",
    "premature_stop_at_junction",
    "splice_site_unresolved",
}


def log(message: str) -> None:
    print(time.strftime("%H:%M:%S"), message, flush=True)


def write_rows(path: Path, rows: list[dict[str, Any]]) -> None:
    columns = list(dict.fromkeys(key for row in rows for key in row))
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_rows(path: Path) -> list[dict[str, Any]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def save_summary(**updates: Any) -> None:
    path = OUT / "summary.json"
    current = json.loads(path.read_text()) if path.exists() else {}
    path.write_text(json.dumps({**current, **updates}, indent=2, sort_keys=True) + "\n")


class PcawgTranscripts(TranscriptModelSource):
    """Bound public batch lookups and preserve missing mappings as unassessed.

    Use the same MSKCC GRCh37 transcript endpoint and build-stamped cache.
    Batch responses are matched by returned hugoSymbols, which may differ
    from the requested symbol. Omissions stay unset until the base class's
    single-gene canonical lookup runs. Only its definitive not-found response
    may write a negative cache marker; request errors never do.
    """

    def prefetch_canonical(self, genes: list[str]) -> None:
        pending = []
        for gene in sorted(set(genes)):
            marker = self._not_found_marker(gene)
            if marker is not None and marker.exists():
                text = marker.read_text().strip()
                if text == "no canonical mapping in successful batch":
                    marker.unlink()  # invalidate the previous runner's poisoned omission marker
                elif text == "404":
                    self._canonical[gene] = None
                else:
                    raise ValueError(f"unrecognized negative transcript cache marker: {marker}")
            if gene not in self._canonical:
                pending.append(gene)
        cache = self._canonical_cache_dir
        chunks = [pending[i : i + 100] for i in range(0, len(pending), 100)]
        errors = {}

        def fetch(chunk):
            try:
                result = batch_fetch_canonical_transcripts(
                    chunk,
                    base_url=self.client.base_url,
                    cache_dir=cache,
                    batch_size=100,
                    max_retries=1,
                    timeout=30,
                    isoform_override_source="mskcc",
                )
                return chunk, result, None
            except requests.RequestException as error:
                return chunk, None, str(error)

        with ThreadPoolExecutor(max_workers=4) as pool:
            for index, (chunk, result, error) in enumerate(pool.map(fetch, chunks), 1):
                found = result.by_gene_symbol if result else {}
                for gene in chunk:
                    if gene in found:
                        self._canonical[gene] = found[gene]
                    if error:
                        errors[gene] = error
                if index % 5 == 0:
                    log(f"canonical batches {index}/{len(chunks)}")
        omitted = [gene for gene in pending if gene not in self._canonical]
        log(f"single-gene fallback for {len(omitted)} batch omissions/unavailable results")
        # canonical() retains the base implementation's positive/definitive
        # negative cache behavior. A failed GET raises without caching None.
        with ThreadPoolExecutor(max_workers=4) as pool:
            for index, _ in enumerate(pool.map(self.canonical, omitted), 1):
                if index % 200 == 0:
                    log(f"single-gene fallback {index}/{len(omitted)}")
        save_summary(
            canonical_batch_errors=errors,
            canonical_lookup_errors={},
            canonical_fallback_requested=len(omitted),
            canonical_fallback_resolved=sum(self._canonical[g] is not None for g in omitted),
            canonical_fallback_not_found=sum(self._canonical[g] is None for g in omitted),
        )


def annotation_problem(row: dict[str, Any], source: TranscriptModelSource) -> str | None:
    """Require canonical identity/strand agreement with the breakend gene model."""
    for side, label in (("five_prime", "five"), ("three_prime", "three")):
        gene = row[f"{side}_gene"]
        transcript = source.canonical(gene)
        if transcript is None:
            return f"no canonical transcript mapping for {gene}"
        directory = source._canonical_cache_dir
        payload = json.loads((directory / f"{gene.upper()}.json").read_text()) if directory else {}
        observed = str(payload.get("geneId") or "").split(".")[0]
        expected = row[f"{label}_gene_id"].split(".")[0]
        if not observed or observed != expected:
            return f"{gene}: canonical gene ID {observed or 'missing'} disagrees with {expected}"
        expected_strand = 1 if row[f"{label}_gene_strand"] == "+" else -1
        if {e.strand for e in transcript.exons} != {expected_strand}:
            return f"{gene}: canonical exon strand missing or inconsistent with breakend gene"
    return None


def model_outcome_category(row: dict[str, Any]) -> str:
    """Separate missing data from negative/ambiguous base-screen product verdicts."""
    if row["reconstitution_class"] != "no_transcript_model":
        return "not_applicable"
    problem = row.get("annotation_problem") or ""
    reason = row.get("pre_annotation_reason") or row.get("class_reason") or ""
    if problem.startswith("no canonical transcript mapping"):
        return "missing_canonical_mapping"
    if problem:
        return "annotation_identity_or_strand_mismatch"
    if "ambiguous genomic breakpoints" in reason:
        return "base_screen_ambiguous_breakpoints"
    if "no site breakpoint within target locus" in reason or "no breakpoint within" in reason:
        return "base_screen_outside_canonical_locus"
    if any(
        text in reason
        for text in ("past its stop codon", "outside the transcript", "before its coding region")
    ):
        return "base_screen_no_coding_product"
    return "missing_or_invalid_transcript_model"


def screen() -> None:
    if hashlib.md5((RAW / "consensus.tgz").read_bytes()).hexdigest() != (
        "650981d737522350503bba74e6d2c65d"
    ):
        raise ValueError("PCAWG archive does not match the published MD5")
    samples = sample_index(RAW / "pcawg_sample_sheet.tsv")
    raw, summary = read_archive(RAW / "consensus.tgz", samples)
    log(f"ingested {len(raw)} SVs")
    genes = read_gencode(RAW / "gencode_comp.txt.gz", RAW / "gencode_attrs.txt.gz")
    pairs, audits = [], []
    for row in raw:
        candidates, audit = pair_genes(row, genes)
        pairs.extend(candidates)
        audits.append(audit)
    counts = Counter(row["pairing_status"] for row in audits)
    write_rows(
        OUT / "breakend_gene_audit.tsv.gz",
        [
            {
                key: row[key]
                for key in (
                    *EXTRA[:17],
                    "chrom1",
                    "start1",
                    "end1",
                    "chrom2",
                    "start2",
                    "end2",
                    "genes1",
                    "genes2",
                    "pair_count",
                    "pairing_status",
                )
            }
            for row in audits
        ],
    )
    donors = {samples[a]["donor_unique_id"]: samples[a] for a in summary["file_row_counts"]}
    eligible = {
        r["donor_unique_id"]: r
        for r in samples.values()
        if r["dcc_project_code"].endswith("-US")
        and r["donor_wgs_exclusion_white_gray"] in {"Whitelist", "Graylist"}
    }
    write_rows(
        OUT / "donors_without_bedpe.tsv",
        [eligible[d] for d in sorted(eligible.keys() - donors.keys())],
    )
    write_rows(
        OUT / "sample_inventory.tsv",
        [
            {**samples[a], "sv_rows": count, "archive_present": True}
            for a, count in sorted(summary["file_row_counts"].items())
        ],
    )
    save_summary(
        ingest=summary,
        donor_status_counts=dict(
            Counter(r["donor_wgs_exclusion_white_gray"] for r in donors.values())
        ),
        gene_pairing_counts=dict(counts),
        sense_gene_pairs=len(pairs),
        donors_with_calls=len({r["donor_unique_id"] for r in raw}),
        eligible_sample_sheet_donors=len(eligible),
        sv_type_counts=dict(Counter(r["svclass"] for r in raw)),
    )
    log(f"{len(pairs)} sense gene pairs; {dict(counts)}")
    source = PcawgTranscripts("GRCh37", cache_dir=CACHE)
    rows, result = screen_cohort(
        pd.DataFrame(pairs), source, cohort=COHORT, passthrough_columns=EXTRA, progress=log
    )
    for row in rows:
        row["cohort"] = COHORT
        row["patient_id"] = row["icgc_donor_id"]
        row["cancer_type"] = row["dcc_project_code"].removesuffix("-US")
        # The shared-family helper labels a missing canonical as 'no_shared_family'.
        # Keep absence of annotation distinct from a negative domain screen.
        row["pre_annotation_class"] = row["reconstitution_class"]
        row["pre_annotation_reason"] = row["class_reason"]
        problem = annotation_problem(row, source)
        row["annotation_problem"] = problem or ""
        if problem:
            row["reconstitution_class"] = "no_transcript_model"
            row["assessment_status"] = "not_assessed"
            row["class_reason"] = problem
            row["no_product_reason"] = problem
        else:
            row["assessment_status"] = (
                "not_assessed"
                if row["reconstitution_class"] == "no_transcript_model"
                else "assessed"
            )
        row["model_outcome_category"] = model_outcome_category(row)
    write_rows(OUT / "domain_screen.tsv.gz", rows)
    counts_by_class = Counter(row["reconstitution_class"] for row in rows)
    result["class_counts"] = {
        cls: counts_by_class[cls] for cls in ("no_shared_family", *RECONSTITUTION_CLASSES)
    }
    shared = [r for r in rows if r["pfam_family"]]
    result["family_rows"] = len(shared)
    result["events_with_a_shared_family"] = len(
        {(r["aliquot_id"], r["sv_id"], r["five_prime_gene"], r["three_prime_gene"]) for r in shared}
    )
    unassessed = [r for r in rows if r["reconstitution_class"] == "no_transcript_model"]
    ambiguous = [
        r
        for r in unassessed
        if r["pfam_family"] and "ambiguous genomic breakpoints" in r["pre_annotation_reason"]
    ]
    save_summary(
        screen=result,
        model_outcome_counts=dict(Counter(r["model_outcome_category"] for r in unassessed)),
        ambiguous_breakpoint_shared_rows=len(ambiguous),
        ambiguous_breakpoint_svs=len({(r["aliquot_id"], r["sv_id"]) for r in ambiguous}),
    )
    log(f"screen complete: {result['class_counts']}")


class SafeSequences(EnsemblSequenceSource):
    """Record failed lookups as unavailable, never as negative homology."""

    def region(self, chromosome: str, start: int, end: int) -> str | None:
        try:
            return super().region(chromosome, start, end)
        except requests.RequestException as error:
            log(f"reference unavailable {chromosome}:{start}-{end}: {error}")
            return None

    def transcript_cds(self, transcript_id: str) -> str | None:
        try:
            return super().transcript_cds(transcript_id)
        except requests.RequestException as error:
            log(f"CDS unavailable {transcript_id}: {error}")
            return None

    def paralog_lookup(self, gene: str, transcript_id: str | None = None) -> ParalogLookup:
        try:
            return super().paralog_lookup(gene, transcript_id)
        except requests.RequestException as error:
            return ParalogLookup({}, False, (f"{gene}: {error}",))


def donor_burdens(rows: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    """Known shared-family SV burden; unannotated pairs make low counts incomplete."""
    calls: dict[str, set] = defaultdict(set)
    unknown: set[str] = set()
    donors = {r["donor_unique_id"] for r in rows}
    for row in rows:
        donor = row["donor_unique_id"]
        if row.get("pfam_family"):
            calls[donor].add((row["aliquot_id"], row["sv_id"]))
        elif row["reconstitution_class"] == "no_transcript_model":
            unknown.add(donor)
    return {
        donor: {
            "sample_shared_family_calls": len(calls[donor]),
            "donor_family_lookup_status": "not_assessed" if donor in unknown else "assessed",
            "sample_burden_flag": (
                "flagged"
                if len(calls[donor]) >= 2
                else "not_assessed"
                if donor in unknown
                else "pass"
            ),
            "sample_burden_reason": (
                f"{len(calls[donor])} distinct known shared-family SVs; "
                + (
                    "additional pairs lack family annotation"
                    if donor in unknown
                    else "family lookups complete"
                )
            ),
        }
        for donor in sorted(donors)
    }


def flags() -> None:
    rows = read_rows(OUT / "domain_screen.tsv.gz")
    candidates = [
        r
        for r in rows
        if r["reconstitution_class"] not in {"no_shared_family", "no_transcript_model"}
    ]
    seq = SafeSequences("GRCh37", WORK / "flags_GRCh37", timeout=30, max_retries=2)
    prefetch_sequences(candidates, {"GRCh37": seq}, workers=4, progress=log)
    burden = donor_burdens(rows)
    paralog_assessed: dict[str, set] = defaultdict(set)
    paralog_burden: dict[str, set] = defaultdict(set)
    for row in candidates:
        key = (row["aliquot_id"], row["sv_id"])
        row.update(paralog_flag(row, seq))
        if row["paralog_sequence_flag"] in {"pass", "flagged"}:
            paralog_assessed[row["donor_unique_id"]].add(key)
        if row["paralog_sequence_flag"] == "flagged":
            paralog_burden[row["donor_unique_id"]].add(key)
    for index, row in enumerate(candidates, 1):
        row.update(
            genomic_flags(
                row,
                seq,
                None,
                burden[row["donor_unique_id"]]["sample_shared_family_calls"],
                "GRCh37 segdup and RepeatMasker tracks not supplied",
            )
        )
        row.update(burden[row["donor_unique_id"]])
        row.update(
            donor_paralog_burden_status=(
                "not_assessed"
                if burden[row["donor_unique_id"]]["donor_family_lookup_status"] == "not_assessed"
                or len(paralog_assessed[row["donor_unique_id"]])
                < burden[row["donor_unique_id"]]["sample_shared_family_calls"]
                else "assessed"
            ),
            normal_recurrence_flag="not_applicable",
            normal_recurrence_reason="RNA recurrence check is not applicable to DNA calls",
            donor_paralog_calls=len(paralog_burden[row["donor_unique_id"]]),
        )
        if index % 50 == 0:
            log(f"flagged {index}/{len(candidates)}")
    write_rows(
        OUT / "donor_burden.tsv",
        [
            {
                **sample,
                **burden.get(sample["donor_unique_id"], {}),
                "donor_paralog_calls": len(paralog_burden[sample["donor_unique_id"]]),
                "paralog_assessed_calls": len(paralog_assessed[sample["donor_unique_id"]]),
                "paralog_burden_status": "not_assessed"
                if burden.get(sample["donor_unique_id"], {}).get("donor_family_lookup_status")
                != "assessed"
                or len(paralog_assessed[sample["donor_unique_id"]])
                < burden.get(sample["donor_unique_id"], {}).get("sample_shared_family_calls", 0)
                else "assessed",
            }
            for sample in read_rows(OUT / "sample_inventory.tsv")
        ],
    )
    write_rows(OUT / "shared_family_flags.tsv.gz", candidates)
    write_rows(
        OUT / "nontrivial_candidates.tsv",
        [r for r in candidates if r["reconstitution_class"] in NONTRIVIAL],
    )
    save_summary(
        flag_rows=len(candidates),
        flag_counts={
            key: dict(Counter(r[key] for r in candidates))
            for key in (
                "paralog_sequence_flag",
                "genomic_similarity_flag",
                "breakpoint_repeat_flag",
                "sample_burden_flag",
                "read_evidence_flag",
                "normal_recurrence_flag",
            )
        },
    )
    log("flags complete")


DNA_REQUIRED_FLAGS = (
    "paralog_sequence_flag",
    "genomic_similarity_flag",
    "read_evidence_flag",
    "breakpoint_repeat_flag",
    "sample_burden_flag",
)


def dna_tier(row: dict[str, Any], sequence_tier: str) -> str:
    """RNA recurrence is separate from the DNA evidence gate."""
    if row.get("gene_pairing_status") != "unique" or any(
        row.get(flag) != "pass" for flag in DNA_REQUIRED_FLAGS
    ):
        return "not_assessed"
    return sequence_tier


def sequence() -> None:
    rows = [
        r
        for r in read_rows(OUT / "shared_family_flags.tsv.gz")
        if r["reconstitution_class"] in NONTRIVIAL
    ]
    records, problems = build_records(
        rows, {"GRCh37": TranscriptModelSource("GRCh37", cache_dir=CACHE)}, progress=log
    )
    log(f"{len(records)} records; problems {dict(problems)}")
    hmms, hmm_problems = load_hmms(
        {r.pfam_family for r in records}, WORK / "pfam_hmm", progress=log
    )
    uniprot = UniProtClient(WORK / "uniprot")
    genes = sorted({g for r in records for g in (r.five.gene, r.three.gene)})
    uniprot.prefetch(genes)
    entries = {gene: uniprot.entry(gene) for gene in genes}
    structures, structure_problems = load_structures(
        {e.accession for e in entries.values() if e}, WORK / "afdb", progress=log
    )
    metrics = assess_all(
        records, hmms, entries, structures, progress=log, structure_problems=structure_problems
    )
    provenance = {
        "cohort": COHORT,
        "genome_build": "GRCh37",
        "input_shared_family_flags_sha256": hashlib.sha256(
            (OUT / "shared_family_flags.tsv.gz").read_bytes()
        ).hexdigest(),
        "databases": {
            **fetch_pfam_release(),
            "uniprot_release": uniprot.release,
            "afdb_model_versions": sorted({s.model_version for s in structures.values() if s}),
        },
    }
    (OUT / "sequence_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    write_records(records, OUT / "chimera_records.json.gz", provenance)
    write_tsv(metrics, OUT / "domain_function_metrics.tsv.gz")
    by_id = {r.record_id: r for r in records}
    lookup = {}
    for metric in metrics:
        for obs in by_id[metric["record_id"]].observations:
            lookup[
                (
                    obs["sample_id"],
                    str(obs["five_prime_breakpoint"]),
                    str(obs["three_prime_breakpoint"]),
                    metric["pfam_family"],
                    metric["fusion"],
                )
            ] = metric
    for row in rows:
        for column in ("translated_frame_status", "alternative_convention_frame_status"):
            row[column] = row.get(column) or "not_assessed"
        key = (
            row["sample_id"],
            row["five_prime_breakpoint"],
            row["three_prime_breakpoint"],
            row["pfam_family"],
            row["fusion_name"],
        )
        metric = lookup.get(key)
        row["sequence_tier"] = metric["tier"] if metric else "not_assessed"
        row["sequence_record_id"] = metric["record_id"] if metric else ""
        row["tier"] = dna_tier(row, row["sequence_tier"])

    write_rows(OUT / "nontrivial_candidates.tsv", rows)
    save_summary(
        candidate_tiers=dict(Counter(r["tier"] for r in rows)),
        sequence_tiers=tier_counts(metrics, records),
        candidate_sequence_tiers=dict(Counter(r["sequence_tier"] for r in rows)),
        sequence_problems=dict(problems),
        hmm_problems=hmm_problems,
        structure_problems=structure_problems,
        sequence_records=len(records),
    )
    log("sequence complete")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("screen", "flags", "sequence"), required=True)
    args = parser.parse_args()
    WORK.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    {"screen": screen, "flags": flags, "sequence": sequence}[args.stage]()
