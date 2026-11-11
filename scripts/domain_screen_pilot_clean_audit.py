"""Audit the domain-screen pilot's clean_reconstitution calls.

Usage::

    python -I scripts/domain_screen_pilot_clean_audit.py RUN_DIR PEP_FASTA HOMOLOGIES_TSV

``RUN_DIR`` holds ``<cohort>.domain_screen.tsv(.gz)`` and ``flag_summary.json``.
``PEP_FASTA`` is Ensembl ``Homo_sapiens.GRCh38.pep.all.fa.gz`` (transcript ->
gene and symbol -> gene) and ``HOMOLOGIES_TSV`` the Ensembl Compara
``protein_default.homologies.tsv.gz`` for homo_sapiens (human paralogue pairs).
Writes ``RUN_DIR/clean_call_audit.json`` and ``RUN_DIR/clean_call_audit.md``.

A clean call requires a truncated copy of the family in both partners, in frame,
and a chimeric domain with exactly one parent's length and no indel at the
junction. By construction, then, the two pieces are stitched at the homologous
position of the two domains. The audit asks how such a stitch arose:

* ``same_gene``: 5' and 3' gene are the same symbol.
* ``readthrough``: same chromosome, same transcribed strand (open tables only)
  and breakpoints within ``READTHROUGH_BP``; a transcriptional read-through or
  local rearrangement of neighbouring genes.
* ``compara_paralog``: the two genes are human paralogues in Ensembl Compara.
* ``ref_splice``: STAR-Fusion ``ONLY_REF_SPLICE`` (both ends at annotated splice
  sites).
* ``minor_side``: one side contributes fewer than ``MIN_SIDE_RESIDUES`` residues
  or less than ``MIN_SIDE_FRACTION`` of its domain, i.e. the domain is in effect
  intact in the other partner.
* ``near_identical_parent``: the chimeric domain is at least
  ``NEAR_IDENTICAL`` identical to a parent's, i.e. at protein level it cannot be
  told apart from an unfused parent.

``strict`` = not same gene, not read-through, no minor side, and all three
artefact flags ``pass``.
"""

from __future__ import annotations

import csv
import gzip
import json
import re
import sys
from collections import Counter, defaultdict
from collections.abc import Iterator
from pathlib import Path
from typing import Any

csv.field_size_limit(sys.maxsize)

READTHROUGH_BP = 200_000
MIN_SIDE_RESIDUES = 10
MIN_SIDE_FRACTION = 0.10
NEAR_IDENTICAL = 0.95
FLAGS = ("paralog_sequence_flag", "genomic_similarity_flag", "normal_recurrence_flag")
_SPAN = re.compile(r"^(\d+)-(\d+) kept (\d+)-(\d+)$")


def _open(path: Path) -> Any:
    return gzip.open(path, "rt", newline="") if path.suffix == ".gz" else path.open(newline="")


def screen_rows(run_dir: Path, cohort: str) -> Iterator[dict[str, str]]:
    plain = run_dir / f"{cohort}.domain_screen.tsv"
    path = plain if plain.exists() else plain.with_suffix(".tsv.gz")
    with _open(path) as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def ensembl_gene_maps(pep: Path) -> tuple[dict[str, str], dict[str, set[str]]]:
    """Transcript (unversioned) -> gene, and symbol -> genes, from peptide headers."""
    transcript_gene: dict[str, str] = {}
    symbol_genes: dict[str, set[str]] = defaultdict(set)
    with gzip.open(pep, "rt") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            fields = dict(part.split(":", 1) for part in line[1:].split() if ":" in part)
            gene = fields.get("gene", "").split(".")[0]
            transcript = fields.get("transcript", "").split(".")[0]
            if gene and transcript:
                transcript_gene[transcript] = gene
            if gene and fields.get("gene_symbol"):
                symbol_genes[fields["gene_symbol"]].add(gene)
    return transcript_gene, symbol_genes


def human_paralogues(homologies: Path) -> set[frozenset[str]]:
    pairs: set[frozenset[str]] = set()
    with gzip.open(homologies, "rt", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if (
                row["species"] == row["homology_species"] == "homo_sapiens"
                and "paralog" in (row["homology_type"])
            ):
                pairs.add(frozenset((row["gene_stable_id"], row["homology_gene_stable_id"])))
    return pairs


def _side(span: str) -> tuple[int, int] | None:
    """(kept residues, domain residues) from ``"4-376 kept 4-329"``."""
    match = _SPAN.match(span or "")
    if not match:
        return None
    start, end, kept_start, kept_end = (int(value) for value in match.groups())
    return kept_end - kept_start + 1, end - start + 1


def audit_row(
    row: dict[str, str],
    transcript_gene: dict[str, str],
    symbol_genes: dict[str, set[str]],
    paralogues: set[frozenset[str]],
) -> dict[str, Any]:
    five, three = row["five_prime_gene"], row["three_prime_gene"]
    same_chrom = row["five_prime_chromosome"] == row["three_prime_chromosome"]
    try:
        distance = abs(int(row["five_prime_breakpoint"]) - int(row["three_prime_breakpoint"]))
    except ValueError:
        distance = None
    strands = (row.get("Open_strand_5p") or "", row.get("Open_strand_3p") or "")
    same_strand = bool(strands[0]) and strands[0] == strands[1]
    readthrough = (
        five != three
        and same_chrom
        and same_strand
        and distance is not None
        and distance <= READTHROUGH_BP
    )

    def genes(side: str, symbol: str) -> set[str]:
        transcript = (row.get(f"{side}_ensembl_transcript") or "").split(".")[0]
        if transcript in transcript_gene:
            return {transcript_gene[transcript]}
        return set(symbol_genes.get(symbol, set()))

    five_genes, three_genes = genes("five_prime", five), genes("three_prime", three)
    if five == three:
        paralog = "same_gene"
    elif not five_genes or not three_genes:
        paralog = "gene_not_in_ensembl_116"
    elif any(frozenset((a, b)) in paralogues for a in five_genes for b in three_genes):
        paralog = "compara_paralog"
    else:
        paralog = "not_paralog_same_family"

    five_side, three_side = (
        _side(row["five_prime_domain_span"]),
        _side(row["three_prime_domain_span"]),
    )
    minor = None
    if five_side and three_side:
        minor = any(
            kept < MIN_SIDE_RESIDUES or kept / total < MIN_SIDE_FRACTION
            for kept, total in (five_side, three_side)
        )
    try:
        identity = float(row["alignment_identity"])
    except ValueError:
        identity = None
    flags_pass = all(row[flag] == "pass" for flag in FLAGS)
    return {
        "same_gene": five == three,
        "readthrough": readthrough,
        "same_chrom_within_200kb": five != three
        and same_chrom
        and distance is not None
        and distance <= READTHROUGH_BP,
        "paralog": paralog,
        "ref_splice": row.get("Open_splice_type") or "not_reported",
        "minor_side": minor,
        "five_kept": five_side,
        "three_kept": three_side,
        "near_identical_parent": identity is not None and identity >= NEAR_IDENTICAL,
        "identity": identity,
        "flags_pass": flags_pass,
        "paralog_seq_flagged": row["paralog_sequence_flag"] == "flagged",
        "strict": five != three and not readthrough and minor is False and flags_pass,
    }


def junction_key(row: dict[str, str]) -> tuple[str, ...]:
    return (
        row["five_prime_gene"],
        row["three_prime_gene"],
        row["five_prime_chromosome"],
        row["five_prime_breakpoint"],
        row["three_prime_chromosome"],
        row["three_prime_breakpoint"],
        row["pfam_family"],
    )


def audit(run_dir: Path, pep: Path, homologies: Path) -> dict[str, Any]:
    transcript_gene, symbol_genes = ensembl_gene_maps(pep)
    paralogues = human_paralogues(homologies)
    cohorts = list(json.loads((run_dir / "flag_summary.json").read_text())["screen_summaries"])
    result: dict[str, Any] = {
        "definitions": {
            "readthrough_bp": READTHROUGH_BP,
            "min_side_residues": MIN_SIDE_RESIDUES,
            "min_side_fraction": MIN_SIDE_FRACTION,
            "near_identical": NEAR_IDENTICAL,
        },
        "cohorts": {},
        "strict_candidates": [],
    }
    for cohort in cohorts:
        rows = [
            r
            for r in screen_rows(run_dir, cohort)
            if r["reconstitution_class"] == "clean_reconstitution"
        ]
        counts: Counter[str] = Counter()
        identity_bins: Counter[str] = Counter()
        pairs: dict[tuple[str, str], dict[str, Any]] = {}
        junctions: set[tuple[str, ...]] = set()
        samples: set[str] = set()
        strict_junctions: dict[tuple[str, ...], list[dict[str, str]]] = defaultdict(list)
        for row in rows:
            info = audit_row(row, transcript_gene, symbol_genes, paralogues)
            junctions.add(junction_key(row))
            samples.add(row["sample_id"])
            counts["rows"] += 1
            for key in (
                "same_gene",
                "readthrough",
                "same_chrom_within_200kb",
                "near_identical_parent",
                "flags_pass",
                "paralog_seq_flagged",
                "strict",
            ):
                counts[key] += bool(info[key])
            counts[f"paralog:{info['paralog']}"] += 1
            counts[f"splice:{info['ref_splice']}"] += 1
            counts[
                "minor_side"
                if info["minor_side"]
                else ("both_sides_substantial" if info["minor_side"] is False else "side_unparsed")
            ] += 1
            identity = info["identity"]
            identity_bins[
                "unknown"
                if identity is None
                else ">=0.95"
                if identity >= 0.95
                else "0.80-0.95"
                if identity >= 0.80
                else "<0.80"
            ] += 1
            pair = (row["five_prime_gene"], row["three_prime_gene"])
            entry = pairs.setdefault(
                pair,
                {
                    "rows": 0,
                    "samples": set(),
                    "junctions": set(),
                    "paralog": info["paralog"],
                    "readthrough": 0,
                    "families": set(),
                },
            )
            entry["rows"] += 1
            entry["samples"].add(row["sample_id"])
            entry["junctions"].add(junction_key(row))
            entry["readthrough"] += info["readthrough"]
            entry["families"].add(row["pfam_family"])
            if info["strict"]:
                strict_junctions[junction_key(row)].append(row)
        top = sorted(pairs.items(), key=lambda item: (-item[1]["rows"], item[0]))[:10]
        result["cohorts"][cohort] = {
            "counts": dict(counts),
            "distinct_junctions": len(junctions),
            "distinct_samples": len(samples),
            "distinct_gene_pairs": len(pairs),
            "chimera_parent_identity": dict(identity_bins),
            "strict_distinct_junctions": len(strict_junctions),
            "top_pairs": [
                {
                    "pair": f"{a}::{b}",
                    "rows": e["rows"],
                    "samples": len(e["samples"]),
                    "junctions": len(e["junctions"]),
                    "paralog": e["paralog"],
                    "readthrough_rows": e["readthrough"],
                    "families": sorted(e["families"]),
                }
                for (a, b), e in top
            ],
        }
        for key, members in strict_junctions.items():
            best = max(members, key=lambda r: int(float(r["tumor_split_reads"] or 0)))
            info = audit_row(best, transcript_gene, symbol_genes, paralogues)
            result["strict_candidates"].append(
                {
                    "cohort": cohort,
                    "pair": f"{key[0]}::{key[1]}",
                    "family": key[6],
                    "breakpoints": f"{key[2]}:{key[3]} / {key[4]}:{key[5]}",
                    "samples": sorted(
                        {
                            f"{r['sample_id']} ({r.get('Open_tissue_class') or ''}"
                            f"{',' + r['Open_tumor_type'] if r.get('Open_tumor_type') else ''})"
                            for r in members
                        }
                    ),
                    "n_samples": len({r["sample_id"] for r in members}),
                    "max_split_reads": int(float(best["tumor_split_reads"] or 0)),
                    "reads": "; ".join(
                        f"{r['sample_id']} {r['tumor_split_reads']}/{r['tumor_paired_end_reads']}"
                        f"{' FFPM ' + r['Open_ffpm'] if r.get('Open_ffpm') else ''}"
                        f"{' FI ' + r['Open_fi_validated'] if r.get('Open_fi_validated') else ''}"
                        for r in members
                    ),
                    "junction_peptide": best["junction_peptide"],
                    "frame": (
                        f"{best['translated_frame_status']} / alt "
                        f"{best['alternative_convention_frame_status']}"
                    ),
                    "five_kept": best["five_prime_domain_span"],
                    "three_kept": best["three_prime_domain_span"],
                    "identity_to_parent": best["alignment_identity"],
                    "alignment_parent": best["alignment_parent"],
                    "paralog_relation": info["paralog"],
                    "splice_type": best.get("Open_splice_type"),
                    "flags": {flag: best[flag] for flag in FLAGS},
                    "flag_reasons": {
                        "paralog": best["paralog_sequence_reason"],
                        "genomic": best["genomic_similarity_reason"],
                        "normal": best["normal_recurrence_reason"],
                    },
                }
            )
    result["strict_candidates"].sort(
        key=lambda c: (-c["n_samples"], -c["max_split_reads"], c["pair"])
    )
    return result


def markdown(result: dict[str, Any]) -> str:
    cohorts = [c for c, v in result["cohorts"].items() if v["counts"].get("rows")]
    out = []
    keys = [
        ("rows", "clean rows"),
        ("distinct_junctions", "distinct junctions (pair + breakpoints + family)"),
        ("distinct_samples", "distinct samples"),
        ("distinct_gene_pairs", "distinct gene pairs"),
        ("same_gene", "same gene (5' == 3')"),
        ("readthrough", f"read-through: same chrom + strand, <= {READTHROUGH_BP // 1000} kb"),
        ("same_chrom_within_200kb", "any different-gene pair on one chrom <= 200 kb"),
        ("paralog:compara_paralog", "Compara human paralogues"),
        ("paralog:not_paralog_same_family", "same Pfam family, not Compara paralogues"),
        ("paralog:gene_not_in_ensembl_116", "gene not in Ensembl 116 (retired symbol)"),
        ("splice:ONLY_REF_SPLICE", "both ends at annotated splice sites (ONLY_REF_SPLICE)"),
        ("splice:INCL_NON_REF_SPLICE", "includes a non-reference splice site"),
        (
            "minor_side",
            f"one side < {MIN_SIDE_RESIDUES} aa or < {MIN_SIDE_FRACTION:.0%} of its domain",
        ),
        ("both_sides_substantial", "both sides contribute substantially"),
        (
            "near_identical_parent",
            f"chimera >= {NEAR_IDENTICAL:.0%} identical to one parent domain",
        ),
        ("paralog_seq_flagged", "paralog_sequence_flag = flagged"),
        ("flags_pass", "all three artefact flags pass"),
        (
            "strict",
            "**strict**: not same gene, not read-through, both sides substantial, flags pass",
        ),
        ("strict_distinct_junctions", "**strict distinct junctions**"),
    ]
    out.append("| Measure | " + " | ".join(cohorts) + " |")
    out.append("|---|" + "---|" * len(cohorts))
    for key, label in keys:
        values = []
        for cohort in cohorts:
            data = result["cohorts"][cohort]
            value = data.get(key, data["counts"].get(key, 0))
            values.append(f"{value:,}")
        out.append(f"| {label} | " + " | ".join(values) + " |")
    out.append("")
    for cohort in cohorts:
        out.append(f"\nTop clean gene pairs, {cohort}:\n")
        out.append(
            "| Pair | Rows | Samples | Junctions | Relation | Read-through rows | Families |"
        )
        out.append("|---|---|---|---|---|---|---|")
        for entry in result["cohorts"][cohort]["top_pairs"]:
            out.append(
                f"| {entry['pair']} | {entry['rows']} | {entry['samples']} | "
                f"{entry['junctions']} | {entry['paralog']} | {entry['readthrough_rows']} | "
                f"{', '.join(entry['families'])} |"
            )
    return "\n".join(out) + "\n"


def main() -> None:
    run_dir, pep, homologies = (Path(arg) for arg in sys.argv[1:4])
    result = audit(run_dir, pep, homologies)
    (run_dir / "clean_call_audit.json").write_text(json.dumps(result, indent=1) + "\n")
    (run_dir / "clean_call_audit.md").write_text(markdown(result))
    print(markdown(result))


TIER2_CLASSES = (
    "clean_reconstitution",
    "ambiguous_junction_register",
    "tandem_duplication_or_insertion",
    "gap_or_deletion",
)


def tier2_candidates(run_dir: Path, pep: Path, homologies: Path) -> list[dict[str, Any]]:
    """Distinct junctions that pass every artefact check except sample burden.

    Kept: class in :data:`TIER2_CLASSES`, in frame under the translated
    convention, paralog sequence ``pass``, normal recurrence ``pass``, no
    breakpoint-window similarity, no segdup/repeat, split reads present,
    not read-through, no minor side. Sample burden is reported, not gated, because
    in the RNA tables almost every sample carries >= 2 shared-family calls.
    """
    transcript_gene, symbol_genes = ensembl_gene_maps(pep)
    paralogues = human_paralogues(homologies)
    cohorts = list(json.loads((run_dir / "flag_summary.json").read_text())["screen_summaries"])
    groups: dict[tuple[str, ...], list[dict[str, str]]] = defaultdict(list)
    for cohort in cohorts:
        for row in screen_rows(run_dir, cohort):
            if row["reconstitution_class"] not in TIER2_CLASSES:
                continue
            if row["translated_frame_status"] != "in_frame":
                continue
            if row["paralog_sequence_flag"] != "pass" or row["normal_recurrence_flag"] != "pass":
                continue
            if row["breakpoint_repeat_flag"] != "pass" or row["read_evidence_flag"] != "pass":
                continue
            similar = False
            for window in (150, 500):
                try:
                    identity = float(row[f"genomic_identity_{window}"] or 0)
                    length = int(float(row[f"genomic_aligned_length_{window}"] or 0))
                except ValueError:
                    continue
                similar |= identity >= 0.80 and length >= 100
            if similar:
                continue
            info = audit_row(row, transcript_gene, symbol_genes, paralogues)
            if info["readthrough"] or info["minor_side"] is not False:
                continue
            groups[(cohort, row["reconstitution_class"], *junction_key(row))].append(row)
    out = []
    for key, members in groups.items():
        best = max(members, key=lambda r: int(float(r["tumor_split_reads"] or 0)))
        info = audit_row(best, transcript_gene, symbol_genes, paralogues)
        out.append(
            {
                "cohort": key[0],
                "class": key[1],
                "pair": f"{key[2]}::{key[3]}",
                "family": key[8],
                "breakpoints": f"{key[4]}:{key[5]} / {key[6]}:{key[7]}",
                "n_samples": len({r["sample_id"] for r in members}),
                "samples": "; ".join(
                    f"{r['sample_id']} ({r.get('Open_tissue_class') or ''}"
                    f"{',' + r['Open_tumor_type'] if r.get('Open_tumor_type') else ''}) "
                    f"{r['tumor_split_reads']}/{r['tumor_paired_end_reads']}"
                    f"{' FFPM ' + r['Open_ffpm'] if r.get('Open_ffpm') else ''}"
                    f"{' FI ' + r['Open_fi_validated'] if r.get('Open_fi_validated') else ''}"
                    f" burden {r['sample_shared_family_calls']}"
                    for r in members
                ),
                "max_split_reads": int(float(best["tumor_split_reads"] or 0)),
                "junction_peptide": best["junction_peptide"],
                "indel": best["junction_indel_sequence"],
                "frame": f"{best['translated_frame_status']} / alt "
                f"{best['alternative_convention_frame_status']}",
                "five_kept": best["five_prime_domain_span"],
                "three_kept": best["three_prime_domain_span"],
                "identity_to_parent": best["alignment_identity"],
                "relation": info["paralog"],
                "same_chrom_within_200kb": info["same_chrom_within_200kb"],
                "splice_type": best.get("Open_splice_type") or "",
                "paralog_reason": best["paralog_sequence_reason"],
                "genomic_reason": best["genomic_similarity_reason"],
                "normal_reason": best["normal_recurrence_reason"],
            }
        )
    out.sort(key=lambda c: (-c["n_samples"], -c["max_split_reads"], c["pair"]))
    return out


def tier2_main() -> None:
    run_dir, pep, homologies = (Path(arg) for arg in sys.argv[2:5])
    rows = tier2_candidates(run_dir, pep, homologies)
    (run_dir / "tier2_candidates.json").write_text(json.dumps(rows, indent=1) + "\n")
    print(len(rows))
    for row in rows[:40]:
        print(
            json.dumps(
                {
                    k: row[k]
                    for k in (
                        "cohort",
                        "class",
                        "pair",
                        "family",
                        "n_samples",
                        "max_split_reads",
                        "junction_peptide",
                        "frame",
                        "five_kept",
                        "three_kept",
                        "identity_to_parent",
                        "relation",
                        "same_chrom_within_200kb",
                    )
                }
            )
        )


if __name__ == "__main__":
    if sys.argv[1:2] == ["--tier2"]:
        tier2_main()
    else:
        main()
