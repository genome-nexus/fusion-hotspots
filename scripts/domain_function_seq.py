"""Run sequence-level functional-domain checks on the domain-screen pilot outputs.

Usage::

    python scripts/domain_function_seq.py \\
        --screen runs/domain_screen_pilot_2026-11-10/combined_candidates.tsv.gz \\
        --grch38-cache <pilot cache>/target_starfusion_GRCh38 \\
        --grch37-cache <pilot cache>/msk_impact_50k_2026_GRCh37 \\
        --work runs/.domain_function_seq_work --out runs/domain_function_seq_<date>

Stages are cached in ``--work`` (records, HMMs, UniProt, AFDB-derived JSON),
so a rerun only repeats what changed.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from collections import Counter
from importlib.metadata import version
from pathlib import Path

from cfh.domain_function import SCHEMA_VERSION
from cfh.domain_function.afdb import (
    BURIED_MAX_RSA,
    CONFIDENT_PLDDT,
    CONTACT_DISTANCE,
    MIN_SEQUENCE_SEPARATION,
)
from cfh.domain_function.metrics import (
    JUNCTION_HALF_WINDOW,
    KEY_RESIDUE_PROXIMITY,
    MAX_SITE_SPAN,
    MIN_UNIPROT_IDENTITY,
)
from cfh.domain_function.pfam_hmm import fetch_pfam_release
from cfh.domain_function.record import FLANK, TRANSLATED_CLASSES, build_records
from cfh.domain_function.run import (
    assess_all,
    load_hmms,
    load_or_build_records,
    load_structures,
    read_screen_rows,
    records_fingerprint,
    shortlist,
    tier_counts,
    write_records,
    write_tsv,
)
from cfh.domain_function.tiers import REQUIRED_INPUTS, RULES, TIER_LABEL
from cfh.domain_function.uniprot import UniProtClient
from cfh.mapping.transcript_models import TranscriptModelSource


def log(message: str) -> None:
    print(time.strftime("%H:%M:%S"), message, flush=True)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--screen", type=Path, required=True)
    parser.add_argument("--grch38-cache", type=Path, required=True)
    parser.add_argument("--grch37-cache", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    args.out.mkdir(parents=True, exist_ok=True)

    log("reading screen rows")
    rows = read_screen_rows(args.screen)
    cache_dirs = {"GRCh38": args.grch38_cache, "GRCh37": args.grch37_cache}

    def build() -> tuple[list, Counter]:
        sources = {
            build: TranscriptModelSource(build, cache_dir=path)
            for build, path in cache_dirs.items()
        }
        return build_records(rows, sources, progress=log)

    records, problems, fingerprint, reused = load_or_build_records(
        args.work / "records.pickle",
        lambda: records_fingerprint(args.screen, cache_dirs),
        build,
        progress=log,
    )
    log(
        f"{len(records)} records ({'reused cache' if reused else 'rebuilt'}); "
        f"rebuild problems {dict(problems)}"
    )

    cache = args.work / "cache"
    hmms, hmm_problems = load_hmms(
        {r.pfam_family for r in records}, cache / "pfam_hmm", progress=log
    )
    log(f"HMMs: {sum(h is not None for h in hmms.values())}/{len(hmms)}")

    uniprot = UniProtClient(cache / "uniprot")
    genes = sorted({g for r in records for g in (r.five.gene, r.three.gene)})
    uniprot.prefetch(genes)
    entries = {gene: uniprot.entry(gene) for gene in genes}
    log(f"UniProt: {sum(e is not None for e in entries.values())}/{len(genes)} genes")

    accessions = sorted({e.accession for e in entries.values() if e is not None})
    structures, afdb_problems = load_structures(accessions, cache / "afdb", progress=log)
    log(f"AFDB: {sum(s is not None for s in structures.values())}/{len(accessions)}")

    metrics = assess_all(
        records, hmms, entries, structures, progress=log, structure_problems=afdb_problems
    )
    counts = tier_counts(metrics, records)
    short = shortlist(metrics)

    release = fetch_pfam_release()
    afdb_versions = sorted({s.model_version for s in structures.values() if s is not None})
    git_sha = subprocess.run(
        ["git", "rev-parse", "HEAD"], capture_output=True, text=True, check=False
    ).stdout.strip()
    dirty = subprocess.run(
        ["git", "status", "--porcelain", "--", "src", "scripts", "pyproject.toml"],
        capture_output=True,
        text=True,
        check=False,
    ).stdout.strip()
    provenance = {
        "generated": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "code_sha": git_sha,
        "code_dirty": bool(dirty),
        "input_screen": str(args.screen),
        "input_screen_sha256": fingerprint["screen_sha256"],
        "records_fingerprint": fingerprint,
        "records_cache_reused": reused,
        "translated_classes": list(TRANSLATED_CLASSES),
    }
    write_records(records, args.out / "chimera_records.json.gz", provenance)
    write_tsv(metrics, args.out / "domain_function_metrics.tsv.gz")
    write_tsv(short, args.out / "shortlist_for_structure.tsv")
    (args.out / "tier_counts.json").write_text(json.dumps(counts, indent=1) + "\n")
    manifest = {
        "schema_version": SCHEMA_VERSION,
        **provenance,
        "tier_kind": TIER_LABEL,
        "tier_rules": RULES,
        "tier_required_inputs_for_likely_intact": list(REQUIRED_INPUTS),
        "records": len(records),
        "distinct_sequences": len({r.sequence_id for r in records}),
        "record_identity": (
            "record_id: family, class, parent transcripts, parent domain spans, junction "
            "residues, hybrid residue, chimeric domain span and flanked sequence; "
            "sequence_id: family, parent transcripts and flanked sequence"
        ),
        "read_support_rule": (
            "samples counted once per biological sample (both TARGET tables share one "
            "namespace); a sample's reads = max over its rows, never summed; "
            "total = sum over samples of per-sample max"
        ),
        "rebuild_problems": dict(problems),
        "screen_rows_in_translated_classes": sum(
            1 for r in rows if r["reconstitution_class"] in TRANSLATED_CLASSES
        ),
        "databases": {
            **release,
            "pfam_hmm_source": "InterPro API /entry/pfam/<acc>/?annotation=hmm",
            "pfam_families": len(hmms),
            "pfam_families_unavailable": hmm_problems,
            "uniprot_release": uniprot.release,
            "uniprot_genes": len(genes),
            "uniprot_genes_without_reviewed_entry": [g for g, e in entries.items() if e is None],
            "afdb_model_versions": afdb_versions,
            "afdb_accessions": len(accessions),
            "afdb_problems": afdb_problems,
            "afdb_missing": sorted(a for a, s in structures.items() if s is None),
        },
        "parameters": {
            "flank": FLANK,
            "contact_distance_A": CONTACT_DISTANCE,
            "min_sequence_separation": MIN_SEQUENCE_SEPARATION,
            "confident_plddt": CONFIDENT_PLDDT,
            "buried_max_rsa": BURIED_MAX_RSA,
            "key_residue_proximity": KEY_RESIDUE_PROXIMITY,
            "max_site_span": MAX_SITE_SPAN,
            "min_uniprot_identity": MIN_UNIPROT_IDENTITY,
            "junction_half_window": JUNCTION_HALF_WINDOW,
            "hmmer": "pyhmmer search_hmm, filters off (--max), Z=domZ=1, all hits reported",
        },
        "tools": {name: version(name) for name in ("pyhmmer", "gemmi", "biopython", "freesasa")},
    }
    (args.out / "run_manifest.json").write_text(json.dumps(manifest, indent=1) + "\n")
    log(f"done: tiers {counts['overall_records']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
