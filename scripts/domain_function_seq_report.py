"""Write REPORT.md tables for a domain_function_seq run directory.

Usage: python scripts/domain_function_seq_report.py runs/domain_function_seq_<date>
"""

from __future__ import annotations

import csv
import gzip
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

TIERS = ("likely_intact", "uncertain", "not_assessed", "likely_broken")
CLASSES = (
    "clean_reconstitution",
    "tandem_duplication_or_insertion",
    "gap_or_deletion",
    "ambiguous_junction_register",
)
NAMED = (
    "FLI1::ERG",
    "STAT5B::STAT5A",
    "GPC4::GPC3",
    "GPC6::GPC5",
    "GPC3::GPC4",
    "LPO::MPO",
    "PIM3::ULK1",
)


def read(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:  # type: ignore[operator]
        return list(csv.DictReader(handle, delimiter="\t"))


def table(header: list[str], rows: list[list[object]]) -> str:
    lines = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    lines += ["| " + " | ".join(str(c) for c in row) + " |" for row in rows]
    return "\n".join(lines)


def counts_table(counts: dict[str, dict[str, int]]) -> str:
    rows = []
    for group in sorted(counts):
        cohort, klass = group.split("|")
        tiers = counts[group]
        rows.append([cohort, klass, *(tiers.get(t, 0) for t in TIERS), sum(tiers.values())])
    return table(["cohort", "class", *TIERS, "total"], rows)


def flag_summary(row: dict[str, str]) -> str:
    flagged = []
    for column in (
        "paralog_sequence_flag",
        "genomic_similarity_flag",
        "breakpoint_repeat_flag",
        "sample_burden_flag",
        "read_evidence_flag",
        "normal_recurrence_flag",
    ):
        counts = dict(
            item.split(":") for item in row.get(f"{column}_counts", "").split(",") if item
        )
        if counts.get("flagged"):
            flagged.append(f"{column.replace('_flag', '')}={counts['flagged']}")
    return ", ".join(flagged) or "none flagged"


def keys(row: dict[str, str]) -> str:
    """Key-feature statuses, or why they could not be assessed."""
    if row.get("key_features_assessed") != "True":
        statuses = [
            f"{p}: {row.get(f'key_annotation_status_{p}', '?')}"
            for p in ("5p", "3p")
            if row.get(f"key_annotation_status_{p}")
            not in ("annotated", "no_annotations_in_domain")
        ]
        return "not assessed (" + ", ".join(statuses) + ")"
    parts = []
    for status in ("lost", "unresolved", "altered", "swapped", "conserved", "retained"):
        value = row.get(f"key_{status}", "")
        if value not in ("", "0"):
            parts.append(f"{status} {value}")
    if row.get("key_point_lost"):
        parts.append(f"POINT LOST {row['key_point_lost']}")
    if row.get("key_point_altered"):
        parts.append(f"point altered {row['key_point_altered']}")
    if row.get("key_span_lost_or_altered"):
        parts.append(f"span lost/altered {row['key_span_lost_or_altered']}")
    return "; ".join(parts) or "no annotated feature in window"


def reads(row: dict[str, str]) -> str:
    """``max one sample / total over samples (samples)``."""
    best = row.get("max_split_reads_one_sample") or "-"
    total = row.get("total_split_reads_over_samples") or "-"
    return f"{best} / {total} ({row.get('n_samples')})"


def bits_vs_ga(row: dict[str, str]) -> str:
    five = row.get("parent5_domain_bits") or "-"
    three = row.get("parent3_domain_bits") or "-"
    chimera = row.get("chimera_domain_bits") or "-"
    return f"{chimera} / GA {row.get('ga_domain', '-')[:5]} (5' {five}, 3' {three})"


def candidate_rows(rows: list[dict[str, str]]) -> list[list[object]]:
    out = []
    for row in rows:
        out.append(
            [
                row.get("shortlist_rank", ""),
                row["fusion"],
                row["pfam_family"] + " " + row.get("hmm_name", ""),
                row["junction_class"].replace("_", " "),
                row["tier"],
                bits_vs_ga(row),
                row.get("chimera_match_coverage", ""),
                f"{row.get('schema_E_5p', '')}/{row.get('schema_E_3p', '')}",
                row.get("indel_ss", ""),
                keys(row),
                reads(row),
                row["cohorts"].replace(";", ", "),
                flag_summary(row),
                row.get("shortlist_reason", ""),
            ]
        )
    return out


CANDIDATE_HEADER = [
    "#",
    "fusion",
    "family",
    "class",
    "tier",
    "bits vs GA (parents)",
    "cov",
    "E 5'/3'",
    "indel SS",
    "key residues",
    "split reads: best sample / sum over samples (samples)",
    "cohorts",
    "artefact flags (rows flagged)",
    "shortlist",
]


def main() -> int:
    run = Path(sys.argv[1])
    metrics = read(run / "domain_function_metrics.tsv.gz")
    short = read(run / "shortlist_for_structure.tsv")
    counts = json.loads((run / "tier_counts.json").read_text())
    manifest = json.loads((run / "run_manifest.json").read_text())
    out: list[str] = []

    out.append("## Tier counts (HEURISTIC)\n")
    by_class: dict[str, Counter[str]] = defaultdict(Counter)
    for row in metrics:
        by_class[row["junction_class"]][row["tier"]] += 1
    out.append(
        f"Distinct chimeric-domain records (sequence + origin), all cohorts; "
        f"{counts['distinct_sequences']} distinct sequences:\n"
    )
    out.append(
        table(
            ["class", *TIERS, "total"],
            [
                [k, *(by_class[k].get(t, 0) for t in TIERS), sum(by_class[k].values())]
                for k in CLASSES
            ]
            + [
                [
                    "**all**",
                    *(counts["overall_records"].get(t, 0) for t in TIERS),
                    len(metrics),
                ]
            ],
        )
    )
    out.append(
        "\nRecords per cohort and class (a record seen in several cohorts counts in each):\n"
    )
    out.append(counts_table(counts["records"]))
    out.append(
        "\nDistinct biological samples per cohort and class, by their record's tier (a "
        "TARGET sample in both TARGET tables appears in both rows here but once in record "
        "totals):\n"
    )
    out.append(counts_table(counts["samples"]))
    out.append("\nScreen rows (observations) per cohort and class, by their record's tier:\n")
    out.append(counts_table(counts["rows"]))

    rules: Counter[str] = Counter()
    missing: Counter[str] = Counter()
    disagreements: Counter[str] = Counter()
    for row in metrics:
        rules.update(r for r in row["tier_rules"].split(";") if r)
        missing.update(m for m in row["tier_missing_inputs"].split(";") if m)
        disagreements.update(d for d in row["metric_disagreements"].split(";") if d)
    out.append("\n## Rules fired, missing inputs, metric disagreements (records)\n")
    out.append(table(["rule", "records"], sorted(rules.items())))
    out.append("")
    out.append(table(["missing input", "records"], sorted(missing.items())))
    out.append("")
    out.append(table(["disagreement", "records"], sorted(disagreements.items())))
    for title, column in (
        ("structure status (AlphaFold DB)", "structure_status"),
        ("key-annotation status (UniProt)", "key_annotation_status"),
    ):
        tally: Counter[tuple[str, str]] = Counter()
        for row in metrics:
            for prefix in ("5p", "3p"):
                tally[(prefix, row.get(f"{column}_{prefix}", ""))] += 1
        out.append("")
        out.append(
            table(
                ["parent", title, "records"],
                [[p, status, n] for (p, status), n in sorted(tally.items())],
            )
        )

    out.append("\n## Shortlist for PR 2 structure prediction\n")
    out.append(
        "One record per fusion+family, ordered by tier, then lowest SCHEMA E (max of the two "
        "parent structures), then GA margin, then read support. `ranked`: a row passes all six "
        "artefact flags; `fill`: a row passes all but the near-universal genomic_similarity and "
        "sample_burden flags (fired flags named); `named_case`: appended named fusions.\n"
    )
    out.append(table(CANDIDATE_HEADER, candidate_rows(short)))

    out.append("\n## Named cases (every record)\n")
    for name in NAMED:
        rows = [r for r in metrics if r["fusion"] == name]
        if not rows:
            out.append(f"\n**{name}**: no translated shared-domain record in the pilot.\n")
            continue
        out.append(f"\n**{name}** ({len(rows)} record(s))\n")
        out.append(table(CANDIDATE_HEADER, candidate_rows(rows)))
        for row in rows:
            out.append(
                f"\n- `{row['record_id']}` {row['junction_class']}: rules "
                f"`{row['tier_rules']}`; missing `{row['tier_missing_inputs'] or '-'}`; "
                f"disagreements `{row['metric_disagreements'] or '-'}`; aligned indel "
                f"+{row['aligned_inserted_residues']} ({row['aligned_inserted_sequence'] or '-'})"
                f"/-{row['aligned_deleted_columns']} cols; junction residues in HMM inserts "
                f"{row.get('junction_residues_in_inserts') or '-'}/"
                f"{row.get('junction_residues')}; delta bits 5'/3' "
                f"{row.get('delta_bits_vs_5p') or '-'}/{row.get('delta_bits_vs_3p') or '-'}; "
                f"HMMvar 5'/3' {row.get('hmmvar_delta_vs_5p') or '-'}/"
                f"{row.get('hmmvar_delta_vs_3p') or '-'}; junction DSSP 5' "
                f"`{row.get('junction_dssp_5p', '-')}` 3' `{row.get('junction_dssp_3p', '-')}`; "
                f"E_core 5'/3' {row.get('schema_E_core_5p')}/{row.get('schema_E_core_3p')}"
            )

    db = manifest["databases"]
    out.append("\n## Data coverage\n")
    out.append(
        f"- Records: {manifest['records']} from "
        f"{manifest['screen_rows_in_translated_classes']} translated-class screen rows; "
        f"rebuild problems: {manifest['rebuild_problems'] or 'none'}.\n"
        f"- Pfam {db.get('pfam_version')} ({db.get('pfam_release_date')}), "
        f"{db['pfam_families']} family HMMs, unavailable: "
        f"{len(db['pfam_families_unavailable'])}.\n"
        f"- UniProt {db['uniprot_release']}: {db['uniprot_genes']} genes, "
        f"{len(db['uniprot_genes_without_reviewed_entry'])} without a reviewed entry.\n"
        f"- AlphaFold DB model versions {db['afdb_model_versions']}: "
        f"{db['afdb_accessions']} accessions, {len(db['afdb_missing'])} without a usable "
        f"model, {len(db['afdb_problems'])} fetch/parse problems."
    )
    (run / "report_tables.md").write_text("\n".join(out) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
