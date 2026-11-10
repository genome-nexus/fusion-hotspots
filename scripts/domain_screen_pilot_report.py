"""Render the generated tables of the domain-screen pilot REPORT.md.

Usage::

    python -I scripts/domain_screen_pilot_report.py RUN_DIR > RUN_DIR/report_tables.md

``RUN_DIR`` holds ``flag_summary.json`` and ``<cohort>.domain_screen.tsv``
written by ``scripts/domain_screen_pilot_flags.py``. The output is Markdown:
class counts per cohort, every non-trivial call with its flags, flag reasons,
and the FLI1::ERG rows. Narrative sections of REPORT.md are written by hand.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path
from typing import Any

csv.field_size_limit(sys.maxsize)

CLASS_ORDER = (
    "clean_reconstitution",
    "ambiguous_junction_register",
    "tandem_duplication_or_insertion",
    "gap_or_deletion",
    "out_of_frame",
    "premature_stop_at_junction",
    "splice_site_unresolved",
    "not_split_on_both",
    "no_transcript_model",
    "no_shared_family",
)
NON_TRIVIAL = CLASS_ORDER[:7]
FLAGS = ("paralog_sequence_flag", "genomic_similarity_flag", "normal_recurrence_flag")
SHORT = {"flagged": "**FLAG**", "pass": "pass", "not_assessed": "n/a", "not_applicable": "-"}


def _cell(value: Any) -> str:
    text = "" if value is None else str(value)
    return text.replace("|", "\\|").replace("\n", " ") or "-"


def _table(header: list[str], rows: list[list[Any]]) -> str:
    lines = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    lines += ["| " + " | ".join(_cell(value) for value in row) + " |" for row in rows]
    return "\n".join(lines)


def _rows(run_dir: Path, cohort: str) -> list[dict[str, str]]:
    with (run_dir / f"{cohort}.domain_screen.tsv").open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _frame(row: dict[str, str]) -> str:
    return f"{row['translated_frame_status']} / alt {row['alternative_convention_frame_status']}"


def _reads(row: dict[str, str]) -> str:
    split, pairs = row.get("tumor_split_reads") or "?", row.get("tumor_paired_end_reads") or "?"
    text = f"{split} split / {pairs} pairs"
    if row.get("breakpoint_precision"):
        text += f" ({row['breakpoint_precision']})"
    if row.get("Open_ffpm"):
        text += f", FFPM {row['Open_ffpm']}"
    if row.get("Open_fi_validated"):
        text += f", FI {row['Open_fi_validated']}"
    return text


def _sample(row: dict[str, str]) -> str:
    sample = row.get("sample_id") or "-"
    detail = [row.get("Open_tissue_class"), row.get("Open_tumor_type")]
    detail_text = ", ".join(part for part in detail if part)
    return f"{sample} ({detail_text})" if detail_text else sample


def render(run_dir: Path) -> str:
    summary = json.loads((run_dir / "flag_summary.json").read_text())
    screens: dict[str, Any] = summary["screen_summaries"]
    cohorts = list(screens)
    out: list[str] = []

    out.append("### Class counts per cohort (family rows)\n")
    header = ["Class", *cohorts]
    body = [
        [name, *(f"{screens[c]['class_counts'].get(name, 0):,}" for c in cohorts)]
        for name in CLASS_ORDER
    ]
    body.append(
        [
            "**screened fusion events**",
            *(f"{screens[c]['screened_fusion_events']:,}" for c in cohorts),
        ]
    )
    body.append(
        [
            "events sharing a family",
            *(f"{screens[c]['events_with_a_shared_family']:,}" for c in cohorts),
        ]
    )
    body.append(["shared-family rows", *(f"{screens[c]['family_rows']:,}" for c in cohorts)])
    out.append(_table(header, body))

    out.append("\n### Flag counts per cohort (shared-family rows)\n")
    sub = (
        "paralog_sequence_flag",
        "genomic_similarity_flag",
        "breakpoint_repeat_flag",
        "sample_burden_flag",
        "read_evidence_flag",
        "normal_recurrence_flag",
    )
    flag_rows = []
    for cohort in cohorts:
        counts = summary["flag_counts"][cohort]
        flag_rows.append(
            [
                cohort,
                *(
                    ", ".join(f"{k} {v}" for k, v in sorted(counts.get(flag, {}).items()))
                    for flag in sub
                ),
            ]
        )
    out.append(_table(["Cohort", *sub], flag_rows))

    all_rows = {cohort: _rows(run_dir, cohort) for cohort in cohorts}
    out.append("\n### Every non-trivial call\n")
    out.append(
        "Frame: `translated` is the screen's convention (`coordinate_convention`); `alt` is the "
        "other convention. Flags: **FLAG** = flagged, pass, n/a = not assessed (reason below).\n"
    )
    detail: list[str] = []
    for name in NON_TRIVIAL:
        hits = [
            (cohort, row)
            for cohort in cohorts
            for row in all_rows[cohort]
            if row["reconstitution_class"] == name
        ]
        out.append(f"\n#### {name} ({len(hits)})\n")
        if not hits:
            out.append("None.\n")
            continue
        table = []
        for cohort, row in hits:
            table.append(
                [
                    f"{row['five_prime_gene']}::{row['three_prime_gene']}",
                    cohort,
                    _sample(row),
                    row["pfam_family"],
                    f"`{row['junction_peptide']}`" if row.get("junction_peptide") else "-",
                    row.get("junction_indel_sequence") or "-",
                    _reads(row),
                    _frame(row),
                    *(SHORT.get(row[flag], row[flag]) for flag in FLAGS),
                ]
            )
            detail.append(
                f"- **{row['five_prime_gene']}::{row['three_prime_gene']}** {cohort} "
                f"{row.get('sample_id')} {row['pfam_family']} ({name}); breakpoints "
                f"{row.get('five_prime_chromosome')}:{row.get('five_prime_breakpoint')} / "
                f"{row.get('three_prime_chromosome')}:{row.get('three_prime_breakpoint')} "
                f"{row['genome_build']}\n"
                f"  - paralog: {row['paralog_sequence_flag']}: {row['paralog_sequence_reason']}\n"
                f"  - genomic: {row['genomic_similarity_flag']}: "
                f"{row['genomic_similarity_reason']} "
                f"(identity/length +/-150 {row.get('genomic_identity_150') or '-'}/"
                f"{row.get('genomic_aligned_length_150') or '-'}, +/-500 "
                f"{row.get('genomic_identity_500') or '-'}/"
                f"{row.get('genomic_aligned_length_500') or '-'}; repeat "
                f"{row['breakpoint_repeat_flag']}: {row['breakpoint_repeat_reason']})\n"
                f"  - normal: {row['normal_recurrence_flag']}: {row['normal_recurrence_reason']}\n"
            )
        out.append(
            _table(
                [
                    "Fusion",
                    "Cohort",
                    "Sample",
                    "Family",
                    "Junction peptide",
                    "Indel",
                    "Reads",
                    "Frame (translated / alt)",
                    "Paralog",
                    "Genomic",
                    "Normal",
                ],
                table,
            )
        )
    out.append("\n### Flag reasons for every non-trivial call\n")
    out.append("".join(detail) or "None.\n")

    out.append("\n### FLI1::ERG rows (all cohorts)\n")
    fli = [
        [
            cohort,
            _sample(row),
            row["pfam_family"],
            row["reconstitution_class"],
            row.get("junction_mechanism"),
            f"`{row['junction_peptide']}`" if row.get("junction_peptide") else "-",
            _reads(row),
            _frame(row),
            f"{row.get('five_prime_break_location')} | {row.get('three_prime_break_location')}",
            *(SHORT.get(row[flag], row[flag]) for flag in FLAGS),
        ]
        for cohort in cohorts
        for row in all_rows[cohort]
        if (row["five_prime_gene"], row["three_prime_gene"]) == ("FLI1", "ERG")
    ]
    out.append(
        _table(
            [
                "Cohort",
                "Sample",
                "Family",
                "Class",
                "Mechanism",
                "Junction peptide",
                "Reads",
                "Frame",
                "Break locations (5' | 3')",
                "Paralog",
                "Genomic",
                "Normal",
            ],
            fli,
        )
        if fli
        else "No FLI1::ERG row in any cohort."
    )
    return "\n".join(out) + "\n"


if __name__ == "__main__":
    sys.stdout.write(render(Path(sys.argv[1])))
