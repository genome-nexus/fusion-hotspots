"""Render the generated tables of the domain-screen pilot REPORT.md.

Usage::

    python -I scripts/domain_screen_pilot_report.py RUN_DIR > RUN_DIR/report_tables.md

``RUN_DIR`` holds ``flag_summary.json`` and ``<cohort>.domain_screen.tsv`` (or
``.tsv.gz``) written by ``scripts/domain_screen_pilot_flags.py``. Output is
Markdown: class and flag counts per cohort, the candidates that pass all three
flags, every non-trivial call grouped by distinct junction (each sample listed
with its read support), and the FLI1::ERG rows. Narrative sections of
REPORT.md are written by hand.
"""

from __future__ import annotations

import csv
import gzip
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

from cfh.cohort.artefact_flags import GENOMIC_MIN_IDENTITY, GENOMIC_MIN_LENGTH

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
HIGHLIGHT = CLASS_ORDER[:4]
FLAGS = ("paralog_sequence_flag", "genomic_similarity_flag", "normal_recurrence_flag")
SUB_FLAGS = (
    "paralog_sequence_flag",
    "genomic_similarity_flag",
    "breakpoint_repeat_flag",
    "sample_burden_flag",
    "read_evidence_flag",
    "normal_recurrence_flag",
)
SHORT = {"flagged": "**FLAG**", "pass": "pass", "not_assessed": "n/a", "not_applicable": "-"}


def _cell(value: Any) -> str:
    text = "" if value is None else str(value)
    return text.replace("|", "\\|").replace("\n", " ") or "-"


def _table(header: list[str], rows: list[list[Any]]) -> str:
    lines = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    lines += ["| " + " | ".join(_cell(value) for value in row) + " |" for row in rows]
    return "\n".join(lines)


def _rows(run_dir: Path, cohort: str) -> list[dict[str, str]]:
    plain = run_dir / f"{cohort}.domain_screen.tsv"
    if plain.exists():
        with plain.open(newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))
    with gzip.open(run_dir / f"{cohort}.domain_screen.tsv.gz", "rt", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _frame(row: dict[str, str]) -> str:
    return f"{row['translated_frame_status']} / alt {row['alternative_convention_frame_status']}"


def _reads(row: dict[str, str]) -> str:
    split, pairs = row.get("tumor_split_reads") or "?", row.get("tumor_paired_end_reads") or "?"
    text = f"{split}/{pairs}"
    if row.get("breakpoint_precision"):
        text += f" {row['breakpoint_precision']}"
    if row.get("Open_ffpm"):
        text += f" FFPM {row['Open_ffpm']}"
    if row.get("Open_fi_validated"):
        text += f" FI {row['Open_fi_validated']}"
    return text


def _sample(row: dict[str, str]) -> str:
    detail = [row.get("Open_tissue_class"), row.get("Open_tumor_type")]
    detail_text = ",".join(part for part in detail if part)
    sample = row.get("sample_id") or "-"
    return f"{sample} ({detail_text})" if detail_text else sample


def _similar(row: dict[str, str], window: int) -> bool:
    try:
        identity = float(row.get(f"genomic_identity_{window}") or 0)
        length = int(float(row.get(f"genomic_aligned_length_{window}") or 0))
    except ValueError:
        return False
    return identity >= GENOMIC_MIN_IDENTITY and length >= GENOMIC_MIN_LENGTH


def _genomic_parts(row: dict[str, str]) -> str:
    """Which genomic sub-checks fired: sim150/sim500, rep, burdenN, reads."""
    parts = [f"sim{w}" for w in (150, 500) if _similar(row, w)]
    if row.get("breakpoint_repeat_flag") == "flagged":
        parts.append("rep")
    elif row.get("breakpoint_repeat_flag") == "not_assessed":
        parts.append("rep n/a")
    if row.get("sample_burden_flag") == "flagged":
        parts.append(f"burden{row.get('sample_shared_family_calls')}")
    if row.get("read_evidence_flag") == "flagged":
        parts.append("reads")
    return ", ".join(parts) or "pass"


def _junction_key(cohort: str, row: dict[str, str]) -> tuple[str, ...]:
    return (
        cohort,
        row["five_prime_gene"],
        row["three_prime_gene"],
        row.get("five_prime_chromosome") or "",
        row.get("five_prime_breakpoint") or "",
        row.get("three_prime_chromosome") or "",
        row.get("three_prime_breakpoint") or "",
        row["pfam_family"],
        row["reconstitution_class"],
    )


def _passes_all(row: dict[str, str]) -> bool:
    return all(row[flag] == "pass" for flag in FLAGS)


def _detail(cohort: str, row: dict[str, str]) -> str:
    return (
        f"- **{row['five_prime_gene']}::{row['three_prime_gene']}** {cohort} "
        f"{row.get('sample_id')} {row['pfam_family']} ({row['reconstitution_class']}); "
        f"{row.get('five_prime_chromosome')}:{row.get('five_prime_breakpoint')} / "
        f"{row.get('three_prime_chromosome')}:{row.get('three_prime_breakpoint')} "
        f"{row['genome_build']}; transcripts {row.get('five_prime_transcript')} / "
        f"{row.get('three_prime_transcript')}; junction `{row.get('junction_peptide') or '-'}`; "
        f"reads {_reads(row)}; frame {_frame(row)}\n"
        f"  - paralog: {row['paralog_sequence_flag']}: {row['paralog_sequence_reason']}\n"
        f"  - genomic: {row['genomic_similarity_flag']}: "
        f"{row['genomic_similarity_reason']} (identity/aligned bp: +/-150 "
        f"{row.get('genomic_identity_150') or '-'}/{row.get('genomic_aligned_length_150') or '-'}"
        f", +/-500 {row.get('genomic_identity_500') or '-'}/"
        f"{row.get('genomic_aligned_length_500') or '-'})\n"
        f"  - normal: {row['normal_recurrence_flag']}: {row['normal_recurrence_reason']}\n"
    )


def _counts_tables(summary: dict[str, Any], cohorts: list[str]) -> list[str]:
    screens = summary["screen_summaries"]
    body = [
        [name, *(f"{screens[c]['class_counts'].get(name, 0):,}" for c in cohorts)]
        for name in CLASS_ORDER
    ]
    for label, key in (
        ("screened fusion events", "screened_fusion_events"),
        ("events sharing a family", "events_with_a_shared_family"),
        ("shared-family rows", "family_rows"),
    ):
        body.append([f"*{label}*", *(f"{screens[c][key]:,}" for c in cohorts)])
    flag_rows = [
        [
            cohort,
            *(
                ", ".join(
                    f"{k} {v:,}"
                    for k, v in sorted(summary["flag_counts"][cohort].get(flag, {}).items())
                )
                for flag in SUB_FLAGS
            ),
        ]
        for cohort in cohorts
    ]
    return [
        "### Class counts per cohort (family rows)\n",
        _table(["Class", *cohorts], body),
        "\n### Flag counts per cohort (shared-family rows)\n",
        "`no_transcript_model` rows never reach junction translation and are "
        "`not_applicable`; `not_assessed` always carries a reason in the TSV.\n",
        _table(["Cohort", *SUB_FLAGS], flag_rows),
    ]


def _survivors(cohorts: list[str], all_rows: dict[str, list[dict[str, str]]]) -> list[str]:
    out = [
        "\n### Clean / ambiguous / indel calls that pass all three flags\n\n"
        "Rows of class clean_reconstitution, ambiguous_junction_register, "
        "tandem_duplication_or_insertion or gap_or_deletion whose paralog, genomic and normal "
        "flags are all `pass`. Full reasons follow the table.\n"
    ]
    survivors = [
        (cohort, row)
        for cohort in cohorts
        for row in all_rows[cohort]
        if row["reconstitution_class"] in HIGHLIGHT and _passes_all(row)
    ]
    if not survivors:
        return [*out, "None.\n"]
    table = [
        [
            f"{row['five_prime_gene']}::{row['three_prime_gene']}",
            cohort,
            _sample(row),
            row["pfam_family"],
            row["reconstitution_class"],
            f"`{row['junction_peptide']}`" if row.get("junction_peptide") else "-",
            _reads(row),
            _frame(row),
        ]
        for cohort, row in survivors
    ]
    header = ["Fusion", "Cohort", "Sample", "Family", "Class", "Junction peptide", "Reads", "Frame"]
    return [*out, _table(header, table), "\n" + "".join(_detail(c, r) for c, r in survivors)]


def _by_junction(cohorts: list[str], all_rows: dict[str, list[dict[str, str]]]) -> list[str]:
    out = [
        "\n### Every non-trivial call, by distinct junction\n",
        "One line per cohort x gene pair x breakpoints x Pfam family x class; every sample "
        "with that junction is listed with its reads (`split/pairs`, MSK precision, open-table "
        "FFPM and FusionInspector validation). Frame: `translated` is the screen's "
        "both-bases-retained convention, `alt` drops the reported 3' base. Genomic column: "
        f"sub-checks that fired (`sim150`/`sim500` = breakpoint windows >="
        f"{GENOMIC_MIN_IDENTITY:.0%} identical over >={GENOMIC_MIN_LENGTH} bp, `rep` = segdup/"
        "RepeatMasker, `burdenN` = sample has N shared-family calls, `reads` = 0 split reads or "
        "imprecise); `burden` and `reads` are per sample, so the column shows the union. "
        "Per-row reasons are in the TSVs.\n",
    ]
    header = [
        "Fusion",
        "Cohort",
        "Family",
        "Breakpoints",
        "Junction peptide",
        "Indel",
        "Frame (translated / alt)",
        "Paralog",
        "Genomic",
        "Normal",
        "Samples: reads",
    ]
    for name in NON_TRIVIAL:
        groups: dict[tuple[str, ...], list[dict[str, str]]] = defaultdict(list)
        for cohort in cohorts:
            for row in all_rows[cohort]:
                if row["reconstitution_class"] == name:
                    groups[_junction_key(cohort, row)].append(row)
        total = sum(len(rows) for rows in groups.values())
        out.append(f"\n#### {name} ({total:,} rows, {len(groups):,} distinct junctions)\n")
        if not groups:
            out.append("None.\n")
            continue
        table = []
        for key, rows in sorted(groups.items(), key=lambda item: (item[0][0], -len(item[1]))):
            first = rows[0]
            genomic = sorted({part for row in rows for part in _genomic_parts(row).split(", ")})
            if len(genomic) > 1 and "pass" in genomic:
                genomic.remove("pass")
            samples = "; ".join(f"{_sample(row)} {_reads(row)}" for row in rows)
            table.append(
                [
                    f"{key[1]}::{key[2]}",
                    key[0],
                    key[7],
                    f"{key[3]}:{key[4]} / {key[5]}:{key[6]}",
                    f"`{first['junction_peptide']}`" if first.get("junction_peptide") else "-",
                    first.get("junction_indel_sequence") or "-",
                    _frame(first),
                    "/".join(sorted({SHORT.get(r["paralog_sequence_flag"], "?") for r in rows})),
                    ", ".join(genomic),
                    "/".join(sorted({SHORT.get(r["normal_recurrence_flag"], "?") for r in rows})),
                    f"{len(rows)}: {samples}",
                ]
            )
        out.append(_table(header, table))
    return out


def _fli1_erg(cohorts: list[str], all_rows: dict[str, list[dict[str, str]]]) -> list[str]:
    out = ["\n### FLI1::ERG rows (all cohorts)\n"]
    fli = [
        (cohort, row)
        for cohort in cohorts
        for row in all_rows[cohort]
        if (row["five_prime_gene"], row["three_prime_gene"]) == ("FLI1", "ERG")
    ]
    if not fli:
        return [*out, "No FLI1::ERG row in any cohort.\n"]
    header = [
        "Cohort",
        "Sample",
        "Family",
        "Class",
        "Mechanism",
        "Break locations (5' | 3')",
        "Junction peptide",
        "Reads",
        "Frame",
        "Paralog",
        "Genomic",
        "Normal",
    ]
    table = [
        [
            cohort,
            _sample(row),
            row["pfam_family"],
            row["reconstitution_class"],
            row.get("junction_mechanism"),
            f"{row.get('five_prime_break_location')} | {row.get('three_prime_break_location')}",
            f"`{row['junction_peptide']}`" if row.get("junction_peptide") else "-",
            _reads(row),
            _frame(row),
            SHORT.get(row["paralog_sequence_flag"]),
            _genomic_parts(row),
            SHORT.get(row["normal_recurrence_flag"]),
        ]
        for cohort, row in fli
    ]
    return [*out, _table(header, table), "\n" + "".join(_detail(c, r) for c, r in fli)]


def render(run_dir: Path) -> str:
    summary = json.loads((run_dir / "flag_summary.json").read_text())
    cohorts = list(summary["screen_summaries"])
    all_rows = {cohort: _rows(run_dir, cohort) for cohort in cohorts}
    out = [
        *_counts_tables(summary, cohorts),
        *_survivors(cohorts, all_rows),
        *_fli1_erg(cohorts, all_rows),
    ]
    msk = [
        ("msk_impact_50k_2026", row)
        for row in all_rows.get("msk_impact_50k_2026", [])
        if row["reconstitution_class"] in NON_TRIVIAL
    ]
    if msk:
        out.append("\n### MSK-IMPACT 50k non-trivial calls with full flag reasons\n")
        out.append("".join(_detail(cohort, row) for cohort, row in msk))
    out.extend(_by_junction(cohorts, all_rows))
    return "\n".join(out) + "\n"


if __name__ == "__main__":
    sys.stdout.write(render(Path(sys.argv[1])))
