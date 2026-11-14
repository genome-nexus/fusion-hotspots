"""Committed outputs: metrics TSV, junction profiles, counts, report tables and manifest.

* ``plm_metrics.tsv.gz`` -- one row per record (every PR #155 record plus the
  PCAWG/TCGA DNA records), PR #155 tier and PR #160 flag copied unchanged next
  to the PLM columns;
* ``junction_profiles.tsv.gz`` -- per-residue masked log-probabilities over the
  chimera's junction window (scored records only);
* ``plm_counts.json`` -- coverage and cross-tabulations;
* ``report_tables.md`` -- the tables REPORT.md quotes;
* ``run_manifest.json`` -- inputs (sha256), model spec, method constants.
"""

from __future__ import annotations

import csv
import gzip
import io
import json
from collections import Counter
from pathlib import Path
from typing import Any, Iterable

from cfh.domain_plm.compare import DISAGREEMENT_RULES, RANK_MEASURES
from cfh.domain_plm.metrics import METRIC_COLUMNS

ID_COLUMNS = (
    "record_id",
    "sequence_id",
    "fusion",
    "source",
    "scoring_group",
    "pfam_family",
    "junction_class",
    "n_samples",
    "chimeric_domain_span",
    "chimeric_domain_length",
    "junction_index",
    "pilot_indel_residues",
    "pr155_tier",
    "pr160_fold_flag",
)
COLUMNS = (
    ID_COLUMNS
    + METRIC_COLUMNS
    + tuple(f"rank_{m}" for m in RANK_MEASURES)
    + ("plm_rank", "disagreements")
)
PROFILE_COLUMNS = ("record_id", "fusion", "chimera_position", "residue", "logp")


def fmt(value: Any) -> str:
    if isinstance(value, float):
        return f"{value:.4f}"
    return "" if value is None else str(value)


def _write_tsv_gz(path: Path, columns: Iterable[str], rows: Iterable[Iterable[Any]]) -> None:
    buffer = io.StringIO()
    writer = csv.writer(buffer, delimiter="\t", lineterminator="\n")
    writer.writerow(list(columns))
    for row in rows:
        writer.writerow([fmt(v) for v in row])
    # mtime=0: byte-identical output for identical content
    path.write_bytes(gzip.compress(buffer.getvalue().encode(), mtime=0))


def counts(rows: list[dict[str, Any]]) -> dict[str, Any]:
    pr155 = [r for r in rows if r["source"] == "pr155"]
    scored = [r for r in pr155 if r["plm_status"] == "scored"]

    def tab(items: list[dict[str, Any]], *keys: str) -> dict[str, int]:
        return dict(sorted(Counter("|".join(str(r[k]) for k in keys) for r in items).items()))

    return {
        "records": len(rows),
        "pr155_records": len(pr155),
        "pcawg_tcga_records": len(rows) - len(pr155),
        "plm_status": tab(rows, "plm_status"),
        "pr155_plm_status_by_group": tab(pr155, "scoring_group", "plm_status"),
        "pr155_plm_status_by_tier": tab(pr155, "pr155_tier", "plm_status"),
        "pr155_plm_status_by_junction_class": tab(pr155, "junction_class", "plm_status"),
        "not_assessed_reasons": dict(
            Counter(
                r["plm_not_assessed_reason"].split(" (")[0]
                for r in rows
                if r["plm_status"] != "scored"
            ).most_common()
        ),
        "scored_pr155_tier_x_plm_rank": tab(scored, "pr155_tier", "plm_rank"),
        "scored_pr155_junction_class_x_plm_rank": tab(scored, "junction_class", "plm_rank"),
        "fold_flag_x_plm_rank": tab(
            [r for r in rows if r["pr160_fold_flag"] != "not_in_structure_run"],
            "pr160_fold_flag",
            "plm_rank",
        ),
        "disagreements": dict(
            sorted(Counter(d for r in rows for d in r["disagreements"].split(";") if d).items())
        ),
    }


def _cell(row: dict[str, Any], key: str) -> str:
    value = row.get(key, "")
    return fmt(value) if value != "" else "–"


TABLE_COLUMNS = (
    ("fusion", "fusion"),
    ("Pfam", "pfam_family"),
    ("class", "junction_class"),
    ("indel", "pilot_indel_residues"),
    ("PR #155 tier", "pr155_tier"),
    ("PR #160 fold", "pr160_fold_flag"),
    ("PLL chim", "chimera_domain_pll"),
    ("PLL 5'", "parent5_domain_pll"),
    ("PLL 3'", "parent3_domain_pll"),
    ("Δ vs 5'", "delta_vs_5p"),
    ("Δ vs 3'", "delta_vs_3p"),
    ("ctx 5'", "context_shift_5p"),
    ("ctx 3'", "context_shift_3p"),
    ("jx dip", "chimera_junction_dip"),
    ("dip excess (worst)", "dip_excess_worst"),
    ("PLM rank", "plm_rank"),
    ("disagreements", "disagreements"),
)


def table(rows: list[dict[str, Any]], extra: tuple[tuple[str, str], ...] = ()) -> str:
    columns = extra + TABLE_COLUMNS
    lines = [
        "| " + " | ".join(h for h, _ in columns) + " |",
        "|" + "---|" * len(columns),
    ]
    for row in rows:
        cells = []
        for _, key in columns:
            cell = _cell(row, key)
            if key == "plm_rank" and row["plm_status"] != "scored":
                cell = f"not assessed: {row['plm_not_assessed_reason']}"
            cells.append(cell.replace("|", "/").replace(";", "; "))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def report_tables(rows: list[dict[str, Any]], named: tuple[str, ...]) -> str:
    shortlist = [r for r in rows if r["scoring_group"] == "shortlist"]
    named_rows = [r for r in rows if r["fusion"] in named or r["source"] == "pcawg_tcga_dna"]
    named_rows.sort(
        key=lambda r: (
            named.index(r["fusion"]) if r["fusion"] in named else len(named),
            r["fusion"],
            r["record_id"],
        )
    )
    disagreeing = [r for r in rows if r["disagreements"]]
    disagreeing.sort(key=lambda r: (r["scoring_group"] != "shortlist", r["fusion"]))
    parts = [
        "## Shortlist (PR #155 shortlist, PR #160 fold flag)\n",
        table(shortlist, (("record", "record_id"),)),
        "\n\n## Named cases and every PCAWG DNA event\n",
        table(named_rows, (("record", "record_id"), ("source", "source"))),
        "\n\n## Disagreements\n",
        "\n".join(f"* `{key}` -- {text}" for key, text in DISAGREEMENT_RULES.items()),
        "\n\n",
        table(disagreeing, (("record", "record_id"), ("group", "scoring_group")))
        if disagreeing
        else "_none_",
        "\n",
    ]
    return "".join(parts)


def write_outputs(
    out: Path,
    rows: list[dict[str, Any]],
    profiles: list[tuple[str, str, int, str, float]],
    manifest: dict[str, Any],
    named: tuple[str, ...],
) -> None:
    out.mkdir(parents=True, exist_ok=True)
    _write_tsv_gz(
        out / "plm_metrics.tsv.gz", COLUMNS, ([row.get(c, "") for c in COLUMNS] for row in rows)
    )
    _write_tsv_gz(out / "junction_profiles.tsv.gz", PROFILE_COLUMNS, profiles)
    (out / "plm_counts.json").write_text(json.dumps(counts(rows), indent=1) + "\n")
    (out / "report_tables.md").write_text(report_tables(rows, named))
    (out / "run_manifest.json").write_text(json.dumps(manifest, indent=1, default=str) + "\n")
