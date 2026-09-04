"""Write ``cfh cohort-scan`` output artifacts to ``runs/<run_id>/cohort_scan/``.

Two kinds of artifact are produced:

* One consolidated summary (TSV + JSON + Markdown + PDF) covering every
  scanned gene -- recurrence, in-frame%, domain-retention%, FDR-adjusted
  p-values, and composite score -- sorted by significance. The summary also
  calls out an "honorable mentions" tier (see :func:`build_honorable_mentions`):
  the top-ranked genes by raw Fisher p-value among genes that did NOT survive
  genome-wide FDR correction, so a human reviewer can see the real
  second-tier signal that a strict significant/not-significant split hides.
* Full detailed per-gene reports (the same ``results.tsv``/``results.json``/
  ``report.md``/``report.pdf``/``manifest.json`` that ``cfh analyze`` already
  produces, via :func:`cfh.real_benchmark.write_outputs`) for every
  hand-curated gene, every FDR-significant (q<0.05) gene, and every
  honorable-mention gene -- never for every scanned gene, so a run with
  hundreds of non-significant genes does not generate hundreds of PDFs.
"""

from __future__ import annotations

import csv
import json
import math
from datetime import datetime, timezone
from pathlib import Path

from cfh.cohort.scan import (
    CohortScanResult,
    genes_needing_full_report,
    per_gene_min_q_value,
)
from cfh.real_benchmark import write_outputs
from cfh.reporting.manhattan import render_manhattan_svg
from cfh.reporting.manuscript_text import (
    render_discussion_bullets,
    render_gene_highlight,
    render_manhattan_caption,
    render_manuscript_abstract,
    render_manuscript_methods,
    render_manuscript_title,
)
from cfh.reporting.pdf import render_cohort_summary_pdf, render_manuscript_pdf

_MANHATTAN_SVG_FILENAME = "manhattan.svg"
_MANUSCRIPT_MARKDOWN_FILENAME = "paper.md"
_MANUSCRIPT_PDF_FILENAME = "paper.pdf"

DEFAULT_HONORABLE_MENTION_COUNT = 15
"""Default size of the "honorable mentions" / highly ranked non-FDR-significant
tier: the top N genes by raw Fisher p-value among genes that did NOT survive
genome-wide FDR correction. Configurable via ``write_cohort_scan_outputs``'s
``honorable_mention_count`` (and the ``cfh cohort-scan
--honorable-mention-count`` CLI flag)."""

_HONORABLE_MENTION_HEADING = (
    "Honorable mentions: highly ranked non-FDR-significant genes worth human review"
)

_HONORABLE_MENTION_NOTE = (
    "Did not survive genome-wide multiple-testing correction (FDR-adjusted "
    "q-value at or above the significance threshold), but ranks highly by "
    "raw p-value among the non-FDR-significant genes and may warrant "
    "targeted follow-up. This is NOT a claim of statistical significance."
)

_HONORABLE_MENTION_NOTE_NO_Q = (
    "No FDR-adjusted q-value was ever computed for this gene (it contributed no "
    "p-value to the genome-wide Benjamini-Hochberg correction), so whether it "
    "would have survived multiple-testing correction is unknown -- this is NOT "
    "a claim about where its (nonexistent) q-value would fall relative to the "
    "significance threshold. It ranks highly by raw p-value among the "
    "non-FDR-significant genes and may warrant targeted follow-up. This is NOT "
    "a claim of statistical significance."
)
"""Used in place of :data:`_HONORABLE_MENTION_NOTE` for a candidate gene that
has a real raw Fisher p-value (so it is eligible for this ranked tier) but
whose ``min_fdr_adjusted_q_value`` is ``None`` -- i.e. it never actually
contributed a hypothesis to the BH correction. Asserting the standard note's
"FDR-adjusted q-value at or above the significance threshold" framing for
such a gene would invent a verdict about a q-value that was never derived."""

_SUMMARY_FIELDNAMES = [
    "gene_symbol",
    "config_source",
    "status",
    "distinct_patient_count",
    "total_sv_count",
    "n_events_analyzed",
    "in_frame_percent",
    "domain_retention_percent",
    "fisher_p_value",
    "permutation_p_value",
    "min_fdr_adjusted_q_value",
    "fdr_significant",
    "top_composite_score",
    "top_composite_partner_gene",
    "error",
]


def _top_composite(run) -> tuple[float | None, str | None]:
    if run is None:
        return None, None
    for result in run.results:
        if result.Algorithm != "composite_score":
            continue
        ranking = (result.Tables or {}).get("composite_evidence_ranking") or []
        scored = [row for row in ranking if row.get("Composite_score") is not None]
        if not scored:
            return None, None
        top = max(scored, key=lambda row: row["Composite_score"])
        return top["Composite_score"], top.get("Partner_gene")
    return None, None


def _algorithm_summary(run, algorithm: str) -> dict:
    if run is None:
        return {}
    return next(
        (result.Summary for result in run.results if result.Algorithm == algorithm),
        {},
    )


def build_summary_rows(result: CohortScanResult) -> list[dict]:
    """Build one row per scanned gene, sorted by significance (most
    significant FDR-adjusted q-value first; genes with no q-value sort
    last, tie-broken by descending recurrence).
    """
    min_q_by_gene = per_gene_min_q_value(result)
    significant = set(result.significant_genes)
    rows: list[dict] = []
    for outcome in result.gene_outcomes:
        summary = outcome.run.summary if outcome.run is not None else {}
        retention = _algorithm_summary(outcome.run, "domain_retention")
        disruption = _algorithm_summary(outcome.run, "domain_disruption")
        top_score, top_partner = _top_composite(outcome.run)
        q_value = min_q_by_gene.get(outcome.gene_symbol)
        rows.append(
            {
                "gene_symbol": outcome.gene_symbol,
                "config_source": outcome.config_source,
                "status": outcome.status,
                "distinct_patient_count": outcome.distinct_patient_count,
                "total_sv_count": outcome.total_sv_count,
                "n_events_analyzed": summary.get("total_fusions"),
                "in_frame_percent": summary.get("in_frame_percent"),
                "domain_retention_percent": summary.get("kinase_retained_percent"),
                "fisher_p_value": summary.get("fisher_p_value"),
                "fisher_odds_ratio": retention.get("fisher_odds_ratio"),
                "key_domains": summary.get("configured_key_domains") or [],
                "disruption_fisher_p_value": disruption.get("fisher_p_value"),
                "disruption_fisher_odds_ratio": disruption.get("fisher_odds_ratio"),
                "disruption_required_domains": summary.get("configured_disruption_required_domains")
                or [],
                "permutation_p_value": summary.get("permutation_p_value"),
                "min_fdr_adjusted_q_value": q_value,
                # ``None`` (never computed -- this gene contributed no
                # p-value to the BH correction) is kept distinct from
                # ``False`` (a real q-value was computed and it was >= the
                # significance threshold); collapsing the two into a single
                # ``False`` would let downstream text claim a specific
                # "did not reach FDR significance" verdict for a gene that
                # was never actually tested.
                "fdr_significant": None if q_value is None else outcome.gene_symbol in significant,
                "top_composite_score": top_score,
                "top_composite_partner_gene": top_partner,
                "error": outcome.error,
            }
        )

    def _sort_key(row: dict) -> tuple:
        q_value = row["min_fdr_adjusted_q_value"]
        return (
            q_value is None,
            q_value if q_value is not None else 1.0,
            -(row["distinct_patient_count"] or 0),
            row["gene_symbol"],
        )

    rows.sort(key=_sort_key)
    return rows


def build_honorable_mentions(
    rows: list[dict], *, limit: int = DEFAULT_HONORABLE_MENTION_COUNT
) -> list[dict]:
    """The "honorable mentions" / highly ranked non-FDR-significant tier: the
    top ``limit`` genes by raw Fisher p-value among genes that are NOT
    FDR-significant.

    This is a distinct, additive ranking -- never a claim that any of these
    genes IS significant, or even that it is close to significant (a run's
    544 genes can rank from a raw p-value of 0.0004 down to 0.9, and this
    tier surfaces the top of that non-significant ranking regardless of how
    far any individual gene actually sits from the threshold). A run can
    genuinely have only one FDR-significant gene (a strict, correct
    genome-wide multiple-testing result) while still having a real,
    biologically sensible second tier of genes that rank highly on raw
    p-value but did not survive correction; without this, that tier is
    invisible to a human reviewer skimming a binary significant/not-significant
    summary.
    """
    candidates = [
        row for row in rows if not row["fdr_significant"] and row.get("fisher_p_value") is not None
    ]
    candidates.sort(key=lambda row: (row["fisher_p_value"], row["gene_symbol"]))
    mentions = []
    for rank, row in enumerate(candidates[:limit], start=1):
        q_value = row["min_fdr_adjusted_q_value"]
        mentions.append(
            {
                "rank": rank,
                "gene_symbol": row["gene_symbol"],
                "fisher_p_value": row["fisher_p_value"],
                "min_fdr_adjusted_q_value": q_value,
                "n_events_analyzed": row["n_events_analyzed"],
                "in_frame_percent": row["in_frame_percent"],
                "domain_retention_percent": row["domain_retention_percent"],
                # A candidate's raw Fisher p-value alone makes it eligible
                # for this tier; it does not guarantee a q-value was ever
                # computed for it (see ``_HONORABLE_MENTION_NOTE_NO_Q``).
                "note": (
                    _HONORABLE_MENTION_NOTE if q_value is not None else _HONORABLE_MENTION_NOTE_NO_Q
                ),
            }
        )
    return mentions


def _json_safe(value):
    if isinstance(value, dict):
        return {key: _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return str(value)
    return value


def _format_cell(value) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.4g}"
    if isinstance(value, bool):
        return "yes" if value else "no"
    return str(value)


_HONORABLE_MENTION_FIELDNAMES = [
    "rank",
    "gene_symbol",
    "fisher_p_value",
    "min_fdr_adjusted_q_value",
    "n_events_analyzed",
    "in_frame_percent",
    "domain_retention_percent",
]


def _write_summary_markdown(
    result: CohortScanResult,
    rows: list[dict],
    honorable_mentions: list[dict],
    path: Path,
) -> None:
    lines = [
        f"# Genome-wide fusion-hotspot cohort scan: {result.study_id}",
        "",
        f"- Total genes with any structural-variant record in the cohort: "
        f"{result.total_genes_before_gating}",
        f"- Genes passing the >= {result.min_distinct_patients}-distinct-patient "
        f"recurrence gate: {result.genes_after_gating}",
        f"- Curated gene configs used: {result.curated_gene_count}",
        f"- Auto-generated gene configs used: {result.auto_config_gene_count}",
        f"- Genes gated in but unresolvable (no Genome Nexus canonical transcript): "
        f"{result.unresolved_gene_count}",
        f"- FDR-significant genes (q < 0.05) after Benjamini-Hochberg correction across "
        f"all {len(rows)} scanned genes: {len(result.significant_genes)}",
        "",
    ]
    if result.warnings:
        lines.extend(["## Warnings", ""])
        lines.extend(f"- {warning}" for warning in result.warnings)
        lines.append("")
    if honorable_mentions:
        lines.extend(
            [
                f"## {_HONORABLE_MENTION_HEADING}",
                "",
                f"The following {len(honorable_mentions)} gene(s) **did not survive genome-wide "
                "multiple-testing correction** (FDR-adjusted q-value at or above the "
                f"q={result.significance_level:g} significance threshold), but rank highest by "
                "raw Fisher p-value among the non-FDR-significant genes and may warrant targeted "
                "follow-up. This section is **not** a claim of statistical significance -- see "
                "the FDR-significant genes above for that.",
                "",
            ]
        )
        no_q_genes = [
            mention["gene_symbol"]
            for mention in honorable_mentions
            if mention["min_fdr_adjusted_q_value"] is None
        ]
        if no_q_genes:
            lines.extend(
                [
                    f"Note: {len(no_q_genes)} of the gene(s) above "
                    f"({', '.join(no_q_genes)}) never had a genome-wide FDR-adjusted "
                    "q-value computed at all (they contributed no p-value to the "
                    "Benjamini-Hochberg correction) -- for these genes only, the framing "
                    'above should be read as "FDR status unknown", not as a confirmed '
                    "q-value-at-or-above-threshold verdict; see that gene's own note "
                    "for precise wording.",
                    "",
                ]
            )
        lines.append("| " + " | ".join(_HONORABLE_MENTION_FIELDNAMES) + " |")
        lines.append("|" + "---|" * len(_HONORABLE_MENTION_FIELDNAMES))
        for mention in honorable_mentions:
            cells = " | ".join(
                _format_cell(mention[field]) for field in _HONORABLE_MENTION_FIELDNAMES
            )
            lines.append(f"| {cells} |")
        lines.append("")
    lines.extend(
        [
            "## Genome-wide summary plot",
            "",
            f"One point per scanned gene with an FDR-adjusted q-value: x-axis is genes "
            "ranked by composite evidence score (descending), y-axis is "
            "-log10(FDR-adjusted q-value), with a dashed line at the "
            f"q={result.significance_level:g} significance threshold.",
            "",
            f"![Genome-wide fusion-hotspot summary plot]({_MANHATTAN_SVG_FILENAME})",
            "",
        ]
    )
    lines.extend(["## Scanned genes (sorted by significance)", ""])
    lines.append("| " + " | ".join(_SUMMARY_FIELDNAMES) + " |")
    lines.append("|" + "---|" * len(_SUMMARY_FIELDNAMES))
    for row in rows:
        cells = " | ".join(_format_cell(row[field]) for field in _SUMMARY_FIELDNAMES)
        lines.append(f"| {cells} |")
    lines.append("")
    path.write_text("\n".join(lines))


def _write_summary_pdf(
    result: CohortScanResult, rows: list[dict], honorable_mentions: list[dict], path: Path
) -> None:
    header = [
        "Gene",
        "Source",
        "Status",
        "Patients",
        "In-frame%",
        "Domain-ret.%",
        "Best q",
        "Sig.",
        "Composite",
    ]
    table_rows = [header]
    for row in rows:
        table_rows.append(
            [
                row["gene_symbol"],
                row["config_source"],
                row["status"],
                _format_cell(row["distinct_patient_count"]),
                _format_cell(row["in_frame_percent"]),
                _format_cell(row["domain_retention_percent"]),
                _format_cell(row["min_fdr_adjusted_q_value"]),
                _format_cell(row["fdr_significant"]),
                _format_cell(row["top_composite_score"]),
            ]
        )

    extra_tables = []
    if honorable_mentions:
        mention_header = [
            "Rank",
            "Gene",
            "Fisher p",
            "Best q",
            "N events",
            "In-frame%",
            "Domain-ret.%",
        ]
        mention_rows = [mention_header] + [
            [
                _format_cell(mention["rank"]),
                mention["gene_symbol"],
                _format_cell(mention["fisher_p_value"]),
                _format_cell(mention["min_fdr_adjusted_q_value"]),
                _format_cell(mention["n_events_analyzed"]),
                _format_cell(mention["in_frame_percent"]),
                _format_cell(mention["domain_retention_percent"]),
            ]
            for mention in honorable_mentions
        ]
        extra_tables.append(
            {
                "heading": _HONORABLE_MENTION_HEADING,
                "note": (
                    "Did NOT survive genome-wide FDR correction -- ranked by raw Fisher "
                    "p-value among non-significant genes. Not a claim of significance."
                ),
                "rows": mention_rows,
            }
        )

    render_cohort_summary_pdf(
        path,
        title=f"Genome-wide fusion-hotspot cohort scan: {result.study_id}",
        subtitle=(
            f"{result.genes_after_gating} of {result.total_genes_before_gating} genes passed the "
            f">= {result.min_distinct_patients}-patient recurrence gate; "
            f"{len(result.significant_genes)} FDR-significant (q<0.05)."
        ),
        notes=[f"Generated {datetime.now(timezone.utc).isoformat()}"],
        rows=table_rows,
        extra_tables=extra_tables,
        figures_dir=path.parent,
    )


def _manuscript_gene_order(
    rows_by_gene: dict[str, dict],
    significant_genes: list[str],
    honorable_mentions: list[dict],
    highlighted_genes: "set[str] | frozenset[str]",
) -> list[str]:
    """Order the "Gene highlights" (and Appendix) section: FDR-significant
    genes first (most significant q-value first), then honorable-mention
    genes in their existing rank order, then any remaining hand-curated
    gene not already covered, alphabetically -- restricted throughout to
    ``highlighted_genes`` (the genes a full per-gene report was actually
    written for, see :func:`cfh.cohort.scan.genes_needing_full_report`), so
    every gene this manuscript highlights is guaranteed to have a real
    on-disk ``gene_reports/`` entry to link to.
    """

    def _q_sort_key(gene: str) -> float:
        q_value = rows_by_gene.get(gene, {}).get("min_fdr_adjusted_q_value")
        return q_value if q_value is not None else 1.0

    ordered: list[str] = []
    seen: set[str] = set()
    for gene in sorted((g for g in significant_genes if g in highlighted_genes), key=_q_sort_key):
        ordered.append(gene)
        seen.add(gene)
    for mention in honorable_mentions:
        gene = mention["gene_symbol"]
        if gene in highlighted_genes and gene not in seen:
            ordered.append(gene)
            seen.add(gene)
    for gene in sorted(highlighted_genes):
        if gene not in seen:
            ordered.append(gene)
            seen.add(gene)
    return ordered


def _gene_badges(
    gene: str,
    rows_by_gene: dict[str, dict],
    significant_genes: list[str],
    honorable_by_gene: dict[str, dict],
) -> list[str]:
    badges = []
    if gene in significant_genes:
        badges.append("FDR-significant")
    if gene in honorable_by_gene:
        badges.append("Honorable mention")
    if rows_by_gene.get(gene, {}).get("config_source") == "curated":
        badges.append("Curated gene config")
    return badges or ["Scanned"]


def _manuscript_key_figure(gene_report_paths: dict[str, Path] | None) -> Path | None:
    """Pick this gene's key figure for the manuscript highlight: the
    fusion-transcript schematic if one was generated for this gene, else
    the domain-retention lollipop/outlier diagram -- both already generated
    by :func:`cfh.real_benchmark.write_outputs` for this gene's own run, so
    nothing is regenerated here."""
    if not gene_report_paths:
        return None
    return gene_report_paths.get("fusion_schematic_svg") or gene_report_paths.get("domain_svg")


def _write_manuscript_markdown(
    payload: dict,
    gene_report_paths: dict[str, dict[str, Path]],
    path: Path,
) -> None:
    """Write the cross-gene manuscript-style synthesis report
    (``paper.md``): title, abstract, methods, results (Manhattan figure,
    FDR-significant/honorable-mention table, per-gene highlights with
    embedded figures), discussion caveats, and an appendix index into each
    highlighted gene's full existing per-gene report. Every sentence comes
    from :mod:`cfh.reporting.manuscript_text`; this function only handles
    Markdown assembly (headings, tables, image/link syntax) -- no numbers
    are computed here.
    """
    rows = payload.get("genes") or []
    honorable_mentions = payload.get("honorable_mentions") or []
    significant_genes = payload.get("significant_genes") or []
    rows_by_gene = {row["gene_symbol"]: row for row in rows}
    honorable_by_gene = {mention["gene_symbol"]: mention for mention in honorable_mentions}
    highlighted_genes = set(gene_report_paths)
    ordered_genes = _manuscript_gene_order(
        rows_by_gene, significant_genes, honorable_mentions, highlighted_genes
    )
    # The same denominator the Methods/Discussion sections report as having
    # entered the Benjamini-Hochberg correction (see
    # ``cfh.reporting.manuscript_text._scan_counts``'s ``n_with_q``) -- kept
    # in sync by deriving it here rather than re-stating a separate number.
    n_fdr_corrected_genes = sum(
        1 for row in rows if row.get("min_fdr_adjusted_q_value") is not None
    )

    lines = [
        f"# {render_manuscript_title(payload)}",
        "",
        "## Abstract",
        "",
        render_manuscript_abstract(payload),
        "",
        "## Methods",
        "",
        render_manuscript_methods(payload),
        "",
        "## Results",
        "",
        "### Genome-wide summary",
        "",
        render_manhattan_caption(payload),
        "",
        f"![Genome-wide fusion-hotspot summary plot]({_MANHATTAN_SVG_FILENAME})",
        "",
    ]

    table_fields = [
        "gene_symbol",
        "tier",
        "n_events_analyzed",
        "in_frame_percent",
        "domain_retention_percent",
        "fisher_p_value",
        "min_fdr_adjusted_q_value",
    ]
    significant_rows = sorted(
        (rows_by_gene[gene] for gene in significant_genes if gene in rows_by_gene),
        key=lambda row: (
            row.get("min_fdr_adjusted_q_value")
            if row.get("min_fdr_adjusted_q_value") is not None
            else 1.0
        ),
    )
    if significant_rows or honorable_mentions:
        lines.extend(["### FDR-significant and honorable-mention genes", ""])
        lines.append("| " + " | ".join(table_fields) + " |")
        lines.append("|" + "---|" * len(table_fields))
        for row in significant_rows:
            cells = {**row, "tier": "FDR-significant"}
            lines.append("| " + " | ".join(_format_cell(cells[f]) for f in table_fields) + " |")
        for mention in honorable_mentions:
            cells = {
                "gene_symbol": mention["gene_symbol"],
                "tier": "Honorable mention",
                "n_events_analyzed": mention["n_events_analyzed"],
                "in_frame_percent": mention["in_frame_percent"],
                "domain_retention_percent": mention["domain_retention_percent"],
                "fisher_p_value": mention["fisher_p_value"],
                "min_fdr_adjusted_q_value": mention["min_fdr_adjusted_q_value"],
            }
            lines.append("| " + " | ".join(_format_cell(cells[f]) for f in table_fields) + " |")
        lines.append("")

    lines.extend(["### Gene highlights", ""])
    for gene in ordered_genes:
        row = rows_by_gene.get(gene, {})
        badges = _gene_badges(gene, rows_by_gene, significant_genes, honorable_by_gene)
        lower = gene.lower()
        report_paths = gene_report_paths.get(gene) or {}
        report_md_rel = f"gene_reports/{lower}/report.md"
        note = honorable_by_gene.get(gene, {}).get("note")

        lines.append(f"#### {gene} ({', '.join(badges)})")
        lines.append("")
        lines.append(
            render_gene_highlight(
                row,
                n_fdr_corrected_genes=n_fdr_corrected_genes,
                honorable_mention_note=note,
            )
        )
        lines.append("")
        figure_path = _manuscript_key_figure(report_paths)
        if figure_path is not None:
            figure_rel = f"gene_reports/{lower}/visualizations/{figure_path.name}"
            lines.append(f"![{gene} key figure]({figure_rel})")
            lines.append("")
        lines.append(f"Full per-gene detail: [{report_md_rel}]({report_md_rel})")
        lines.append("")

    lines.extend(["## Discussion", ""])
    for bullet in render_discussion_bullets(payload):
        lines.append(f"- {bullet}")
    lines.append("")

    lines.extend(["## Appendix: per-gene report index", ""])
    appendix_fields = [
        "gene_symbol",
        "tier",
        "config_source",
        "n_events_analyzed",
        "min_fdr_adjusted_q_value",
        "report",
    ]
    lines.append("| " + " | ".join(appendix_fields) + " |")
    lines.append("|" + "---|" * len(appendix_fields))
    for gene in ordered_genes:
        row = rows_by_gene.get(gene, {})
        badges = _gene_badges(gene, rows_by_gene, significant_genes, honorable_by_gene)
        lower = gene.lower()
        report_md_rel = f"gene_reports/{lower}/report.md"
        cells = {
            "gene_symbol": gene,
            "tier": "; ".join(badges),
            "config_source": row.get("config_source"),
            "n_events_analyzed": row.get("n_events_analyzed"),
            "min_fdr_adjusted_q_value": row.get("min_fdr_adjusted_q_value"),
            "report": f"[{report_md_rel}]({report_md_rel})",
        }
        row_cells = [
            cells[field] if field == "report" else _format_cell(cells[field])
            for field in appendix_fields
        ]
        lines.append("| " + " | ".join(row_cells) + " |")
    lines.append("")

    path.write_text("\n".join(lines))


def _write_manuscript_pdf(
    payload: dict, gene_report_paths: dict[str, dict[str, Path]], path: Path
) -> None:
    """Write the cross-gene manuscript-style synthesis report
    (``paper.pdf``) via :func:`cfh.reporting.pdf.render_manuscript_pdf`,
    assembling the same figures/tables as :func:`_write_manuscript_markdown`
    from the same real data -- no numbers computed here."""
    rows = payload.get("genes") or []
    honorable_mentions = payload.get("honorable_mentions") or []
    significant_genes = payload.get("significant_genes") or []
    rows_by_gene = {row["gene_symbol"]: row for row in rows}
    honorable_by_gene = {mention["gene_symbol"]: mention for mention in honorable_mentions}
    highlighted_genes = set(gene_report_paths)
    ordered_genes = _manuscript_gene_order(
        rows_by_gene, significant_genes, honorable_mentions, highlighted_genes
    )
    n_fdr_corrected_genes = sum(
        1 for row in rows if row.get("min_fdr_adjusted_q_value") is not None
    )

    table_header = [
        "Gene",
        "Tier",
        "N events",
        "In-frame%",
        "Domain-ret.%",
        "Fisher p",
        "Best q",
    ]
    results_table_rows = [table_header]
    significant_rows = sorted(
        (rows_by_gene[gene] for gene in significant_genes if gene in rows_by_gene),
        key=lambda row: (
            row.get("min_fdr_adjusted_q_value")
            if row.get("min_fdr_adjusted_q_value") is not None
            else 1.0
        ),
    )
    for row in significant_rows:
        results_table_rows.append(
            [
                row["gene_symbol"],
                "FDR-significant",
                _format_cell(row["n_events_analyzed"]),
                _format_cell(row["in_frame_percent"]),
                _format_cell(row["domain_retention_percent"]),
                _format_cell(row["fisher_p_value"]),
                _format_cell(row["min_fdr_adjusted_q_value"]),
            ]
        )
    for mention in honorable_mentions:
        results_table_rows.append(
            [
                mention["gene_symbol"],
                "Honorable mention",
                _format_cell(mention["n_events_analyzed"]),
                _format_cell(mention["in_frame_percent"]),
                _format_cell(mention["domain_retention_percent"]),
                _format_cell(mention["fisher_p_value"]),
                _format_cell(mention["min_fdr_adjusted_q_value"]),
            ]
        )

    gene_highlights = []
    for gene in ordered_genes:
        row = rows_by_gene.get(gene, {})
        badges = _gene_badges(gene, rows_by_gene, significant_genes, honorable_by_gene)
        lower = gene.lower()
        report_paths = gene_report_paths.get(gene) or {}
        note = honorable_by_gene.get(gene, {}).get("note")
        figure_path = _manuscript_key_figure(report_paths)
        gene_highlights.append(
            {
                "heading": f"{gene} ({', '.join(badges)})",
                "paragraph": render_gene_highlight(
                    row,
                    n_fdr_corrected_genes=n_fdr_corrected_genes,
                    honorable_mention_note=note,
                ),
                "figure_path": figure_path,
                "figure_caption": (
                    f"Reused from the {gene} individual gene report ({figure_path.name})."
                    if figure_path is not None
                    else None
                ),
                "report_note": f"Full per-gene detail: gene_reports/{lower}/report.md",
            }
        )

    appendix_header = ["Gene", "Tier", "Config source", "N events", "Best q", "Full report"]
    appendix_rows = [appendix_header]
    for gene in ordered_genes:
        row = rows_by_gene.get(gene, {})
        badges = _gene_badges(gene, rows_by_gene, significant_genes, honorable_by_gene)
        lower = gene.lower()
        appendix_rows.append(
            [
                gene,
                "; ".join(badges),
                _format_cell(row.get("config_source")),
                _format_cell(row.get("n_events_analyzed")),
                _format_cell(row.get("min_fdr_adjusted_q_value")),
                f"gene_reports/{lower}/report.md",
            ]
        )

    render_manuscript_pdf(
        path,
        title=render_manuscript_title(payload),
        abstract=render_manuscript_abstract(payload),
        methods=render_manuscript_methods(payload),
        manhattan_svg_path=path.parent / _MANHATTAN_SVG_FILENAME,
        manhattan_caption=render_manhattan_caption(payload),
        results_table_rows=results_table_rows if len(results_table_rows) > 1 else [],
        gene_highlights=gene_highlights,
        discussion_bullets=render_discussion_bullets(payload),
        appendix_rows=appendix_rows if len(appendix_rows) > 1 else [],
    )


def write_cohort_scan_outputs(
    result: CohortScanResult,
    output_dir: str | Path,
    *,
    run_id: str | None = None,
    pdf: bool = True,
    honorable_mention_count: int = DEFAULT_HONORABLE_MENTION_COUNT,
) -> dict[str, Path]:
    """Write the consolidated summary and full per-gene reports for a
    completed cohort scan to ``<output_dir>/<run_id>/cohort_scan/``."""
    generated_at = datetime.now(timezone.utc).isoformat()
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    run_id = run_id or f"cohort-scan_{result.study_id}_{timestamp}"
    destination = Path(output_dir) / run_id / "cohort_scan"
    destination.mkdir(parents=True, exist_ok=True)

    rows = build_summary_rows(result)
    honorable_mentions = build_honorable_mentions(rows, limit=honorable_mention_count)
    # Derived from this run's own results -- never a live import of
    # ``cfh.algorithms.registry`` -- so it stays correct for a programmatic
    # run made with a restricted ``algorithm_names`` list, and so
    # re-rendering this run's summary.json later, after the registry has
    # changed, still describes what THIS run actually did.
    algorithms_run = sorted(
        {
            algorithm_result.Algorithm
            for outcome in result.gene_outcomes
            if outcome.run is not None
            for algorithm_result in outcome.run.results
        }
    )

    curated_genes = {
        outcome.gene_symbol
        for outcome in result.gene_outcomes
        if outcome.config_source == "curated"
    }
    honorable_mention_genes = {mention["gene_symbol"] for mention in honorable_mentions}
    priority_genes = curated_genes | honorable_mention_genes

    manhattan_svg_path = destination / _MANHATTAN_SVG_FILENAME
    manhattan_svg_path.write_text(
        render_manhattan_svg(
            rows, significance_level=result.significance_level, priority_genes=priority_genes
        )
        + "\n"
    )

    tsv_path = destination / "summary.tsv"
    with tsv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=_SUMMARY_FIELDNAMES,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)

    json_path = destination / "summary.json"
    payload = _json_safe(
        {
            "study_id": result.study_id,
            "min_distinct_patients": result.min_distinct_patients,
            "total_genes_before_gating": result.total_genes_before_gating,
            "genes_after_gating": result.genes_after_gating,
            "curated_gene_count": result.curated_gene_count,
            "auto_config_gene_count": result.auto_config_gene_count,
            "unresolved_gene_count": result.unresolved_gene_count,
            "significant_genes": result.significant_genes,
            "significance_level": result.significance_level,
            "generated_at": generated_at,
            "algorithms_run": algorithms_run,
            "warnings": result.warnings,
            "genes": rows,
            "honorable_mentions": honorable_mentions,
        }
    )
    json_path.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")

    markdown_path = destination / "summary.md"
    _write_summary_markdown(result, rows, honorable_mentions, markdown_path)

    paths = {
        "run_directory": destination,
        "summary_tsv": tsv_path,
        "summary_json": json_path,
        "summary_markdown": markdown_path,
        "manhattan_svg": manhattan_svg_path,
    }

    if pdf:
        pdf_path = destination / "summary.pdf"
        _write_summary_pdf(result, rows, honorable_mentions, pdf_path)
        paths["summary_pdf"] = pdf_path

    gene_reports_dir = destination / "gene_reports"
    outcomes_by_gene = {outcome.gene_symbol: outcome for outcome in result.gene_outcomes}
    full_report_genes = genes_needing_full_report(result, honorable_mention_genes)
    gene_report_paths: dict[str, dict[str, Path]] = {}
    for gene_symbol in full_report_genes:
        outcome = outcomes_by_gene.get(gene_symbol)
        if outcome is None or outcome.run is None:
            continue
        gene_report_paths[gene_symbol] = write_outputs(
            outcome.run,
            gene_reports_dir,
            run_id=gene_symbol.lower(),
            pdf=pdf,
            cli_args=["cohort-scan", result.study_id],
        )
    paths["gene_reports"] = gene_report_paths

    manuscript_payload = {
        "study_id": result.study_id,
        "min_distinct_patients": result.min_distinct_patients,
        "total_genes_before_gating": result.total_genes_before_gating,
        "genes_after_gating": result.genes_after_gating,
        "curated_gene_count": result.curated_gene_count,
        "auto_config_gene_count": result.auto_config_gene_count,
        "unresolved_gene_count": result.unresolved_gene_count,
        "significant_genes": result.significant_genes,
        "significance_level": result.significance_level,
        "generated_at": generated_at,
        "algorithms_run": algorithms_run,
        "genes": rows,
        "honorable_mentions": honorable_mentions,
    }
    manuscript_markdown_path = destination / _MANUSCRIPT_MARKDOWN_FILENAME
    _write_manuscript_markdown(manuscript_payload, gene_report_paths, manuscript_markdown_path)
    paths["manuscript_markdown"] = manuscript_markdown_path

    if pdf:
        manuscript_pdf_path = destination / _MANUSCRIPT_PDF_FILENAME
        _write_manuscript_pdf(manuscript_payload, gene_report_paths, manuscript_pdf_path)
        paths["manuscript_pdf"] = manuscript_pdf_path

    return paths
