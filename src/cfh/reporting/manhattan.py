"""Static, gene-agnostic genome-wide Manhattan/volcano summary plot for a
``cfh cohort-scan`` run, drawn as inline raw SVG in the same style already
used for the per-gene lollipop/domain-track figures in
:mod:`cfh.real_benchmark` (``<rect>``/``<circle>``/``<line>``/``<text>``
only -- no new plotting-library dependency).

Input is exactly the list of row dicts :func:`cfh.cohort.outputs.build_summary_rows`
already produces (or the equivalent ``genes`` list read back from a written
``summary.json``) -- no gene names are hardcoded here, and no statistics are
computed here; every value plotted was already computed by the cohort-scan
pipeline and the WP13 cross-gene Benjamini-Hochberg correction
(:mod:`cfh.stats.multiple_testing`).

Axes:

* y: ``-log10`` of each gene's ``min_fdr_adjusted_q_value`` (its smallest
  BH-adjusted q-value across whichever algorithms/tests it produced), with a
  horizontal dashed line at ``-log10(significance_level)`` (q=0.05 by
  default). Genes with no q-value at all -- i.e. they never contributed a
  usable p-value to the cross-gene correction family -- have no y-value and
  are not plotted; every other scanned gene appears as exactly one point.
* x: genes ranked left-to-right by effect size, strongest evidence first.
  The effect-size proxy used is each gene's ``top_composite_score`` (its
  best partner-level score from the ``composite_score`` algorithm) --
  already computed per gene by ``build_summary_rows``, already available
  for exactly the same genes that have a q-value (a raw log-odds-ratio
  would work equally well as the ranking axis; composite score was chosen
  because it combines evidence across algorithms rather than depending on a
  single test, and needed no new statistic to be computed for this plot).

Points at/above the FDR line (q < ``significance_level``) are drawn as
filled red circles; points below it as smaller blue circles -- color and
marker size both encode the same significance call already present in each
row's ``fdr_significant`` flag. Each plotted circle carries a
``data-gene``/``data-significant`` attribute so the placement of a specific
gene can be checked programmatically from the rendered SVG, and the
threshold line carries ``id="fdr-threshold-line"`` for the same reason.

Gene-symbol labels are drawn for every FDR-significant point, plus (up to
``max_labels`` genes total) the remaining most-significant points --
avoiding label clutter on a genome-wide scan of hundreds of genes. A
genome-wide scan analyzes fusion PARTNER genes as their own independent
scanned genes (e.g. ``EML4``, ``KIAA1549``, ``TACC3`` each get their own row,
alongside the primary driver genes ``ALK``, ``BRAF``, ``NTRK1``), and a
partner gene's row can happen to reach a smaller raw q-value than a
biologically more interesting primary driver simply because they share the
same underlying breakpoint events. Left unweighted, filling the label budget
by raw q-value alone lets those partner-gene rows crowd out hand-curated
driver genes and top near-miss genes from the limited label slots. Callers
should therefore pass ``priority_genes`` (hand-curated genes plus the
"honorable mentions" highly ranked non-FDR-significant tier -- see
:func:`cfh.cohort.outputs.build_honorable_mentions`) so those genes' labels
are filled first, ahead of any other non-significant gene, regardless of
where they land in the raw q-value ranking. Each point's underlying
``top_composite_partner_gene`` (if any) is still preserved -- as an SVG
``<title>`` tooltip on its circle -- even when that gene's own label loses
its slot.
"""

from __future__ import annotations

import math

from cfh.reporting.svg_utils import escape_xml_text

_WIDTH = 960
_HEIGHT = 520
_MARGIN_LEFT = 70
_MARGIN_RIGHT = 30
_MARGIN_TOP = 50
_MARGIN_BOTTOM = 70

_SIGNIFICANT_COLOR = "#d62728"
_NOT_SIGNIFICANT_COLOR = "#2878b5"
_LABEL_Y_OFFSETS = (-10, -22, -34)

DEFAULT_MAX_LABELS = 15


def _plottable_points(rows: list[dict]) -> list[dict]:
    """Genes with a real FDR-adjusted q-value, ranked by descending
    composite evidence score (the plot's x axis)."""
    plottable = [
        row
        for row in rows
        if row.get("min_fdr_adjusted_q_value") is not None
        and row.get("top_composite_score") is not None
    ]
    plottable.sort(key=lambda row: row["top_composite_score"], reverse=True)
    return plottable


def _label_indices(
    points: list[dict], max_labels: int, priority_gene_symbols: frozenset[str] = frozenset()
) -> set[int]:
    """Every FDR-significant point is always labeled. Remaining label slots
    (``max_labels`` minus the significant count) go first to non-significant
    points whose gene is in ``priority_gene_symbols`` (ranked by q-value
    among themselves), then to the next most-significant points overall --
    so a hand-curated or near-miss "honorable mention" gene never loses its
    label slot to an incidentally-significant fusion-partner gene."""
    significant_idx = [index for index, point in enumerate(points) if point["significant"]]
    remaining = [index for index in range(len(points)) if not points[index]["significant"]]

    def by_neg_log_q(index: int) -> float:
        return points[index]["neg_log_q"]

    priority_idx = sorted(
        (index for index in remaining if points[index]["gene_symbol"] in priority_gene_symbols),
        key=by_neg_log_q,
        reverse=True,
    )
    other_idx = sorted(
        (index for index in remaining if points[index]["gene_symbol"] not in priority_gene_symbols),
        key=by_neg_log_q,
        reverse=True,
    )
    slots = max(max_labels - len(significant_idx), 0)
    ordered_candidates = priority_idx + other_idx
    return set(significant_idx) | set(ordered_candidates[:slots])


def render_manhattan_svg(
    rows: list[dict],
    *,
    significance_level: float = 0.05,
    max_labels: int = DEFAULT_MAX_LABELS,
    priority_genes: "frozenset[str] | set[str]" = frozenset(),
) -> str:
    """Render the genome-wide Manhattan/volcano-style summary SVG for one
    cohort scan's already-built summary rows. Returns the raw SVG markup
    (no trailing newline)."""
    if not 0.0 < significance_level < 1.0:
        raise ValueError(f"significance_level must be in (0, 1); got {significance_level!r}")

    plot_left, plot_right = _MARGIN_LEFT, _WIDTH - _MARGIN_RIGHT
    plot_top, plot_bottom = _MARGIN_TOP, _HEIGHT - _MARGIN_BOTTOM
    plot_width = plot_right - plot_left
    plot_height = plot_bottom - plot_top

    ranked = _plottable_points(rows)
    points = [
        {
            "gene_symbol": row["gene_symbol"],
            "neg_log_q": -math.log10(max(row["min_fdr_adjusted_q_value"], 1e-300)),
            "significant": bool(row.get("fdr_significant")),
            "partner_gene": row.get("top_composite_partner_gene"),
        }
        for row in ranked
    ]

    threshold_neg_log_q = -math.log10(significance_level)
    y_max = max([threshold_neg_log_q, *(point["neg_log_q"] for point in points)]) * 1.1
    y_max = max(y_max, 0.1)

    def x_for_rank(index: int) -> float:
        if len(points) <= 1:
            return plot_left + plot_width / 2
        return plot_left + plot_width * index / (len(points) - 1)

    def y_for_value(value: float) -> float:
        return plot_bottom - min(value, y_max) / y_max * plot_height

    elements = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{_WIDTH}" height="{_HEIGHT}" '
        f'viewBox="0 0 {_WIDTH} {_HEIGHT}">',
        f'<rect width="{_WIDTH}" height="{_HEIGHT}" fill="white"/>',
        f'<text x="{plot_left}" y="20" font-family="sans-serif" font-size="16">'
        "Genome-wide fusion-hotspot summary: FDR significance vs. composite evidence rank"
        "</text>",
        f'<line x1="{plot_left}" y1="{plot_top}" x2="{plot_left}" y2="{plot_bottom}" '
        'stroke="#444" stroke-width="1.5"/>',
        f'<line x1="{plot_left}" y1="{plot_bottom}" x2="{plot_right}" y2="{plot_bottom}" '
        'stroke="#444" stroke-width="1.5"/>',
        f'<text x="{(plot_left + plot_right) / 2:.1f}" y="{_HEIGHT - 18}" '
        'font-family="sans-serif" font-size="12" text-anchor="middle">'
        "Genes ranked by composite evidence score (highest first)</text>",
        f'<text x="18" y="{(plot_top + plot_bottom) / 2:.1f}" font-family="sans-serif" '
        'font-size="12" text-anchor="middle" '
        f'transform="rotate(-90 18 {(plot_top + plot_bottom) / 2:.1f})">'
        "-log10(FDR-adjusted q-value)</text>",
    ]

    if not points:
        elements.append(
            f'<text x="{(plot_left + plot_right) / 2:.1f}" y="{(plot_top + plot_bottom) / 2:.1f}" '
            'font-family="sans-serif" font-size="13" text-anchor="middle" fill="#888">'
            "No scanned gene produced an FDR-adjusted q-value in this run.</text>"
        )

    threshold_y = y_for_value(threshold_neg_log_q)
    elements.append(
        f'<line id="fdr-threshold-line" x1="{plot_left}" y1="{threshold_y:.2f}" '
        f'x2="{plot_right}" y2="{threshold_y:.2f}" stroke="#888" stroke-width="1.2" '
        'stroke-dasharray="6,4"/>'
    )
    elements.append(
        f'<text x="{plot_right}" y="{threshold_y - 5:.1f}" font-family="sans-serif" '
        f'font-size="11" text-anchor="end" fill="#666">q = {significance_level:g} threshold</text>'
    )

    for index, point in enumerate(points):
        x = x_for_rank(index)
        y = y_for_value(point["neg_log_q"])
        if point["significant"]:
            color, radius = _SIGNIFICANT_COLOR, 4.0
        else:
            color, radius = _NOT_SIGNIFICANT_COLOR, 2.6
        circle_attrs = (
            f'cx="{x:.2f}" cy="{y:.2f}" r="{radius}" fill="{color}" '
            f'data-gene="{point["gene_symbol"]}" '
            f'data-significant="{"true" if point["significant"] else "false"}"'
        )
        if point["partner_gene"]:
            title = escape_xml_text(
                f"{point['gene_symbol']}: top composite-evidence partner gene "
                f"{point['partner_gene']}"
            )
            elements.append(f"<circle {circle_attrs}><title>{title}</title></circle>")
        else:
            elements.append(f"<circle {circle_attrs}/>")

    priority_gene_symbols = frozenset(priority_genes)
    labeled_indices = sorted(
        _label_indices(points, max_labels, priority_gene_symbols), key=x_for_rank
    )
    for order, index in enumerate(labeled_indices):
        point = points[index]
        x = x_for_rank(index)
        y = y_for_value(point["neg_log_q"])
        dy = _LABEL_Y_OFFSETS[order % len(_LABEL_Y_OFFSETS)]
        elements.append(
            f'<text x="{x:.2f}" y="{y + dy:.2f}" font-family="sans-serif" font-size="9" '
            f'text-anchor="middle">{point["gene_symbol"]}</text>'
        )

    elements.extend(
        [
            f'<circle cx="{plot_left + 10}" cy="38" r="4" fill="{_SIGNIFICANT_COLOR}"/>',
            f'<text x="{plot_left + 20}" y="42" font-family="sans-serif" '
            'font-size="11">FDR-significant (q &lt; significance threshold)</text>',
            f'<circle cx="{plot_left + 280}" cy="38" r="3" fill="{_NOT_SIGNIFICANT_COLOR}"/>',
            f'<text x="{plot_left + 290}" y="42" font-family="sans-serif" '
            'font-size="11">Not significant</text>',
        ]
    )

    elements.append("</svg>")
    return "\n".join(elements)
