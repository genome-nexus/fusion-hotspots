"""Static SVG fusion-transcript schematics, in the style of Figure 4 panels
B and C of Zehir et al. (PMC5461196): one horizontal row per (recurrent
group of) fusion event, all rows sharing one amino-acid-position x-axis for
the target gene's full protein.

Pure rendering only, same discipline as :mod:`cfh.reporting.pdf`/
:mod:`cfh.reporting.text`: every function here is a function of an
already-computed run payload (the same dict shape written to and read back
from ``results.json``) plus its ``gene_track``/``intragenic_deletions``
fields -- nothing here makes a network call or computes a new statistic.
Gene-agnostic by construction: nothing below references a gene symbol
literally; every gene-specific fact (protein length, domain boundaries,
exon boundaries, partner names, breakpoints) comes from the payload.

Colors for domain-retention status (``RETAINED_COLOR``/``TRUNCATED_COLOR``)
and the breakpoint marker (``BREAKPOINT_COLOR``) are imported from
:mod:`cfh.reporting.palette`, the single source of truth also used by the
existing per-event domain-retention lollipop track
(``cfh.real_benchmark._domain_track_svg``) -- this module does not define
its own copies of those hex values, so the two renderers cannot drift
apart on what a color means. Partner-gene colors and the neutral domain
backbone are new (the lollipop track has neither per-partner nor
backbone-fill colors to reuse), computed deterministically so a given
partner always renders the same color within and across runs.
"""

from __future__ import annotations

from cfh.reporting.domain_names import domain_label_for_accession
from cfh.reporting.palette import (
    BREAKPOINT_COLOR,
    RETAINED_COLOR,
    TRUNCATED_COLOR,
    deterministic_color,
)

__all__ = [
    "RETAINED_COLOR",
    "TRUNCATED_COLOR",
    "BREAKPOINT_COLOR",
    "BACKBONE_COLOR",
    "AXIS_COLOR",
    "CONNECTOR_COLOR",
    "partner_color",
    "exon_boundary_ticks_svg",
    "render_position_axis_svg",
    "render_fusion_schematic_svg",
    "render_intragenic_deletion_schematic_svg",
]

# New to this module (the lollipop track has no backbone/partner fill).
BACKBONE_COLOR = "#e2e2e2"
AXIS_COLOR = "#444444"
CONNECTOR_COLOR = "#999999"

_MAX_ROWS_DEFAULT = 28
"""Matches the source paper's own Figure 4B row count (~28 rows) rather
than attempting to draw one row per raw event on a cohort with hundreds."""

_ROW_HEIGHT = 22
_ROW_GAP = 4
_AXIS_LEFT = 60
_AXIS_WIDTH = 560
_LABEL_LEFT = _AXIS_LEFT + _AXIS_WIDTH + 40
_SVG_WIDTH = 1000
_EXON_TICK_LABEL_HEIGHT = 20
"""Extra vertical room reserved below the shared position axis for its
rotated exon-number labels (see :func:`exon_boundary_ticks_svg`)."""
_TOP_MARGIN = 56
_LEGEND_HEIGHT = 30
_BOTTOM_MARGIN = 40


def partner_color(partner_gene: str) -> str:
    """Deterministic, arbitrary-but-stable color for a partner gene name.

    Same partner always gets the same color within a run and across runs
    (a pure function of the name), so a reader can visually track one
    partner across rows. Thin wrapper around the shared
    :func:`cfh.reporting.palette.deterministic_color` hash so this
    module's own hue-from-name formula can't silently diverge from
    :mod:`cfh.real_benchmark`'s per-domain highlight coloring.
    """
    return deterministic_color(partner_gene)


def _clip(value: float, low: float, high: float) -> float:
    return max(low, min(high, value))


def _domain_color_segments(
    domains: list[dict], interval_start: float, interval_end: float
) -> list[tuple[float, float, str, str]]:
    """Sub-intervals of ``domains`` overlapping ``[interval_start, interval_end]``.

    A domain wholly inside the interval is drawn ``RETAINED_COLOR`` in
    full; a domain straddling one edge of the interval is drawn
    ``TRUNCATED_COLOR`` for only its overlapping (retained) portion. A
    domain entirely outside the interval contributes nothing -- it belongs
    to whatever is on the other side of the breakpoint/deletion, which is
    not this gene's own sequence in the resulting fusion/deletion protein.
    """
    segments: list[tuple[float, float, str, str]] = []
    for domain in domains:
        start, end = domain.get("start_aa"), domain.get("end_aa")
        if start is None or end is None:
            continue
        overlap_start = max(start, interval_start)
        overlap_end = min(end, interval_end)
        if overlap_start > overlap_end:
            continue
        fully_contained = start >= interval_start and end <= interval_end
        color = RETAINED_COLOR if fully_contained else TRUNCATED_COLOR
        segments.append((overlap_start, overlap_end, color, domain.get("name") or "domain"))
    return segments


def _exon_ticks(
    exon_boundaries: list[dict], interval_start: float, interval_end: float
) -> list[float]:
    """Exon-start positions from ``exon_boundaries`` that fall strictly
    inside ``[interval_start, interval_end]`` (the interval's own edges
    already get a breakpoint/axis line, so they don't need a duplicate
    tick)."""
    ticks = []
    for boundary in exon_boundaries:
        start = boundary.get("start_aa")
        if start is None:
            continue
        if interval_start < start < interval_end:
            ticks.append(start)
    return sorted(set(ticks))


def exon_boundary_ticks_svg(
    exon_boundaries: list[dict],
    *,
    axis_left: float,
    scale: float,
    y: float,
    max_position: float | None = None,
    tick_length: float = 4.0,
) -> list[str]:
    """One small tick + rotated ``E<rank>`` label per exon boundary in
    ``exon_boundaries`` (gene_track's ``exon_boundaries_aa`` field, which
    carries each exon's own ``exon_rank`` -- the real ordinal exon number
    from the canonical transcript, not something derived/guessed here),
    against a position axis already scaled the same way as the rest of a
    diagram (``x = axis_left + position * scale``).

    Drawn once along a diagram's shared axis rather than once per row:
    exon boundaries are a property of the target gene's canonical
    transcript, not of any individual fusion/breakpoint row, so repeating
    the same label on every row (up to 28 in the fusion schematic) would
    only add clutter, not information. Exported so
    :func:`cfh.real_benchmark._domain_track_svg`'s own x-axis draws exon
    boundaries with this exact convention instead of maintaining a second,
    potentially-drifting copy of it.
    """
    elements: list[str] = []
    seen_positions: set[float] = set()
    for boundary in exon_boundaries:
        start = boundary.get("start_aa")
        rank = boundary.get("exon_rank")
        if start is None or rank is None or start in seen_positions:
            continue
        if max_position is not None and not (0 <= start <= max_position):
            continue
        seen_positions.add(start)
        x = axis_left + start * scale
        elements.append(
            f'<line x1="{x:.1f}" y1="{y:.1f}" x2="{x:.1f}" y2="{y + tick_length:.1f}" '
            'stroke="#777777" stroke-width="0.7"/>'
        )
        label_y = y + tick_length + 1
        elements.append(
            f'<text x="{x:.1f}" y="{label_y:.1f}" font-family="sans-serif" font-size="5.5" '
            f'fill="#666666" transform="rotate(55 {x:.1f} {label_y:.1f})">E{rank}</text>'
        )
    return elements


def _fit_domain_label(name: str, available_width: float, *, font_size: float = 7.0) -> str | None:
    """Best-effort short label for a domain-colored segment: the full
    domain name when it fits the segment's rendered pixel width, a
    truncated ``…``-suffixed prefix when it doesn't, or ``None`` (skip
    entirely) when the segment is too narrow for even a short
    abbreviation -- stacking unreadable overlapping glyphs would be worse
    than no label at all."""
    avg_char_width = font_size * 0.58
    max_chars = int(available_width // avg_char_width)
    if max_chars < 3:
        return None
    if len(name) <= max_chars:
        return name
    return name[: max_chars - 1].rstrip() + "…"


def _status_word(status: str | None) -> str:
    return {
        "retained": "retained",
        "disrupted": "truncated",
        "lost": "lost",
    }.get(status or "", "unknown")


def _fusion_groups(payload: dict) -> list[dict]:
    """Group per-event rows into one schematic row per (partner, breakpoint,
    role): identical repeats (same partner gene at the exact same
    protein-position breakpoint) collapse into one row with a count, same
    as the paper's own ``PARTNER (xN)`` convention. Distinct breakpoints
    for the same partner stay separate rows -- they are visually different
    fusions, not repeats of one."""
    groups: dict[tuple[str, int, str], dict] = {}
    for row in payload.get("events") or []:
        partner = row.get("partner_gene")
        breakpoint_aa = row.get("breakpoint_protein_position")
        role = row.get("target_role")
        if not partner or breakpoint_aa is None or role not in {"five_prime", "three_prime"}:
            continue
        key = (partner, int(breakpoint_aa), role)
        group = groups.setdefault(
            key,
            {
                "partner_gene": partner,
                "breakpoint_aa": int(breakpoint_aa),
                "role": role,
                "count": 0,
                "sample_ids": [],
                "domain_status": row.get("domain_status"),
                "tumor_types": [],
            },
        )
        group["count"] += 1
        if row.get("sample_id"):
            group["sample_ids"].append(row["sample_id"])
        tumor_type = row.get("tumor_type") or row.get("Tumor_type")
        if tumor_type:
            group["tumor_types"].append(tumor_type)
    return list(groups.values())


def _tumor_type_clause(tumor_types: list[str]) -> str | None:
    if not tumor_types:
        return None
    counts: dict[str, int] = {}
    for tumor_type in tumor_types:
        counts[tumor_type] = counts.get(tumor_type, 0) + 1
    distinct = sorted(counts, key=lambda t: (-counts[t], t))
    if len(distinct) > 3:
        return "multiple tumor types"
    return ", ".join(distinct)


def _row_label(
    group: dict,
    domain_label: str | None,
) -> str:
    partner = group["partner_gene"]
    count = group["count"]
    name = partner if count == 1 else f"{partner} (x{count})"
    status = _status_word(group.get("domain_status"))
    parts = [name]
    if domain_label:
        parts.append(f"{domain_label} {status}")
    tumor_clause = _tumor_type_clause(group.get("tumor_types") or [])
    if tumor_clause:
        parts.append(tumor_clause)
    return " – ".join(parts)


def _axis_scale(protein_length: int) -> float:
    return _AXIS_WIDTH / max(protein_length, 1)


def _svg_open(width: int, height: int) -> str:
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">\n'
        f'<rect width="{width}" height="{height}" fill="white"/>'
    )


def render_position_axis_svg(
    y: float,
    protein_length: int,
    scale: float,
    *,
    axis_left: float = _AXIS_LEFT,
    exon_boundaries: list[dict] | None = None,
) -> list[str]:
    """The shared amino-acid-position axis drawn once at the bottom of a
    schematic: the 1/protein-length endpoint labels (existing behavior),
    plus -- when ``exon_boundaries`` (gene_track's ``exon_boundaries_aa``)
    is given -- one tick and ``E<rank>`` label per exon boundary, via
    :func:`exon_boundary_ticks_svg`. Exported so
    :func:`cfh.real_benchmark._domain_track_svg` can draw its own x-axis
    with the exact same exon-labeling convention.
    """
    x0 = axis_left
    x1 = axis_left + protein_length * scale
    elements = [
        f'<line x1="{x0:.1f}" y1="{y:.1f}" x2="{x1:.1f}" y2="{y:.1f}" '
        f'stroke="{AXIS_COLOR}" stroke-width="1"/>',
        f'<text x="{x0:.1f}" y="{y + 14:.1f}" font-family="sans-serif" font-size="9">1</text>',
        f'<text x="{x1:.1f}" y="{y + 14:.1f}" font-family="sans-serif" font-size="9" '
        f'text-anchor="end">{protein_length} aa</text>',
    ]
    elements.extend(
        exon_boundary_ticks_svg(
            exon_boundaries or [],
            axis_left=axis_left,
            scale=scale,
            y=y + 18,
            max_position=protein_length,
        )
    )
    return elements


def _legend_svg(y: float) -> list[str]:
    entries = [
        (RETAINED_COLOR, "domain retained"),
        (TRUNCATED_COLOR, "domain truncated"),
        (BACKBONE_COLOR, "target backbone (no domain)"),
        (BREAKPOINT_COLOR, "breakpoint"),
    ]
    elements = []
    x = _AXIS_LEFT
    for color, label in entries:
        elements.append(
            f'<rect x="{x:.1f}" y="{y - 9:.1f}" width="12" height="10" fill="{color}"/>'
        )
        elements.append(
            f'<text x="{x + 16:.1f}" y="{y:.1f}" font-family="sans-serif" font-size="9">'
            f"{label}</text>"
        )
        x += 18 + 9 * len(label) + 14
    return elements


def render_fusion_schematic_svg(payload: dict, *, max_rows: int = _MAX_ROWS_DEFAULT) -> str | None:
    """Render the Figure-4B-style fusion schematic for one run's payload.

    Returns ``None`` (never raises, never fabricates data) when the run has
    no ``gene_track`` (no protein length to build a shared axis from) or no
    mappable fusion rows to draw.
    """
    gene_track = payload.get("gene_track")
    if not gene_track or not gene_track.get("protein_length"):
        return None
    protein_length = gene_track["protein_length"]
    domains = gene_track.get("domains") or []
    exon_boundaries = gene_track.get("exon_boundaries_aa") or []

    groups = _fusion_groups(payload)
    if not groups:
        return None
    groups.sort(key=lambda g: (-g["count"], g["breakpoint_aa"], g["partner_gene"]))
    total_groups = len(groups)
    shown = groups[:max_rows]
    truncated_count = total_groups - len(shown)

    summary = payload.get("summary") or {}
    domain_label = domain_label_for_accession(domains, summary.get("domain_accession"))

    scale = _axis_scale(protein_length)
    n_rows = len(shown)
    plot_height = n_rows * (_ROW_HEIGHT + _ROW_GAP)
    truncation_note_height = 16 if truncated_count else 0
    height = (
        _TOP_MARGIN
        + plot_height
        + _EXON_TICK_LABEL_HEIGHT
        + _LEGEND_HEIGHT
        + truncation_note_height
        + _BOTTOM_MARGIN
    )

    gene_symbol = payload.get("gene_symbol") or "target gene"
    elements = [_svg_open(_SVG_WIDTH, height)]
    elements.append(
        f'<text x="{_AXIS_LEFT}" y="24" font-family="sans-serif" font-size="15">'
        f"{gene_symbol} fusion-transcript schematic</text>"
    )
    elements.append(
        f'<text x="{_AXIS_LEFT}" y="40" font-family="sans-serif" font-size="10" fill="#555">'
        f"partner block → breakpoint → retained {gene_symbol} portion "
        "(domain-colored, exon ticks)</text>"
    )

    y = _TOP_MARGIN
    for group in shown:
        row_top = y
        row_mid = y + _ROW_HEIGHT / 2
        breakpoint_aa = _clip(group["breakpoint_aa"], 0, protein_length)
        if group["role"] == "three_prime":
            partner_span = (0, breakpoint_aa)
            target_span = (breakpoint_aa, protein_length)
        else:
            target_span = (0, breakpoint_aa)
            partner_span = (breakpoint_aa, protein_length)

        color = partner_color(group["partner_gene"])
        p_x0 = _AXIS_LEFT + partner_span[0] * scale
        p_x1 = _AXIS_LEFT + partner_span[1] * scale
        elements.append(
            f'<rect x="{p_x0:.1f}" y="{row_top:.1f}" width="{max(0.5, p_x1 - p_x0):.1f}" '
            f'height="{_ROW_HEIGHT}" fill="{color}"/>'
        )

        t_x0 = _AXIS_LEFT + target_span[0] * scale
        t_x1 = _AXIS_LEFT + target_span[1] * scale
        elements.append(
            f'<rect x="{t_x0:.1f}" y="{row_top:.1f}" width="{max(0.5, t_x1 - t_x0):.1f}" '
            f'height="{_ROW_HEIGHT}" fill="{BACKBONE_COLOR}"/>'
        )
        for seg_start, seg_end, seg_color, seg_name in _domain_color_segments(
            domains, target_span[0], target_span[1]
        ):
            sx0 = _AXIS_LEFT + seg_start * scale
            sx1 = _AXIS_LEFT + seg_end * scale
            seg_width = max(0.5, sx1 - sx0)
            elements.append(
                f'<rect x="{sx0:.1f}" y="{row_top:.1f}" width="{seg_width:.1f}" '
                f'height="{_ROW_HEIGHT}" fill="{seg_color}"/>'
            )
            seg_label = _fit_domain_label(seg_name, seg_width)
            if seg_label:
                elements.append(
                    f'<text x="{(sx0 + sx1) / 2:.1f}" y="{row_mid + 3:.1f}" '
                    'font-family="sans-serif" font-size="7" text-anchor="middle" '
                    'fill="#111111" stroke="white" stroke-width="2" paint-order="stroke">'
                    f"{seg_label}</text>"
                )
        for tick_aa in _exon_ticks(exon_boundaries, target_span[0], target_span[1]):
            tx = _AXIS_LEFT + tick_aa * scale
            tick_bottom = row_top + _ROW_HEIGHT
            elements.append(
                f'<line x1="{tx:.1f}" y1="{row_top:.1f}" x2="{tx:.1f}" y2="{tick_bottom:.1f}" '
                'stroke="#555555" stroke-width="0.6" opacity="0.6"/>'
            )

        bx = _AXIS_LEFT + breakpoint_aa * scale
        breakpoint_top, breakpoint_bottom = row_top - 1, row_top + _ROW_HEIGHT + 1
        elements.append(
            f'<line x1="{bx:.1f}" y1="{breakpoint_top:.1f}" x2="{bx:.1f}" '
            f'y2="{breakpoint_bottom:.1f}" stroke="{BREAKPOINT_COLOR}" stroke-width="1.6"/>'
        )

        label = _row_label(group, domain_label)
        elements.append(
            f'<text x="{_LABEL_LEFT}" y="{row_mid + 3.5:.1f}" font-family="sans-serif" '
            f'font-size="9.5">{label}</text>'
        )
        y += _ROW_HEIGHT + _ROW_GAP

    elements.extend(
        render_position_axis_svg(y + 4, protein_length, scale, exon_boundaries=exon_boundaries)
    )
    y += 24 + _EXON_TICK_LABEL_HEIGHT
    elements.extend(_legend_svg(y))
    y += _LEGEND_HEIGHT - 6
    if truncated_count:
        elements.append(
            f'<text x="{_AXIS_LEFT}" y="{y + 12:.1f}" font-family="sans-serif" font-size="9" '
            f'fill="#555">Showing the top {len(shown)} of {total_groups} partner/breakpoint '
            "groups by recurrence; see results.tsv for the complete list.</text>"
        )
    elements.append("</svg>")
    return "\n".join(elements)


def _deletion_groups(payload: dict) -> list[dict]:
    groups: dict[tuple[int, int], dict] = {}
    for record in payload.get("intragenic_deletions") or []:
        retained_up_to = record.get("retained_up_to_aa")
        resumed_from = record.get("resumed_from_aa")
        if retained_up_to is None or resumed_from is None:
            continue
        key = (int(retained_up_to), int(resumed_from))
        group = groups.setdefault(
            key,
            {
                "retained_up_to_aa": int(retained_up_to),
                "resumed_from_aa": int(resumed_from),
                "count": 0,
                "n_exons_deleted": record.get("n_exons_deleted"),
                "frame_status": record.get("frame_status"),
            },
        )
        group["count"] += 1
    return list(groups.values())


def render_intragenic_deletion_schematic_svg(
    payload: dict, *, max_rows: int = _MAX_ROWS_DEFAULT
) -> str | None:
    """Render the Figure-4C-style intragenic-deletion schematic: a retained
    N-terminal block, a plain (uncolored) connector line for the deleted
    span, and a resumed domain-colored block, for same-gene
    (``Site1_gene == Site2_gene == target``) deletion-style SV records.

    Returns ``None`` when this run has no such records (most genes/runs
    won't) or no ``gene_track`` to build an axis from -- this is never
    fabricated when the underlying data doesn't have this record type.
    """
    gene_track = payload.get("gene_track")
    if not gene_track or not gene_track.get("protein_length"):
        return None
    groups = _deletion_groups(payload)
    if not groups:
        return None

    protein_length = gene_track["protein_length"]
    domains = gene_track.get("domains") or []
    exon_boundaries = gene_track.get("exon_boundaries_aa") or []

    groups.sort(key=lambda g: (-g["count"], g["retained_up_to_aa"]))
    total_groups = len(groups)
    shown = groups[:max_rows]
    truncated_count = total_groups - len(shown)

    scale = _axis_scale(protein_length)
    n_rows = len(shown)
    plot_height = n_rows * (_ROW_HEIGHT + _ROW_GAP)
    truncation_note_height = 16 if truncated_count else 0
    height = (
        _TOP_MARGIN
        + plot_height
        + _EXON_TICK_LABEL_HEIGHT
        + _LEGEND_HEIGHT
        + truncation_note_height
        + _BOTTOM_MARGIN
    )

    gene_symbol = payload.get("gene_symbol") or "target gene"
    elements = [_svg_open(_SVG_WIDTH, height)]
    elements.append(
        f'<text x="{_AXIS_LEFT}" y="24" font-family="sans-serif" font-size="15">'
        f"{gene_symbol} intragenic-deletion schematic</text>"
    )
    elements.append(
        f'<text x="{_AXIS_LEFT}" y="40" font-family="sans-serif" font-size="10" fill="#555">'
        "retained N-terminal block → deleted span (plain connector) → resumed "
        "C-terminal block</text>"
    )

    y = _TOP_MARGIN
    for group in shown:
        row_top = y
        row_mid = y + _ROW_HEIGHT / 2
        retained_up_to = _clip(group["retained_up_to_aa"], 0, protein_length)
        resumed_from = _clip(group["resumed_from_aa"], 0, protein_length)

        for interval_start, interval_end in ((0, retained_up_to), (resumed_from, protein_length)):
            ix0 = _AXIS_LEFT + interval_start * scale
            ix1 = _AXIS_LEFT + interval_end * scale
            elements.append(
                f'<rect x="{ix0:.1f}" y="{row_top:.1f}" width="{max(0.5, ix1 - ix0):.1f}" '
                f'height="{_ROW_HEIGHT}" fill="{BACKBONE_COLOR}"/>'
            )
            for seg_start, seg_end, seg_color, seg_name in _domain_color_segments(
                domains, interval_start, interval_end
            ):
                sx0 = _AXIS_LEFT + seg_start * scale
                sx1 = _AXIS_LEFT + seg_end * scale
                seg_width = max(0.5, sx1 - sx0)
                elements.append(
                    f'<rect x="{sx0:.1f}" y="{row_top:.1f}" width="{seg_width:.1f}" '
                    f'height="{_ROW_HEIGHT}" fill="{seg_color}"/>'
                )
                seg_label = _fit_domain_label(seg_name, seg_width)
                if seg_label:
                    elements.append(
                        f'<text x="{(sx0 + sx1) / 2:.1f}" y="{row_mid + 3:.1f}" '
                        'font-family="sans-serif" font-size="7" text-anchor="middle" '
                        'fill="#111111" stroke="white" stroke-width="2" paint-order="stroke">'
                        f"{seg_label}</text>"
                    )
            for tick_aa in _exon_ticks(exon_boundaries, interval_start, interval_end):
                tx = _AXIS_LEFT + tick_aa * scale
                elements.append(
                    f'<line x1="{tx:.1f}" y1="{row_top:.1f}" x2="{tx:.1f}" '
                    f'y2="{row_top + _ROW_HEIGHT:.1f}" stroke="#555555" stroke-width="0.6" '
                    'opacity="0.6"/>'
                )

        connector_x0 = _AXIS_LEFT + retained_up_to * scale
        connector_x1 = _AXIS_LEFT + resumed_from * scale
        elements.append(
            f'<line x1="{connector_x0:.1f}" y1="{row_mid:.1f}" x2="{connector_x1:.1f}" '
            f'y2="{row_mid:.1f}" stroke="{CONNECTOR_COLOR}" stroke-width="1.5" '
            'stroke-dasharray="3,2"/>'
        )

        count = group["count"]
        n_exons = group.get("n_exons_deleted")
        frame_status = group.get("frame_status") or "unknown frame"
        base_label = f"{n_exons}-exon deletion ({frame_status})" if n_exons else "deletion"
        label = base_label if count == 1 else f"{base_label} (x{count})"
        elements.append(
            f'<text x="{_LABEL_LEFT}" y="{row_mid + 3.5:.1f}" font-family="sans-serif" '
            f'font-size="9.5">{label}</text>'
        )
        y += _ROW_HEIGHT + _ROW_GAP

    elements.extend(
        render_position_axis_svg(y + 4, protein_length, scale, exon_boundaries=exon_boundaries)
    )
    y += 24 + _EXON_TICK_LABEL_HEIGHT
    elements.extend(_legend_svg(y))
    y += _LEGEND_HEIGHT - 6
    if truncated_count:
        elements.append(
            f'<text x="{_AXIS_LEFT}" y="{y + 12:.1f}" font-family="sans-serif" font-size="9" '
            f'fill="#555">Showing the top {len(shown)} of {total_groups} distinct deletion '
            "breakpoint groups by recurrence; see results.tsv for the complete list.</text>"
        )
    elements.append("</svg>")
    return "\n".join(elements)
