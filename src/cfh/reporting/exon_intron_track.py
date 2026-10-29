"""Exon/intron view of a gene's fusion breakpoints.

The domain-retention lollipop track places each breakpoint at a protein
residue, which folds an intronic breakpoint onto the edge of the nearest
exon. This view keeps the gene's structure instead: exons as boxes, introns
as gaps between them, and three stacked panels --

1. exon inclusion: the share of fusions that keep each exon whole (and,
   lighter, the share that keep only part of it after a mid-exon break);
2. the gene structure, with the configured key domains drawn over the
   exons they span;
3. breakpoints counted per exon or intron, stems stacked by the event's
   key-domain status.

Everything is derived from ``results.json`` alone (``gene_track`` exon
boundaries plus each event's ``breakpoint_exon``/``is_intronic_breakpoint``
/``target_role``), so it renders for runs written before this view existed.

Exon numbers are transcript ranks. An intronic breakpoint is reported
against the exon it was snapped to: for a 5' gene the last exon it keeps,
for a 3' gene the first (see
:func:`cfh.mapping.genome_nexus_source.map_genomic_breakpoint_to_protein_position`).
So a 5' break at exon N lies in the intron after N and a 3' break at exon N
in the intron before it. A break in sequence outside the coding exons
snaps to the first or last coding exon, and lands in the slot upstream of
the first or downstream of the last.
"""

from __future__ import annotations

import math
from collections import Counter
from dataclasses import dataclass, field

from cfh.reporting.palette import (
    AXIS_COLOR,
    BACKBONE_COLOR,
    LOST_COLOR,
    MUTED_TEXT_COLOR,
    RETAINED_COLOR,
    TEXT_COLOR,
    TRUNCATED_COLOR,
)
from cfh.reporting.svg_utils import escape_xml_text

UNKNOWN_STATUS_COLOR = "#c8c8c8"
_STATUS_ORDER = ("retained", "truncated", "lost", "unknown")
_STATUS_COLORS = {
    "retained": RETAINED_COLOR,
    "truncated": TRUNCATED_COLOR,
    "lost": LOST_COLOR,
    "unknown": UNKNOWN_STATUS_COLOR,
}

_AXIS_LEFT = 60.0
_AXIS_WIDTH = 960.0
_INTRON_WIDTH = 16.0
"""Intron width when its genomic length is unknown (runs written before
``gene_track`` carried ``exons_genomic``, and the slots beyond the first
and last coding exon)."""
_MIN_INTRON_WIDTH = 10.0
_INTRON_SHARE = 0.4
"""Share of the axis given to introns when their lengths are known."""
INTRON_FILL = "#eef1f4"
INTRON_HIT_FILL = "#fde7d6"
_MIN_EXON_WIDTH = 5.0
_INCLUSION_HEIGHT = 60.0
_STEM_MAX = 70.0


@dataclass
class _Slot:
    kind: str  # "exon" | "intron"
    rank: int
    """Exon rank; for an intron, the rank of the exon upstream of it (0 =
    upstream of the first coding exon)."""
    x0: float = 0.0
    x1: float = 0.0
    statuses: Counter = field(default_factory=Counter)
    partners: Counter = field(default_factory=Counter)
    roles: Counter = field(default_factory=Counter)

    @property
    def center(self) -> float:
        return (self.x0 + self.x1) / 2

    @property
    def total(self) -> int:
        return sum(self.statuses.values())


def _status(event: dict) -> str:
    """Same reading as the protein-coordinate lollipop track."""
    status = event.get("domain_status")
    fraction = event.get("domain_retained_fraction")
    if event.get("domain_is_truncated") or status == "disrupted":
        return "truncated"
    if fraction == 0.0 or status == "lost":
        return "lost"
    if fraction == 1.0 or status == "retained":
        return "retained"
    return "unknown"


def _key_domains(payload: dict) -> list[dict]:
    summary = payload.get("summary") or {}
    domains = summary.get("key_domains")
    if domains:
        return [d for d in domains if d.get("start_aa") is not None and d.get("end_aa")]
    start, end = summary.get("domain_start_aa"), summary.get("domain_end_aa")
    if start is None or end is None:
        return []
    return [
        {"name": summary.get("domain_accession") or "key domain", "start_aa": start, "end_aa": end}
    ]


def _slot_key(event: dict, ranks: list[int]) -> tuple[str, int] | None:
    rank = event.get("breakpoint_exon")
    role = event.get("target_role")
    if rank is None or role not in {"five_prime", "three_prime"}:
        return None
    if not event.get("is_intronic_breakpoint"):
        return ("exon", rank)
    if role == "five_prime":
        return ("intron", rank)
    upstream = [r for r in ranks if r < rank]
    return ("intron", upstream[-1] if upstream else 0)


def _kept(event: dict, ranks: list[int]) -> tuple[set[int], set[int]]:
    """(exons kept whole, exons kept in part) by this event's gene."""
    rank = event.get("breakpoint_exon")
    role = event.get("target_role")
    if rank is None or role not in {"five_prime", "three_prime"}:
        return set(), set()
    mid_exon = not event.get("is_intronic_breakpoint")
    if role == "five_prime":
        whole = {r for r in ranks if r < rank or (r == rank and not mid_exon)}
    else:
        whole = {r for r in ranks if r > rank or (r == rank and not mid_exon)}
    partial = {rank} if mid_exon and rank in ranks else set()
    return whole, partial


def _intron_label(rank: int, ranks: list[int]) -> str:
    if rank == 0:
        return f"upstream of exon {ranks[0]}, outside this transcript's coding exons"
    if rank == ranks[-1]:
        return f"downstream of exon {rank}, outside this transcript's coding exons"
    return f"intron {rank}"


def _intron_lengths(gene_track: dict, ranks: list[int]) -> dict[int, int]:
    """Genomic length of the intron after each coding exon, keyed like the
    intron slots (by the upstream exon's rank); empty without coordinates."""
    exons = {e["exon_rank"]: e for e in gene_track.get("exons_genomic") or []}
    lengths = {}
    for upstream, downstream in zip(ranks, ranks[1:]):
        a, b = exons.get(upstream), exons.get(downstream)
        if not a or not b:
            continue
        gap = max(b["start"] - a["end"], a["start"] - b["end"]) - 1
        if gap > 0:
            lengths[upstream] = gap
    return lengths


def _format_length(bases: int) -> str:
    return f"{bases / 1000:.1f} kb" if bases >= 1000 else f"{bases} bp"


def render_exon_intron_track_svg(payload: dict) -> str | None:
    """Render the exon/intron view, or ``None`` without exon boundaries or
    mappable breakpoints."""
    boundaries = sorted(
        (payload.get("gene_track") or {}).get("exon_boundaries_aa") or [],
        key=lambda b: b["exon_rank"],
    )
    events = payload.get("events") or []
    if not boundaries or not events:
        return None
    ranks = [b["exon_rank"] for b in boundaries]
    by_rank = {b["exon_rank"]: b for b in boundaries}

    slot_keys = [_slot_key(event, ranks) for event in events]
    used = {key for key in slot_keys if key is not None}
    if not used:
        return None

    # Layout: an intron slot before every exon (the first only when used),
    # each exon sized by coding length, a trailing slot only when used.
    slots: dict[tuple[str, int], _Slot] = {}
    order: list[tuple[str, int]] = []
    if ("intron", 0) in used:
        order.append(("intron", 0))
    for index, rank in enumerate(ranks):
        if index:
            order.append(("intron", ranks[index - 1]))
        order.append(("exon", rank))
    if ("intron", ranks[-1]) in used:
        order.append(("intron", ranks[-1]))
    intron_lengths = _intron_lengths(payload.get("gene_track") or {}, ranks)
    intron_keys = [rank for kind, rank in order if kind == "intron"]
    sized = [rank for rank in intron_keys if rank in intron_lengths]
    fixed = len(intron_keys) - len(sized)
    intron_widths: dict[int, float] = {rank: _INTRON_WIDTH for rank in intron_keys}
    if sized:
        # Introns span 100 bp to 100+ kb, so width follows log length: long
        # introns read as long without flattening the exons to slivers.
        budget = _AXIS_WIDTH * _INTRON_SHARE - fixed * _INTRON_WIDTH
        weights = {rank: math.log10(max(intron_lengths[rank], 10)) for rank in sized}
        total_weight = sum(weights.values())
        for rank in sized:
            intron_widths[rank] = max(_MIN_INTRON_WIDTH, budget * weights[rank] / total_weight)
    exon_space = max(_AXIS_WIDTH - sum(intron_widths.values()), len(ranks) * _MIN_EXON_WIDTH)
    lengths = {r: max(1, by_rank[r]["end_aa"] - by_rank[r]["start_aa"] + 1) for r in ranks}
    total_length = sum(lengths.values())
    x = _AXIS_LEFT
    for kind, rank in order:
        width = (
            intron_widths[rank]
            if kind == "intron"
            else max(_MIN_EXON_WIDTH, exon_space * lengths[rank] / total_length)
        )
        slots[(kind, rank)] = _Slot(kind, rank, x, x + width)
        x += width
    right_edge = x

    whole_counts: Counter = Counter()
    partial_counts: Counter = Counter()
    for event, key in zip(events, slot_keys, strict=True):
        whole, partial = _kept(event, ranks)
        whole_counts.update(whole)
        partial_counts.update(partial)
        if key is not None and key in slots:
            slots[key].statuses[_status(event)] += 1
            slots[key].roles[event["target_role"]] += 1
            if event.get("partner_gene"):
                slots[key].partners[event["partner_gene"]] += 1
    n_events = len(events)

    elements: list[str] = []
    gene = escape_xml_text(str(payload.get("gene_symbol") or ""))

    # 1. Exon inclusion.
    inclusion_top = 30.0
    inclusion_base = inclusion_top + _INCLUSION_HEIGHT
    elements.append(
        f'<text x="{_AXIS_LEFT:.1f}" y="{inclusion_top - 12:.1f}" font-family="sans-serif" '
        f'font-size="11" fill="{TEXT_COLOR}">Share of {n_events} {gene} fusions keeping each exon '
        "(darker: whole exon; lighter: part, after a mid-exon break)</text>"
    )
    for fraction in (0.0, 0.5, 1.0):
        y = inclusion_base - fraction * _INCLUSION_HEIGHT
        elements.append(
            f'<line x1="{_AXIS_LEFT - 4:.1f}" y1="{y:.1f}" x2="{right_edge:.1f}" y2="{y:.1f}" '
            f'stroke="{BACKBONE_COLOR}" stroke-width="0.8"/>'
            f'<text x="{_AXIS_LEFT - 8:.1f}" y="{y + 3:.1f}" font-family="sans-serif" '
            f'font-size="9" text-anchor="end" fill="{MUTED_TEXT_COLOR}">{fraction:.0%}</text>'
        )
    for rank in ranks:
        slot = slots[("exon", rank)]
        whole_share = whole_counts[rank] / n_events
        partial_share = partial_counts[rank] / n_events
        whole_height = whole_share * _INCLUSION_HEIGHT
        partial_height = partial_share * _INCLUSION_HEIGHT
        title = escape_xml_text(
            f"exon {rank}: kept whole by {whole_counts[rank]} of {n_events} fusions "
            f"({whole_share:.0%})"
            + (f", in part by {partial_counts[rank]}" if partial_counts[rank] else "")
        )
        elements.append(
            f"<g><title>{title}</title>"
            f'<rect x="{slot.x0 + 0.5:.1f}" y="{inclusion_base - whole_height:.1f}" '
            f'width="{max(1.0, slot.x1 - slot.x0 - 1):.1f}" height="{whole_height:.1f}" '
            f'fill="{RETAINED_COLOR}"/>'
            f'<rect x="{slot.x0 + 0.5:.1f}" '
            f'y="{inclusion_base - whole_height - partial_height:.1f}" '
            f'width="{max(1.0, slot.x1 - slot.x0 - 1):.1f}" height="{partial_height:.1f}" '
            f'fill="{RETAINED_COLOR}" opacity="0.35"/></g>'
        )

    # 2. Gene structure with key domains.
    structure_y = inclusion_base + 34.0
    exon_height = 14.0
    elements.append(
        f'<line x1="{_AXIS_LEFT:.1f}" y1="{structure_y:.1f}" x2="{right_edge:.1f}" '
        f'y2="{structure_y:.1f}" stroke="{AXIS_COLOR}" stroke-width="1"/>'
    )
    intron_height = exon_height * 0.55
    for kind, rank in order:
        if kind != "intron":
            continue
        slot = slots[(kind, rank)]
        length = intron_lengths.get(rank)
        title = escape_xml_text(
            _intron_label(rank, ranks)
            + (f", {_format_length(length)}" if length else "")
            + (f": {slot.total} breakpoint(s)" if slot.total else ": no breakpoints")
        )
        elements.append(
            f"<g><title>{title}</title>"
            f'<rect x="{slot.x0 + 1:.1f}" y="{structure_y - intron_height / 2:.1f}" '
            f'width="{max(1.0, slot.x1 - slot.x0 - 2):.1f}" height="{intron_height:.1f}" '
            f'fill="{INTRON_HIT_FILL if slot.total else INTRON_FILL}" stroke="{AXIS_COLOR}" '
            f'stroke-width="0.5" stroke-dasharray="2 1.5"/></g>'
        )
    for rank in ranks:
        slot = slots[("exon", rank)]
        elements.append(
            f'<rect x="{slot.x0:.1f}" y="{structure_y - exon_height / 2:.1f}" '
            f'width="{slot.x1 - slot.x0:.1f}" height="{exon_height:.1f}" fill="{BACKBONE_COLOR}" '
            f'stroke="{AXIS_COLOR}" stroke-width="0.6"><title>exon {rank} '
            f"(aa {by_rank[rank]['start_aa']}-{by_rank[rank]['end_aa']})</title></rect>"
        )
    for domain in _key_domains(payload):
        spans = []
        for rank in ranks:
            b, slot = by_rank[rank], slots[("exon", rank)]
            lo, hi = max(b["start_aa"], domain["start_aa"]), min(b["end_aa"], domain["end_aa"])
            if lo > hi:
                continue
            scale = (slot.x1 - slot.x0) / lengths[rank]
            spans.append(
                (slot.x0 + (lo - b["start_aa"]) * scale, slot.x0 + (hi - b["start_aa"] + 1) * scale)
            )
        if not spans:
            continue
        name = escape_xml_text(str(domain.get("name") or "key domain"))
        for x0, x1 in spans:
            elements.append(
                f'<rect x="{x0:.1f}" y="{structure_y - exon_height / 2:.1f}" '
                f'width="{max(1.0, x1 - x0):.1f}" height="{exon_height:.1f}" '
                f'fill="#62b36f" opacity="0.75"><title>{name} '
                f"(aa {domain['start_aa']}-{domain['end_aa']})</title></rect>"
            )
        elements.append(
            f'<text x="{spans[0][0]:.1f}" y="{structure_y - exon_height / 2 - 15:.1f}" '
            f'font-family="sans-serif" font-size="9.5" fill="{TEXT_COLOR}">{name}</text>'
        )
    label_every = max(1, round(len(ranks) / 15))
    for index, rank in enumerate(ranks):
        if index % label_every and rank != ranks[-1]:
            continue
        slot = slots[("exon", rank)]
        elements.append(
            f'<text x="{slot.center:.1f}" y="{structure_y + exon_height / 2 + 11:.1f}" '
            f'font-family="sans-serif" font-size="9" text-anchor="middle" '
            f'fill="{MUTED_TEXT_COLOR}">{rank}</text>'
        )
    elements.append(
        f'<text x="{_AXIS_LEFT - 8:.1f}" y="{structure_y + 3:.1f}" font-family="sans-serif" '
        f'font-size="9" text-anchor="end" fill="{MUTED_TEXT_COLOR}">exon</text>'
    )
    # Name the introns that hold breakpoints, above the structure line.
    busiest = max((slot.total for slot in slots.values() if slot.kind == "intron"), default=0)
    for kind, rank in order:
        slot = slots[(kind, rank)]
        if kind != "intron" or not slot.total or slot.total < max(2, busiest * 0.05):
            continue
        label = "up" if rank == 0 else "down" if rank == ranks[-1] else f"i{rank}"
        elements.append(
            f'<text x="{slot.center:.1f}" y="{structure_y - exon_height / 2 - 4:.1f}" '
            f'font-family="sans-serif" font-size="9" font-weight="600" text-anchor="middle" '
            f'fill="{TEXT_COLOR}">{label}</text>'
        )

    # 3. Breakpoints per exon or intron.
    stems_top = structure_y + exon_height / 2 + 22.0
    peak = max(slot.total for slot in slots.values()) or 1
    elements.append(
        f'<text x="{_AXIS_LEFT - 8:.1f}" y="{stems_top + 8:.1f}" font-family="sans-serif" '
        f'font-size="9" text-anchor="end" fill="{MUTED_TEXT_COLOR}">breaks</text>'
        f'<text x="{_AXIS_LEFT - 8:.1f}" y="{stems_top + 19:.1f}" font-family="sans-serif" '
        f'font-size="8" text-anchor="end" fill="{MUTED_TEXT_COLOR}">(√ scale)</text>'
    )
    for key in order:
        slot = slots[key]
        if not slot.total:
            continue
        label = (
            f"exon {slot.rank} (mid-exon)"
            if slot.kind == "exon"
            else _intron_label(slot.rank, ranks)
        )
        breakdown = ", ".join(
            f"{status} {slot.statuses[status]}" for status in _STATUS_ORDER if slot.statuses[status]
        )
        partners = ", ".join(f"{p} {n}" for p, n in slot.partners.most_common(3))
        roles = ", ".join(
            f"{label} {slot.roles[role]}"
            for role, label in (("five_prime", "as 5' gene"), ("three_prime", "as 3' gene"))
            if slot.roles[role]
        )
        title = escape_xml_text(
            f"{label}: {slot.total} breakpoint(s) ({breakdown}; {roles})"
            + (f"; top partners {partners}" if partners else "")
        )
        # Square-root scale so a handful of breaks stays visible beside a
        # dominant junction (ALK intron 19 holds 236 of 271).
        height = max(4.0, _STEM_MAX * math.sqrt(slot.total / peak))
        y = stems_top
        parts = []
        for status in _STATUS_ORDER:
            count = slot.statuses[status]
            if not count:
                continue
            segment = height * count / slot.total
            parts.append(
                f'<rect x="{slot.center - 2:.1f}" y="{y:.1f}" width="4" height="{segment:.2f}" '
                f'fill="{_STATUS_COLORS[status]}"/>'
            )
            y += segment
        head_color = "#ffffff" if slot.kind == "intron" else AXIS_COLOR
        parts.append(
            f'<circle cx="{slot.center:.1f}" cy="{y + 3:.1f}" r="3.2" fill="{head_color}" '
            f'stroke="{AXIS_COLOR}" stroke-width="1"/>'
        )
        if slot.total >= max(2, peak * 0.05):
            parts.append(
                f'<text x="{slot.center:.1f}" y="{y + 16:.1f}" font-family="sans-serif" '
                f'font-size="8.5" text-anchor="middle" fill="{TEXT_COLOR}">{slot.total}</text>'
            )
        elements.append(
            f'<g data-slot="{slot.kind}-{slot.rank}"><title>{title}</title>{"".join(parts)}</g>'
        )

    legend_y = stems_top + _STEM_MAX + 34.0
    legend = []
    lx = _AXIS_LEFT
    for status in _STATUS_ORDER:
        legend.append(
            f'<rect x="{lx:.1f}" y="{legend_y - 8:.1f}" width="10" height="10" '
            f'fill="{_STATUS_COLORS[status]}"/>'
            f'<text x="{lx + 14:.1f}" y="{legend_y:.1f}" font-family="sans-serif" font-size="10" '
            f'fill="{TEXT_COLOR}">key domain {status}</text>'
        )
        lx += 130
    legend.append(
        f'<circle cx="{lx + 5:.1f}" cy="{legend_y - 3:.1f}" r="3.2" fill="#ffffff" '
        f'stroke="{AXIS_COLOR}"/><text x="{lx + 14:.1f}" y="{legend_y:.1f}" '
        f'font-family="sans-serif" font-size="10" fill="{TEXT_COLOR}">intron</text>'
        f'<circle cx="{lx + 65:.1f}" cy="{legend_y - 3:.1f}" r="3.2" fill="{AXIS_COLOR}" '
        f'stroke="{AXIS_COLOR}"/><text x="{lx + 74:.1f}" y="{legend_y:.1f}" '
        f'font-family="sans-serif" font-size="10" fill="{TEXT_COLOR}">mid-exon</text>'
    )
    # Second legend row: what the gene-structure shapes mean.
    row2 = legend_y + 18
    sizing = (
        "intron width follows log genomic length"
        if intron_lengths
        else "intron lengths not recorded in this run; drawn at a fixed width"
    )
    legend.append(
        f'<rect x="{_AXIS_LEFT:.1f}" y="{row2 - 9:.1f}" width="16" height="11" '
        f'fill="{BACKBONE_COLOR}" stroke="{AXIS_COLOR}" stroke-width="0.6"/>'
        f'<text x="{_AXIS_LEFT + 21:.1f}" y="{row2:.1f}" font-family="sans-serif" '
        f'font-size="10" fill="{TEXT_COLOR}">exon (width: coding length)</text>'
        f'<rect x="{_AXIS_LEFT + 190:.1f}" y="{row2 - 7:.1f}" width="16" height="7" '
        f'fill="{INTRON_FILL}" stroke="{AXIS_COLOR}" stroke-width="0.5" stroke-dasharray="2 1.5"/>'
        f'<rect x="{_AXIS_LEFT + 210:.1f}" y="{row2 - 7:.1f}" width="16" height="7" '
        f'fill="{INTRON_HIT_FILL}" stroke="{AXIS_COLOR}" stroke-width="0.5" '
        f'stroke-dasharray="2 1.5"/>'
        f'<text x="{_AXIS_LEFT + 231:.1f}" y="{row2:.1f}" font-family="sans-serif" '
        f'font-size="10" fill="{TEXT_COLOR}">intron, without / with breakpoints ({sizing})</text>'
    )
    elements.extend(legend)

    height = legend_y + 30
    width = max(right_edge + 20, _AXIS_LEFT + _AXIS_WIDTH + 20)
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width:.0f}" height="{height:.0f}" '
        f'viewBox="0 0 {width:.0f} {height:.0f}" role="img" '
        f'aria-label="{gene} breakpoints by exon and intron">'
        f'<rect width="100%" height="100%" fill="#ffffff"/>' + "".join(elements) + "</svg>"
    )
