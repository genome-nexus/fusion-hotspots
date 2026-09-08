"""Tests for the gene-agnostic fusion-transcript schematic renderer.

These are pure structural tests over synthetic ``results.json``-shaped
payloads (no network, no gene-specific literals baked into assertions
beyond what the payload itself declares) -- see
:mod:`cfh.reporting.fusion_schematic`.
"""

from __future__ import annotations

import re

import pytest

from cfh.reporting.fusion_schematic import (
    BACKBONE_COLOR,
    _domain_color_segments,
    _fit_domain_label,
    _fusion_groups,
    exon_boundary_ticks_svg,
    partner_color,
    render_fusion_schematic_svg,
    render_intragenic_deletion_schematic_svg,
    render_position_axis_svg,
)

# Imported from the shared palette module (not re-exported from
# fusion_schematic) so these assertions catch the two renderers actually
# drifting apart, not just a coincidental match today.
from cfh.reporting.palette import (
    BREAKPOINT_COLOR,
    PARTNER_COLOR,
    RETAINED_COLOR,
    TRUNCATED_COLOR,
)

_KINASE = {
    "name": "Protein kinase domain",
    "accession": "PF07714",
    "start_aa": 458,
    "end_aa": 712,
}
_RAS_BINDING = {
    "name": "RAS-binding domain",
    "accession": "PF02196",
    "start_aa": 156,
    "end_aa": 227,
}
_CYS_RICH = {
    "name": "Cysteine-rich domain",
    "accession": "PF00130",
    "start_aa": 235,
    "end_aa": 280,
}
_DOMAINS = [_KINASE, _RAS_BINDING, _CYS_RICH]
_PROTEIN_LENGTH = 766


def _gene_track(domains=None, exon_boundaries=None):
    return {
        "protein_length": _PROTEIN_LENGTH,
        "domains": _DOMAINS if domains is None else domains,
        "exon_boundaries_aa": exon_boundaries or [],
    }


def _event(partner, breakpoint_aa, role, *, status="retained", sample_id=None):
    return {
        "partner_gene": partner,
        "breakpoint_protein_position": breakpoint_aa,
        "target_role": role,
        "domain_status": status,
        "sample_id": sample_id or f"{partner}-{breakpoint_aa}",
    }


_UNSET = object()


def _payload(events, *, gene_track=_UNSET, intragenic_deletions=None, domain_accession="PF07714"):
    return {
        "gene_symbol": "BRAF",
        "summary": {"domain_accession": domain_accession},
        "gene_track": _gene_track() if gene_track is _UNSET else gene_track,
        "events": events,
        "intragenic_deletions": intragenic_deletions or [],
    }


def _rects(svg):
    return re.findall(
        r'<rect x="([\d.]+)" y="([\d.]+)" width="([\d.]+)" height="([\d.]+)" fill="([^"]+)"',
        svg,
    )


def _breakpoint_lines_x(svg):
    return [
        float(m.group(1))
        for m in re.finditer(
            rf'<line x1="([\d.]+)"[^>]*stroke="{re.escape(BREAKPOINT_COLOR)}"', svg
        )
    ]


# --- honest degradation -----------------------------------------------------


def test_returns_none_without_gene_track():
    payload = _payload([_event("AGK", 380, "three_prime")], gene_track=None)
    assert render_fusion_schematic_svg(payload) is None


def test_returns_none_without_protein_length():
    payload = _payload([_event("AGK", 380, "three_prime")], gene_track={"protein_length": None})
    assert render_fusion_schematic_svg(payload) is None


def test_returns_none_with_no_mappable_events():
    payload = _payload([])
    assert render_fusion_schematic_svg(payload) is None


def test_deletion_schematic_returns_none_with_no_records():
    payload = _payload([_event("AGK", 380, "three_prime")], intragenic_deletions=[])
    assert render_intragenic_deletion_schematic_svg(payload) is None


def test_deletion_schematic_returns_none_without_gene_track():
    payload = _payload(
        [],
        gene_track=None,
        intragenic_deletions=[
            {
                "retained_up_to_aa": 150,
                "resumed_from_aa": 450,
                "n_exons_deleted": 5,
                "frame_status": "in-frame",
            }
        ],
    )
    assert render_intragenic_deletion_schematic_svg(payload) is None


# --- grouping / PARTNER (xN) labeling ---------------------------------------


def test_identical_partner_breakpoint_pairs_collapse_with_count():
    events = [
        _event("AGK", 380, "three_prime", sample_id="S1"),
        _event("AGK", 380, "three_prime", sample_id="S2"),
        _event("AGK", 380, "three_prime", sample_id="S3"),
    ]
    groups = _fusion_groups(_payload(events))
    assert len(groups) == 1
    assert groups[0]["count"] == 3


def test_same_partner_different_breakpoints_stay_separate_rows():
    events = [
        _event("AGK", 380, "three_prime"),
        _event("AGK", 500, "three_prime"),
    ]
    groups = _fusion_groups(_payload(events))
    assert len(groups) == 2


def test_multiplicity_label_rendered_for_repeated_group():
    events = [_event("AGK", 380, "three_prime", sample_id=f"S{i}") for i in range(5)]
    svg = render_fusion_schematic_svg(_payload(events))
    assert "AGK (x5)" in svg


def test_single_event_label_has_no_multiplicity_suffix():
    svg = render_fusion_schematic_svg(_payload([_event("AGK", 380, "three_prime")]))
    assert "AGK (x1)" not in svg
    assert re.search(r">AGK\b", svg) or "AGK –" in svg


# --- structural properties: row cap, breakpoint bounds, domain continuity --


def test_row_count_is_capped_and_truncation_is_noted():
    events = [_event(f"PARTNER{i}", 300 + i, "three_prime") for i in range(40)]
    svg = render_fusion_schematic_svg(_payload(events), max_rows=28)
    breakpoint_lines = _breakpoint_lines_x(svg)
    assert len(breakpoint_lines) == 28
    assert "Showing the top 28 of 40" in svg


def test_row_count_not_truncated_below_cap():
    events = [_event(f"PARTNER{i}", 300 + i, "three_prime") for i in range(5)]
    svg = render_fusion_schematic_svg(_payload(events), max_rows=28)
    assert len(_breakpoint_lines_x(svg)) == 5
    assert "Showing the top" not in svg


@pytest.mark.parametrize("role", ["three_prime", "five_prime"])
def test_breakpoint_markers_fall_within_valid_protein_bounds(role):
    events = [
        _event("AGK", 1, role),
        _event("CUL1", _PROTEIN_LENGTH, role),
        _event("SND1", 380, role),
    ]
    svg = render_fusion_schematic_svg(_payload(events))
    axis_left = 60.0
    axis_width = 560.0
    for x in _breakpoint_lines_x(svg):
        aa = (x - axis_left) / axis_width * _PROTEIN_LENGTH
        # +/-0.5aa tolerance absorbs the SVG coordinate's 1-decimal rounding.
        assert 1 - 0.5 <= aa <= _PROTEIN_LENGTH + 0.5


def test_domain_color_segments_stay_within_requested_interval_and_are_non_overlapping():
    segments = _domain_color_segments(_DOMAINS, 0, 380)
    for start, end, color, _name in segments:
        assert 0 <= start <= end <= 380
        assert color in {RETAINED_COLOR, TRUNCATED_COLOR}
    ordered = sorted(segments)
    for (start_a, end_a, *_rest), (start_b, *_rest2) in zip(ordered, ordered[1:]):
        assert end_a <= start_b + 1e-9  # no two domain segments overlap


def test_domain_fully_outside_retained_span_is_not_drawn():
    # Kinase domain (458-712) has no overlap with [0, 380].
    segments = _domain_color_segments(_DOMAINS, 0, 380)
    assert all(name != _KINASE["name"] for *_rest, name in segments)


def test_domain_straddling_span_boundary_is_truncated_color():
    # RAS-binding (156-227) is fully inside [0, 380] -> retained.
    # A span of [0, 200] truncates it (156-200 retained, 200-227 lost).
    segments = _domain_color_segments(_DOMAINS, 0, 200)
    ras = next(s for s in segments if s[3] == _RAS_BINDING["name"])
    assert ras[2] == TRUNCATED_COLOR
    assert ras[0] == 156
    assert ras[1] == 200


def test_kinase_retained_or_lost_matches_role_and_breakpoint():
    # three_prime role retains [breakpoint, protein_length]; kinase domain
    # (458-712) should be drawn (retained) when breakpoint <= 458.
    svg_retained = render_fusion_schematic_svg(
        _payload([_event("AGK", 400, "three_prime", status="retained")])
    )
    assert RETAINED_COLOR in svg_retained
    assert "retained" in svg_retained

    # five_prime role retains [0, breakpoint]; kinase domain should be
    # entirely absent (lost) when breakpoint < 458.
    svg_lost = render_fusion_schematic_svg(
        _payload([_event("AGK", 400, "five_prime", status="lost")])
    )
    assert "lost" in svg_lost


def test_svg_has_valid_dimensions_and_is_well_formed_xml():
    svg = render_fusion_schematic_svg(_payload([_event("AGK", 380, "three_prime")]))
    match = re.match(r'<svg xmlns="[^"]+" width="(\d+)" height="(\d+)"', svg)
    assert match is not None
    width, height = int(match.group(1)), int(match.group(2))
    assert width > 0 and height > 0
    assert svg.strip().endswith("</svg>")


@pytest.mark.parametrize(
    ("role", "five_prime_block", "three_prime_block"),
    [
        ("three_prime", "partner", "target"),
        ("five_prime", "target", "partner"),
    ],
)
def test_fusion_end_labels_follow_the_role_dependent_blocks(
    role, five_prime_block, three_prime_block
):
    svg = render_fusion_schematic_svg(_payload([_event("AGK", 380, role)]))
    assert f'data-fusion-end="5-prime" data-block="{five_prime_block}">5\'' in svg
    assert f'data-fusion-end="3-prime" data-block="{three_prime_block}">3\'' in svg


def test_shared_position_axis_labels_coding_sequence_ends():
    svg = "\n".join(render_position_axis_svg(100, _PROTEIN_LENGTH, 1.0))
    assert 'data-axis-end="5-prime" text-anchor="start">5\' / N-terminus<' in svg
    assert 'data-axis-end="3-prime" text-anchor="end">3\' / C-terminus<' in svg


def test_header_defers_to_each_rows_orientation_labels():
    svg = render_fusion_schematic_svg(_payload([_event("AGK", 380, "five_prime")]))
    assert "Fusion order varies by row; follow each row's 5' and 3' labels" in svg
    assert "partner block → breakpoint" not in svg


# --- partner coloring --------------------------------------------------------


def test_partner_color_is_deterministic_and_a_valid_hex_color():
    color1 = partner_color("AGK")
    color2 = partner_color("AGK")
    assert color1 == color2
    assert re.fullmatch(r"#[0-9a-f]{6}", color1)


def test_every_partner_gets_the_same_fixed_color():
    """Partner-gene block color is decorative chrome, not a per-partner
    identifier -- every partner must render in the one shared
    ``PARTNER_COLOR``, not a color derived from its name."""
    colors = {partner_color(name) for name in ["AGK", "CUL1", "SND1", "KIAA1549", "TRIM24"]}
    assert colors == {PARTNER_COLOR}


# --- intragenic deletion (panel-C style) schematic --------------------------


def test_deletion_schematic_draws_plain_connector_for_deleted_span():
    payload = _payload(
        [],
        intragenic_deletions=[
            {
                "retained_up_to_aa": 150,
                "resumed_from_aa": 450,
                "n_exons_deleted": 5,
                "frame_status": "in-frame",
            }
        ],
    )
    svg = render_intragenic_deletion_schematic_svg(payload)
    assert svg is not None
    assert "stroke-dasharray" in svg  # the plain/uncolored connector line
    assert "5-exon deletion (in-frame)" in svg


def test_deletion_schematic_collapses_identical_records_with_count():
    record = {
        "retained_up_to_aa": 150,
        "resumed_from_aa": 450,
        "n_exons_deleted": 5,
        "frame_status": "in-frame",
    }
    payload = _payload([], intragenic_deletions=[dict(record), dict(record), dict(record)])
    svg = render_intragenic_deletion_schematic_svg(payload)
    assert "(x3)" in svg


def test_deletion_schematic_domain_wholly_inside_deleted_span_is_not_drawn():
    # Cysteine-rich (235-280) sits entirely inside the deleted span [200, 400].
    payload = _payload(
        [],
        intragenic_deletions=[
            {
                "retained_up_to_aa": 200,
                "resumed_from_aa": 400,
                "n_exons_deleted": 3,
                "frame_status": "out-of-frame",
            }
        ],
    )
    svg = render_intragenic_deletion_schematic_svg(payload)
    rects = _rects(svg)
    # Row rects are drawn at the row height (22px); the legend swatches use
    # a different height and share the same fill colors, so exclude them.
    domain_colored = [
        r for r in rects if r[4] in {RETAINED_COLOR, TRUNCATED_COLOR} and r[3] == "22"
    ]
    assert domain_colored  # sanity: the kinase domain in the resumed block is drawn
    for x, _y, width, _height, _fill in domain_colored:
        x0 = float(x)
        x1 = x0 + float(width)
        # Cysteine-rich domain would map to roughly aa 235-280 on this axis;
        # none of the drawn domain rects should land inside the deleted span.
        axis_left = 60.0
        scale = 560.0 / _PROTEIN_LENGTH
        deleted_x0 = axis_left + 200 * scale
        deleted_x1 = axis_left + 400 * scale
        assert not (deleted_x0 < x0 < deleted_x1 and deleted_x0 < x1 < deleted_x1)


def test_deletion_schematic_row_cap_and_truncation_note():
    records = [
        {
            "retained_up_to_aa": 100 + i,
            "resumed_from_aa": 500 + i,
            "n_exons_deleted": 2,
            "frame_status": "in-frame",
        }
        for i in range(35)
    ]
    payload = _payload([], intragenic_deletions=records)
    svg = render_intragenic_deletion_schematic_svg(payload, max_rows=28)
    assert "Showing the top 28 of 35" in svg


def test_deletion_schematic_backbone_color_present():
    payload = _payload(
        [],
        intragenic_deletions=[
            {
                "retained_up_to_aa": 150,
                "resumed_from_aa": 450,
                "n_exons_deleted": 5,
                "frame_status": "in-frame",
            }
        ],
    )
    svg = render_intragenic_deletion_schematic_svg(payload)
    assert BACKBONE_COLOR in svg


# --- exon and domain-name labels -------------------------------------------

_EXON_BOUNDARIES = [
    {"exon_rank": 1, "start_aa": 1, "end_aa": 120},
    {"exon_rank": 2, "start_aa": 121, "end_aa": 300},
    {"exon_rank": 3, "start_aa": 301, "end_aa": 458},
    {"exon_rank": 4, "start_aa": 459, "end_aa": 712},
    {"exon_rank": 5, "start_aa": 713, "end_aa": _PROTEIN_LENGTH},
]


def test_exon_boundary_ticks_svg_derives_real_exon_numbers_from_gene_track():
    """The exon number comes from gene_track's own ``exon_rank`` field
    (the real ordinal exon from the canonical transcript), not something
    guessed from position order."""
    elements = exon_boundary_ticks_svg(
        _EXON_BOUNDARIES, axis_left=60.0, scale=1.0, y=100.0, max_position=_PROTEIN_LENGTH
    )
    text = "\n".join(elements)
    for rank in (1, 2, 3, 4, 5):
        assert f">E{rank}<" in text


def test_fusion_schematic_shared_axis_carries_exon_number_labels():
    payload = _payload(
        [_event("AGK", 380, "three_prime")],
        gene_track=_gene_track(exon_boundaries=_EXON_BOUNDARIES),
    )
    svg = render_fusion_schematic_svg(payload)
    for rank in (2, 3, 4, 5):
        assert f">E{rank}<" in svg


def test_deletion_schematic_shared_axis_carries_exon_number_labels():
    payload = _payload(
        [],
        gene_track=_gene_track(exon_boundaries=_EXON_BOUNDARIES),
        intragenic_deletions=[
            {
                "retained_up_to_aa": 150,
                "resumed_from_aa": 450,
                "n_exons_deleted": 5,
                "frame_status": "in-frame",
            }
        ],
    )
    svg = render_intragenic_deletion_schematic_svg(payload)
    for rank in (2, 3, 4, 5):
        assert f">E{rank}<" in svg


def test_fit_domain_label_returns_full_name_when_it_fits():
    assert _fit_domain_label("RAS-binding domain", 200) == "RAS-binding domain"


def test_fit_domain_label_truncates_with_ellipsis_when_narrow():
    label = _fit_domain_label("Protein kinase domain", 40)
    assert label is not None
    assert label.endswith("…")
    assert label != "Protein kinase domain"


def test_fit_domain_label_skips_unreadably_narrow_segments():
    assert _fit_domain_label("Protein kinase domain", 5) is None


def test_fusion_schematic_domain_segment_carries_its_real_name():
    # Kinase domain (458-712) is fully inside the retained [380, 766] span.
    payload = _payload([_event("AGK", 380, "three_prime", status="retained")])
    svg = render_fusion_schematic_svg(payload)
    assert "Protein kinase domain" in svg


def test_deletion_schematic_domain_segment_carries_its_real_name():
    # Kinase domain (458-712) sits inside the resumed [400, 766] block.
    payload = _payload(
        [],
        intragenic_deletions=[
            {
                "retained_up_to_aa": 150,
                "resumed_from_aa": 400,
                "n_exons_deleted": 5,
                "frame_status": "in-frame",
            }
        ],
    )
    svg = render_intragenic_deletion_schematic_svg(payload)
    assert "Protein kinase domain" in svg
