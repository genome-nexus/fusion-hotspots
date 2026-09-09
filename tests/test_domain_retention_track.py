"""Tests for the per-event domain-retention track
(``cfh.real_benchmark._domain_track_svg``), which produces the committed
``domain_retention_outliers.svg`` artifact.

These parse the actual rendered SVG text -- real exon labels, real domain
names and amino-acid boundaries, real position-axis numbers -- rather than
just asserting a file/element exists, so a future regression that silently
drops a label (as the pre-fix version of this renderer did) is caught.
"""

from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from types import SimpleNamespace

from cfh.real_benchmark import _domain_track_svg

_EXON_BOUNDARIES = [
    {"exon_rank": 1, "start_aa": 1, "end_aa": 50},
    {"exon_rank": 2, "start_aa": 51, "end_aa": 120},
    {"exon_rank": 3, "start_aa": 121, "end_aa": 300},
    {"exon_rank": 4, "start_aa": 301, "end_aa": 458},
    {"exon_rank": 5, "start_aa": 459, "end_aa": 600},
    {"exon_rank": 6, "start_aa": 601, "end_aa": 712},
    {"exon_rank": 7, "start_aa": 713, "end_aa": 766},
]


_KINASE_DOMAIN = {
    "name": "Protein kinase domain",
    "accession": "PF07714",
    "start_aa": 458,
    "end_aa": 712,
}


def _gene_track(*, domains=None, protein_length=766):
    return {
        "protein_length": protein_length,
        "domains": domains or [_KINASE_DOMAIN],
        "exon_boundaries_aa": _EXON_BOUNDARIES,
    }


def _row(
    event_id,
    position,
    *,
    status="retained",
    fraction=1.0,
    truncated=False,
    sample_id=None,
    partner_gene=None,
):
    return {
        "event_id": event_id,
        "breakpoint_protein_position": position,
        "domain_status": status,
        "domain_retained_fraction": fraction,
        "domain_is_truncated": truncated,
        "sample_id": sample_id,
        "partner_gene": partner_gene,
    }


def _run(*, summary, rows, gene_track=None, gene_symbol="BRAF"):
    return SimpleNamespace(
        gene_symbol=gene_symbol,
        summary=summary,
        rows=rows,
        gene_track=gene_track if gene_track is not None else _gene_track(),
    )


def _text_labels(svg: str) -> list[str]:
    return re.findall(r">([^<]+)</text>", svg)


def test_domain_highlight_is_labeled_with_name_and_boundaries():
    run = _run(
        summary={
            "domain_accession": "PF07714",
            "domain_start_aa": 458,
            "domain_end_aa": 712,
            "key_domains": [_KINASE_DOMAIN],
        },
        rows=[_row("E1", 400)],
    )
    svg = _domain_track_svg(run, set())
    assert "Protein kinase domain (458-712)" in svg


def test_falls_back_to_legacy_single_domain_summary_fields_when_key_domains_missing():
    """A results.json written before the ``key_domains`` summary field
    existed should still get a labeled highlight, by looking the
    accession's pretty name up in ``gene_track["domains"]``."""
    run = _run(
        summary={"domain_accession": "PF07714", "domain_start_aa": 458, "domain_end_aa": 712},
        rows=[_row("E1", 400)],
    )
    svg = _domain_track_svg(run, set())
    assert "Protein kinase domain (458-712)" in svg


def test_multiple_key_domains_each_get_their_own_labeled_segment():
    ras_binding_domain = {
        "name": "RAS-binding domain",
        "accession": "PF02196",
        "start_aa": 156,
        "end_aa": 227,
    }
    run = _run(
        summary={
            "domain_accession": "PF07714",
            "domain_start_aa": 458,
            "domain_end_aa": 712,
            "key_domains": [_KINASE_DOMAIN, ras_binding_domain],
        },
        rows=[_row("E1", 400)],
    )
    svg = _domain_track_svg(run, set())
    assert "Protein kinase domain (458-712)" in svg
    assert "RAS-binding domain (156-227)" in svg
    # Distinct fill colors for the two domain highlight rects.
    fills = set(re.findall(r'<rect[^>]*fill="(#[0-9a-f]{6})"[^>]*opacity="0.55"', svg))
    assert len(fills) == 2


def test_exon_boundary_ticks_show_real_exon_numbers():
    run = _run(
        summary={"domain_accession": "PF07714", "domain_start_aa": 458, "domain_end_aa": 712},
        rows=[_row("E1", 400)],
    )
    svg = _domain_track_svg(run, set())
    labels = _text_labels(svg)
    for expected in ("E1", "E2", "E3", "E4", "E5", "E6", "E7"):
        assert expected in labels


def test_position_axis_includes_transcript_endpoints_and_hundred_step_ticks():
    run = _run(
        summary={"domain_accession": "PF07714", "domain_start_aa": 458, "domain_end_aa": 712},
        rows=[_row("E1", 700)],
        gene_track=_gene_track(protein_length=766),
    )
    svg = _domain_track_svg(run, set())
    labels = _text_labels(svg)
    assert "1" in labels
    assert "766 aa" in labels
    assert "5' / N-terminus" in labels
    assert "3' / C-terminus" in labels
    for hundred in ("100", "200", "300", "400", "500", "600", "700"):
        assert hundred in labels


def test_no_key_domain_still_renders_axis_and_dots_without_a_highlight():
    run = _run(
        summary={"domain_accession": None, "domain_start_aa": None, "domain_end_aa": None},
        rows=[_row("E1", 400)],
    )
    svg = _domain_track_svg(run, set())
    assert svg.strip().endswith("</svg>")
    assert "(458-712)" not in svg  # nothing to derive a highlight from


def test_outlier_events_keep_reference_discrepancy_stroke():
    run = _run(
        summary={"domain_accession": "PF07714", "domain_start_aa": 458, "domain_end_aa": 712},
        rows=[_row("E1", 400), _row("E2", 420)],
    )
    svg = _domain_track_svg(run, {"E1"})
    assert 'stroke="#d62728" stroke-width="1.5"' in svg


# --- <title> tooltips --------------------------------------------------------

_SVG_NS = "{http://www.w3.org/2000/svg}"


def _title_texts(svg: str) -> list[str]:
    """Parse ``svg`` as real XML (raising ``ParseError`` if any embedded
    text was left unescaped) and return every ``<title>`` element's text."""
    root = ET.fromstring(svg)
    return [el.text or "" for el in root.iter(f"{_SVG_NS}title")]


def test_dot_title_carries_event_sample_partner_position_and_status():
    run = _run(
        summary={"domain_accession": "PF07714", "domain_start_aa": 458, "domain_end_aa": 712},
        rows=[_row("E1", 400, status="retained", sample_id="S1", partner_gene="AGK")],
    )
    svg = _domain_track_svg(run, set())
    titles = _title_texts(svg)
    assert any(
        "event E1" in t
        and "sample S1" in t
        and "partner AGK" in t
        and "breakpoint aa 400" in t
        and "domain retained" in t
        for t in titles
    )


def test_dot_title_notes_reference_discrepancy_for_outlier_events():
    run = _run(
        summary={"domain_accession": "PF07714", "domain_start_aa": 458, "domain_end_aa": 712},
        rows=[_row("E1", 400, sample_id="S1")],
    )
    svg = _domain_track_svg(run, {"E1"})
    titles = _title_texts(svg)
    assert any("reference discrepancy" in t for t in titles)
    svg_no_outlier = _domain_track_svg(run, set())
    titles_no_outlier = _title_texts(svg_no_outlier)
    assert not any("reference discrepancy" in t for t in titles_no_outlier)


def test_dot_title_reflects_truncated_and_lost_status():
    run = _run(
        summary={"domain_accession": "PF07714", "domain_start_aa": 458, "domain_end_aa": 712},
        rows=[
            _row("E1", 400, status="disrupted", truncated=True, fraction=0.5, sample_id="S1"),
            _row("E2", 420, status="lost", fraction=0.0, sample_id="S2"),
        ],
    )
    svg = _domain_track_svg(run, set())
    titles = _title_texts(svg)
    assert any("event E1" in t and "domain truncated" in t for t in titles)
    assert any("event E2" in t and "domain lost" in t for t in titles)


def test_dot_title_round_trips_xml_unsafe_sample_id():
    """A sample ID with XML-significant characters must be escaped in the
    tooltip -- otherwise the whole SVG document would not be well-formed
    XML at all. ``_title_texts`` parses the rendered SVG with a real XML
    parser, so an unescaped ``&``/``<`` here raises
    ``xml.etree.ElementTree.ParseError`` rather than just looking wrong.
    """
    unsafe_sample_id = "S&1<2>samp"
    run = _run(
        summary={"domain_accession": "PF07714", "domain_start_aa": 458, "domain_end_aa": 712},
        rows=[_row("E1", 400, sample_id=unsafe_sample_id)],
    )
    svg = _domain_track_svg(run, set())
    titles = _title_texts(svg)
    assert any(unsafe_sample_id in t for t in titles)
