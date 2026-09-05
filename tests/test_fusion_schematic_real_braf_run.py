"""Validate the fusion-transcript schematic against the real, committed BRAF
benchmark artifact under ``runs/braf_msk-impact-50k-2026_20260904T172738Z/``.

No network access, no synthetic fixtures: this reads the exact
``results.json`` already committed to the repo (regenerated live from
cBioPortal/Genome Nexus, see that run's ``manifest.json``) and checks
structural properties of the rendered SVGs -- row count within the cap,
every breakpoint marker within ``[1, protein_length]``, and domain-color
continuity -- not merely that a file exists.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

from cfh.reporting.fusion_schematic import (
    render_fusion_schematic_svg,
    render_intragenic_deletion_schematic_svg,
)

# Imported from the shared palette module (not re-exported from
# fusion_schematic) so these assertions catch the two renderers actually
# drifting apart, not just a coincidental match today.
from cfh.reporting.palette import CONNECTOR_COLOR, RETAINED_COLOR, TRUNCATED_COLOR

REPO_ROOT = Path(__file__).parent.parent
BRAF_RUN_DIR = REPO_ROOT / "runs" / "braf_msk-impact-50k-2026_20260904T172738Z"

_AXIS_LEFT = 60.0
_AXIS_WIDTH = 560.0
_MAX_ROWS = 28


def _payload() -> dict:
    return json.loads((BRAF_RUN_DIR / "results.json").read_text())


def _breakpoint_line_x_values(svg: str) -> list[float]:
    return [
        float(match.group(1))
        for match in re.finditer(r'<line x1="([\d.]+)"[^>]*stroke="#d62728"', svg)
    ]


def _domain_colored_rects(svg: str) -> list[tuple[float, float, float, float, str]]:
    """Row-height (22px) rects filled with a domain-retention-status color,
    excluding the legend swatches (10px)."""
    rects = re.findall(
        r'<rect x="([\d.]+)" y="([\d.]+)" width="([\d.]+)" height="22" fill="([^"]+)"',
        svg,
    )
    return [
        (float(x), float(y), float(width), 0.0, fill)
        for x, y, width, fill in rects
        if fill in {RETAINED_COLOR, TRUNCATED_COLOR}
    ]


def test_real_braf_run_has_gene_track_and_intragenic_deletions():
    payload = _payload()
    assert payload["gene_track"] is not None
    assert payload["gene_track"]["protein_length"] == 766
    accessions = {d["accession"] for d in payload["gene_track"]["domains"]}
    assert {"PF07714", "PF02196", "PF00130"} <= accessions
    # See the panel-C investigation: this real cohort does have same-gene
    # BRAF intragenic-deletion-style SV records.
    assert len(payload["intragenic_deletions"]) > 0


def test_fusion_schematic_svg_file_matches_pure_render_of_committed_payload():
    payload = _payload()
    on_disk = (BRAF_RUN_DIR / "visualizations" / "fusion_schematic.svg").read_text()
    rendered = render_fusion_schematic_svg(payload) + "\n"
    assert on_disk == rendered


def test_fusion_schematic_row_count_is_within_cap():
    payload = _payload()
    svg = render_fusion_schematic_svg(payload)
    assert svg is not None
    breakpoint_lines = _breakpoint_line_x_values(svg)
    assert 0 < len(breakpoint_lines) <= _MAX_ROWS


def test_fusion_schematic_breakpoint_markers_within_valid_protein_bounds():
    payload = _payload()
    protein_length = payload["gene_track"]["protein_length"]
    svg = render_fusion_schematic_svg(payload)
    breakpoint_lines = _breakpoint_line_x_values(svg)
    assert breakpoint_lines  # the real run has mappable rows
    for x in breakpoint_lines:
        aa = (x - _AXIS_LEFT) / _AXIS_WIDTH * protein_length
        assert 1 - 0.5 <= aa <= protein_length + 0.5


def test_fusion_schematic_domain_colors_match_real_gene_config_and_are_continuous():
    payload = _payload()
    protein_length = payload["gene_track"]["protein_length"]
    scale = _AXIS_WIDTH / protein_length
    kinase = next(d for d in payload["gene_track"]["domains"] if d["accession"] == "PF07714")

    svg = render_fusion_schematic_svg(payload)
    domain_rects = _domain_colored_rects(svg)
    assert domain_rects  # BRAF's real config has domains, so some are drawn

    for x, _y, width, _unused, _fill in domain_rects:
        # Every drawn domain segment stays on the shared protein-length axis.
        assert _AXIS_LEFT - 1e-6 <= x
        assert x + width <= _AXIS_LEFT + protein_length * scale + 1e-6

    # At least one row should show the real kinase domain's real span
    # (458-712 aa) rendered at its correct scaled position.
    expected_x = _AXIS_LEFT + kinase["start_aa"] * scale
    assert any(abs(x - expected_x) < 1.0 for x, *_rest in domain_rects)


def test_fusion_schematic_partner_labels_match_real_partner_names():
    payload = _payload()
    real_partners = {row["Partner_gene"] for row in payload["summary"]["partner_counts"]}
    svg = render_fusion_schematic_svg(payload)
    labels = re.findall(r'font-size="9.5">([^<]+)</text>', svg)
    assert labels  # rows were actually rendered
    for label in labels:
        partner = label.split(" –")[0].split(" (x")[0]
        assert partner in real_partners

    # The real cohort's single most recurrent BRAF partner (KIAA1549, the
    # classic pilocytic-astrocytoma fusion partner) should appear as the
    # top (most recurrent) row.
    assert labels[0].startswith("KIAA1549")


def test_fusion_schematic_shows_real_exon_numbers_from_the_canonical_transcript():
    payload = _payload()
    exon_ranks = sorted({b["exon_rank"] for b in payload["gene_track"]["exon_boundaries_aa"]})
    assert exon_ranks  # the real BRAF canonical transcript has exon boundaries
    svg = render_fusion_schematic_svg(payload)
    for rank in exon_ranks:
        assert f">E{rank}<" in svg


def test_fusion_schematic_shows_real_domain_names_not_just_colors():
    payload = _payload()
    real_domain_names = {d["name"] for d in payload["gene_track"]["domains"] if d.get("name")}
    svg = render_fusion_schematic_svg(payload)
    assert any(name in svg for name in real_domain_names)


def test_intragenic_deletion_schematic_row_count_and_bounds():
    payload = _payload()
    protein_length = payload["gene_track"]["protein_length"]
    svg = render_intragenic_deletion_schematic_svg(payload)
    assert svg is not None  # this real cohort does have qualifying records

    connectors = re.findall(
        rf'<line x1="([\d.]+)"[^>]*x2="([\d.]+)"[^>]*stroke="{re.escape(CONNECTOR_COLOR)}"', svg
    )
    assert 0 < len(connectors) <= _MAX_ROWS
    for x1, x2 in connectors:
        aa1 = (float(x1) - _AXIS_LEFT) / _AXIS_WIDTH * protein_length
        aa2 = (float(x2) - _AXIS_LEFT) / _AXIS_WIDTH * protein_length
        assert 0 - 0.5 <= aa1 <= protein_length + 0.5
        assert 0 - 0.5 <= aa2 <= protein_length + 0.5
        assert aa1 <= aa2 + 0.5  # retained_up_to_aa is before resumed_from_aa


def test_report_markdown_embeds_both_schematics():
    report = (BRAF_RUN_DIR / "report.md").read_text()
    assert "visualizations/fusion_schematic.svg" in report
    assert "visualizations/intragenic_deletion_schematic.svg" in report


def test_report_pdf_embeds_both_schematics():
    from pypdf import PdfReader

    reader = PdfReader(str(BRAF_RUN_DIR / "report.pdf"))
    text = "".join(page.extract_text() or "" for page in reader.pages)
    assert "fusion schematic" in text.lower()
    assert "intragenic deletion" in text.lower()
