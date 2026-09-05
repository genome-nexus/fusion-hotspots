"""Role-orientation regression coverage from the committed ETV6 cohort scan."""

from __future__ import annotations

import json
import xml.etree.ElementTree as ET
from pathlib import Path

from cfh.reporting.fusion_schematic import _fusion_groups, render_fusion_schematic_svg

REPO_ROOT = Path(__file__).parent.parent
ETV6_RUN_DIR = (
    REPO_ROOT
    / "runs"
    / "cohort-scan_msk_impact_50k_2026_20260904T144201Z"
    / "cohort_scan"
    / "gene_reports"
    / "etv6"
)


def _payload() -> dict:
    return json.loads((ETV6_RUN_DIR / "results.json").read_text())


def test_real_etv6_svg_file_matches_pure_render_of_committed_payload():
    rendered = render_fusion_schematic_svg(_payload()) + "\n"
    on_disk = (ETV6_RUN_DIR / "visualizations" / "fusion_schematic.svg").read_text()
    assert on_disk == rendered


def test_real_etv6_rows_label_both_role_dependent_orientations():
    payload = _payload()
    role_counts = {
        role: sum(event.get("target_role") == role for event in payload["events"])
        for role in ("five_prime", "three_prime")
    }
    assert role_counts == {"five_prime": 67, "three_prime": 21}
    groups = _fusion_groups(payload)
    groups.sort(key=lambda group: (-group["count"], group["breakpoint_aa"], group["partner_gene"]))
    shown_groups = groups[:28]
    assert {group["role"] for group in shown_groups} == {"five_prime", "three_prime"}

    svg = render_fusion_schematic_svg(payload)
    root = ET.fromstring(svg)
    labels_by_y: dict[float, dict[str, ET.Element]] = {}
    for text in root.findall("{http://www.w3.org/2000/svg}text"):
        fusion_end = text.get("data-fusion-end")
        if fusion_end:
            labels_by_y.setdefault(float(text.attrib["y"]), {})[fusion_end] = text

    assert len(labels_by_y) == len(shown_groups)
    for group, (_y, labels) in zip(shown_groups, sorted(labels_by_y.items())):
        five_prime = labels["5-prime"]
        three_prime = labels["3-prime"]
        assert five_prime.text == "5'"
        assert three_prime.text == "3'"
        assert float(five_prime.attrib["x"]) < float(three_prime.attrib["x"])
        if group["role"] == "three_prime":
            assert five_prime.get("data-block") == "partner"
            assert three_prime.get("data-block") == "target"
        else:
            assert five_prime.get("data-block") == "target"
            assert three_prime.get("data-block") == "partner"
