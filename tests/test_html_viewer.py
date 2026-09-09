"""Tests for the Milestone-4 backend-less static-site HTML viewer
(:mod:`cfh.reporting.html_viewer`): both synthetic-payload unit tests and
real-run tests against the already-committed BRAF standalone run and the
genome-wide cohort-scan run.

Every real-run test reads exactly the ``results.json``/``summary.json``/
``*.svg`` artifacts already committed under ``runs/`` -- no network access,
nothing here recomputes a statistic, and the viewer builder must never
modify any pre-existing file (see
``test_building_the_viewer_never_modifies_any_pre_existing_run_file``).
"""

from __future__ import annotations

import hashlib
import json
import re
import shutil
from pathlib import Path

import pytest

from cfh.reporting.html_viewer import build_run_viewer
from conftest import latest_run_dir

BRAF_RUN_DIR = latest_run_dir("braf_msk-impact-50k-2026")
COHORT_SCAN_RUN_DIR = latest_run_dir("cohort-scan_msk_impact_50k_2026")


def _minimal_gene_payload(**overrides) -> dict:
    payload = {
        "gene_symbol": "BRAF",
        "study_id": "some_study",
        "summary": {
            "total_fusions": 10,
            "mapped_fusions": 9,
            "fisher_p_value": 0.0123,
            "partner_counts": [{"Partner_gene": "AGK", "Event_count": 3}],
        },
        "events": [
            {
                "event_id": "E1",
                "sample_id": "S1",
                "tumor_type": "Melanoma",
                "oncotree_code": "SKCM",
                "partner_gene": "AGK",
                "target_role": "three_prime",
                "breakpoint_protein_position": 327,
                "domain_status": "retained",
            },
            {
                "event_id": "E2",
                "sample_id": "S2",
                "tumor_type": "Thyroid Cancer",
                "oncotree_code": "THPA",
                "partner_gene": "KIAA1549",
                "target_role": "three_prime",
                "breakpoint_protein_position": 381,
                "domain_status": "lost",
            },
        ],
    }
    payload.update(overrides)
    return payload


def _write_gene_run(tmp_path: Path, payload: dict, *, with_svgs: bool = True) -> Path:
    run_dir = tmp_path / "some_run"
    run_dir.mkdir()
    (run_dir / "results.json").write_text(json.dumps(payload))
    if with_svgs:
        viz_dir = run_dir / "visualizations"
        viz_dir.mkdir()
        (viz_dir / "domain_retention_outliers.svg").write_text(
            '<svg xmlns="http://www.w3.org/2000/svg"><circle r="1"/></svg>'
        )
        (viz_dir / "fusion_schematic.svg").write_text(
            '<svg xmlns="http://www.w3.org/2000/svg">'
            "<rect><title>Partner gene AGK; BRAF–AGK fusion, breakpoint BRAF aa 327; "
            "samples S1</title></rect>"
            "<rect><title>Partner gene KIAA1549; BRAF–KIAA1549 fusion, breakpoint BRAF "
            "aa 381; samples S2</title></rect>"
            "</svg>"
        )
    return run_dir


class TestBuildRunViewerDispatch:
    def test_returns_none_for_a_directory_with_neither_artifact(self, tmp_path):
        empty_dir = tmp_path / "nothing_here"
        empty_dir.mkdir()
        assert build_run_viewer(empty_dir) is None

    def test_returns_none_for_a_nonexistent_directory(self, tmp_path):
        assert build_run_viewer(tmp_path / "does_not_exist") is None

    def test_single_gene_run_writes_viewer_index_html(self, tmp_path):
        run_dir = _write_gene_run(tmp_path, _minimal_gene_payload())
        result = build_run_viewer(run_dir)
        assert result == run_dir / "viewer" / "index.html"
        assert result.exists()

    def test_cohort_scan_run_is_detected_over_single_gene_shape(self, tmp_path):
        run_dir = tmp_path / "cohort_run"
        cohort_dir = run_dir / "cohort_scan"
        cohort_dir.mkdir(parents=True)
        (cohort_dir / "summary.json").write_text(
            json.dumps({"study_id": "s", "genes": [], "significant_genes": []})
        )
        (cohort_dir / "manhattan.svg").write_text('<svg xmlns="http://www.w3.org/2000/svg"></svg>')
        result = build_run_viewer(run_dir)
        assert result == run_dir / "viewer" / "index.html"
        assert result.exists()


class TestSingleGenePageContent:
    def test_page_embeds_gene_symbol_and_study_id(self, tmp_path):
        run_dir = _write_gene_run(tmp_path, _minimal_gene_payload())
        html = build_run_viewer(run_dir).read_text()
        assert "BRAF" in html
        assert "some_study" in html

    def test_svgs_are_inlined_not_linked_or_fetched(self, tmp_path):
        run_dir = _write_gene_run(tmp_path, _minimal_gene_payload())
        html = build_run_viewer(run_dir).read_text()
        # The actual <svg ...> markup is part of the page's own DOM ...
        assert html.count("<svg") >= 2
        # ... and nothing in the page ever fetches a sibling file over
        # file://, which is blocked by CORS in some browsers.
        assert "fetch(" not in html
        assert 'src="' not in html
        assert "XMLHttpRequest" not in html

    def test_events_are_embedded_as_inline_json_not_a_separate_file(self, tmp_path):
        run_dir = _write_gene_run(tmp_path, _minimal_gene_payload())
        html = build_run_viewer(run_dir).read_text()
        assert 'id="events-data"' in html
        assert "Melanoma" in html
        assert "Thyroid Cancer" in html

    def test_events_json_round_trips_and_is_not_truncated(self, tmp_path):
        run_dir = _write_gene_run(tmp_path, _minimal_gene_payload())
        html = build_run_viewer(run_dir).read_text()
        match = re.search(
            r'<script type="application/json" id="events-data">(.*?)</script>', html, re.DOTALL
        )
        assert match is not None
        events = json.loads(match.group(1))
        assert len(events) == 2
        assert {event["tumor_type"] for event in events} == {"Melanoma", "Thyroid Cancer"}

    def test_summary_scalar_stats_are_shown(self, tmp_path):
        run_dir = _write_gene_run(tmp_path, _minimal_gene_payload())
        html = build_run_viewer(run_dir).read_text()
        assert "0.0123" in html
        # A structured (list-valued) summary field must not leak into the
        # simple stat grid -- it already has a dedicated home in report.md.
        assert "Partner_gene" not in html

    def test_missing_svgs_degrade_gracefully_to_an_events_only_page(self, tmp_path):
        run_dir = _write_gene_run(tmp_path, _minimal_gene_payload(), with_svgs=False)
        html = build_run_viewer(run_dir).read_text()
        assert "<svg" not in html
        assert 'id="events-data"' in html
        assert "<!doctype html>" in html.lower()

    def test_gene_pair_shaped_payload_with_no_gene_track_does_not_crash(self, tmp_path):
        payload = _minimal_gene_payload(
            gene_symbol="EML4-ALK", summary={"gene_pair": ["EML4", "ALK"], "observed_count": 6}
        )
        run_dir = _write_gene_run(tmp_path, payload, with_svgs=False)
        result = build_run_viewer(run_dir)
        assert result is not None
        assert "EML4-ALK" in result.read_text()

    def test_html_is_escaped_against_malicious_or_odd_field_values(self, tmp_path):
        payload = _minimal_gene_payload(gene_symbol="<script>alert(1)</script>")
        run_dir = _write_gene_run(tmp_path, payload, with_svgs=False)
        html = build_run_viewer(run_dir).read_text()
        assert "<script>alert(1)</script>" not in html
        assert "&lt;script&gt;" in html

    def test_json_script_embedding_cannot_be_broken_out_of_by_a_field_value(self, tmp_path):
        payload = _minimal_gene_payload()
        payload["events"][0]["sample_id"] = "S1</script><script>evil()</script>"
        run_dir = _write_gene_run(tmp_path, payload, with_svgs=False)
        html = build_run_viewer(run_dir).read_text()
        assert "</script><script>evil()" not in html


class TestCohortScanLandingPage:
    def _write_cohort_run(self, tmp_path: Path) -> Path:
        run_dir = tmp_path / "cohort_run"
        cohort_dir = run_dir / "cohort_scan"
        cohort_dir.mkdir(parents=True)
        summary = {
            "study_id": "msk_impact_50k_2026",
            "total_genes_before_gating": 500,
            "genes_after_gating": 20,
            "significant_genes": ["BRAF"],
            "significance_level": 0.05,
            "curated_gene_count": 5,
            "auto_config_gene_count": 15,
            "genes": [],
        }
        (cohort_dir / "summary.json").write_text(json.dumps(summary))
        (cohort_dir / "manhattan.svg").write_text(
            '<svg xmlns="http://www.w3.org/2000/svg">'
            '<circle data-gene="BRAF" data-significant="true" cx="1" cy="1" r="1"/>'
            '<circle data-gene="NOTANALYZED" data-significant="false" cx="2" cy="2" r="1"/>'
            "</svg>"
        )
        gene_reports_dir = cohort_dir / "gene_reports"
        (gene_reports_dir / "braf").mkdir(parents=True)
        (gene_reports_dir / "braf" / "results.json").write_text(json.dumps(_minimal_gene_payload()))
        return run_dir

    def test_landing_page_embeds_manhattan_svg_and_links_available_genes(self, tmp_path):
        run_dir = self._write_cohort_run(tmp_path)
        index_path = build_run_viewer(run_dir)
        html = index_path.read_text()
        assert "<svg" in html
        assert 'data-gene="BRAF"' in html
        assert 'href="genes/braf.html"' in html
        assert (run_dir / "viewer" / "genes" / "braf.html").exists()

    def test_gene_without_a_full_report_is_not_linked(self, tmp_path):
        run_dir = self._write_cohort_run(tmp_path)
        html = build_run_viewer(run_dir).read_text()
        assert "notanalyzed.html" not in html.lower()

    def test_gene_page_has_a_back_link_to_the_landing_page(self, tmp_path):
        run_dir = self._write_cohort_run(tmp_path)
        build_run_viewer(run_dir)
        gene_html = (run_dir / "viewer" / "genes" / "braf.html").read_text()
        assert 'href="../index.html"' in gene_html


@pytest.mark.parametrize(
    "run_dir", [BRAF_RUN_DIR, COHORT_SCAN_RUN_DIR], ids=["braf-standalone", "cohort-scan"]
)
def test_viewer_builds_successfully_against_the_real_committed_run(tmp_path, run_dir):
    """No network, no synthetic fixtures: builds the viewer straight from
    the exact run directories already committed to the repo, into a scratch
    copy so the real ``runs/`` tree is untouched by this test."""
    scratch = tmp_path / run_dir.name
    shutil.copytree(run_dir, scratch)
    result = build_run_viewer(scratch)
    assert result is not None
    assert result.exists()
    html = result.read_text()
    assert "<!doctype html>" in html.lower()
    assert "fetch(" not in html


def test_building_the_viewer_never_modifies_any_pre_existing_run_file(tmp_path):
    """Purely additive: building the viewer for the real committed BRAF run
    must not change a single byte of results.json/results.tsv/report.md/
    the existing SVGs -- only new files under viewer/ may appear."""
    scratch = tmp_path / BRAF_RUN_DIR.name
    shutil.copytree(BRAF_RUN_DIR, scratch)

    def _hashes() -> dict[str, str]:
        return {
            str(path.relative_to(scratch)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted(scratch.rglob("*"))
            if path.is_file()
        }

    before = _hashes()
    build_run_viewer(scratch)
    after = _hashes()

    for relative_path, digest in before.items():
        assert after.get(relative_path) == digest, f"{relative_path} was modified"


def test_real_braf_run_viewer_populates_the_real_tumor_types_and_domain_statuses(tmp_path):
    scratch = tmp_path / BRAF_RUN_DIR.name
    shutil.copytree(BRAF_RUN_DIR, scratch)
    payload = json.loads((BRAF_RUN_DIR / "results.json").read_text())
    real_tumor_types = {
        event["tumor_type"] for event in payload["events"] if event.get("tumor_type")
    }

    html = build_run_viewer(scratch).read_text()
    match = re.search(
        r'<script type="application/json" id="events-data">(.*?)</script>', html, re.DOTALL
    )
    embedded_tumor_types = {event["tumor_type"] for event in json.loads(match.group(1))}
    assert embedded_tumor_types == real_tumor_types
    assert "Melanoma" in html  # a real BRAF fusion tumor type, sanity-checked directly


def test_real_cohort_scan_viewer_links_every_gene_report_directory(tmp_path):
    scratch = tmp_path / COHORT_SCAN_RUN_DIR.name
    shutil.copytree(COHORT_SCAN_RUN_DIR, scratch)
    real_gene_reports_dir = COHORT_SCAN_RUN_DIR / "cohort_scan" / "gene_reports"
    expected_genes = sorted(p.name for p in real_gene_reports_dir.iterdir() if p.is_dir())
    assert expected_genes  # the real committed cohort scan does have gene reports

    build_run_viewer(scratch)
    viewer_dir = scratch / "viewer"
    for gene_lower in expected_genes:
        assert (viewer_dir / "genes" / f"{gene_lower}.html").exists()
    index_html = (viewer_dir / "index.html").read_text()
    assert f'href="genes/{expected_genes[0]}.html"' in index_html
