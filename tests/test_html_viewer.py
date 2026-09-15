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
EML4_ALK_RUN_DIR = latest_run_dir("eml4-alk_msk-impact-50k-2026")


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

    def test_lollipop_points_are_wired_to_their_events_table_row(self, tmp_path):
        """A lollipop-track ``<circle data-event-id="...">`` (as emitted by
        ``cfh.real_benchmark._domain_track_svg``) must have a matching
        ``data-event-id`` on its row in the events table below, plus the
        click-to-scroll wiring in the page script, so clicking a point
        navigates to that sample's row."""
        run_dir = _write_gene_run(tmp_path, _minimal_gene_payload())
        viz_dir = run_dir / "visualizations"
        (viz_dir / "domain_retention_outliers.svg").write_text(
            '<svg xmlns="http://www.w3.org/2000/svg">'
            '<circle data-event-id="E2" r="3"><title>event E2</title></circle>'
            "</svg>"
        )
        html = build_run_viewer(run_dir).read_text()
        assert 'id="domain-track-svg"' in html
        assert 'data-event-id="E2"' in html.split('id="events-table"', 1)[1]
        assert "circle[data-event-id]" in html
        assert "row-highlight" in html

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


class TestCompositeScoreTable:
    def _payload_with_composite_score(self, **overrides) -> dict:
        payload = _minimal_gene_payload(**overrides)
        payload["algorithm_results"] = [
            {
                "Algorithm": "composite_score",
                "Algorithm_version": "0.1.0",
                "Parameters": {},
                "Summary": {"n_partners_ranked": 2},
                "Tables": {
                    "composite_evidence_ranking": [
                        {
                            "Partner_gene": "AGK",
                            "Event_count": 3,
                            "Composite_score": 0.91,
                            "Recurrence_score": 0.5,
                            "Domain_retention_score": 0.5,
                            "Domain_disruption_score": 0.1,
                            "Cutpoint_proximity_score": 0.9,
                            "Confidence_certainty_score": 0.7,
                            "Rank": 1,
                        },
                        {
                            "Partner_gene": "KIAA1549",
                            "Event_count": 1,
                            "Composite_score": 0.2,
                            "Recurrence_score": 0.1,
                            "Domain_retention_score": 0.1,
                            "Domain_disruption_score": 0.0,
                            "Cutpoint_proximity_score": 0.2,
                            "Confidence_certainty_score": 0.3,
                            "Rank": 2,
                        },
                    ]
                },
                "Warnings": [],
            }
        ]
        return payload

    def test_composite_score_table_is_rendered_with_sortable_headers(self, tmp_path):
        run_dir = _write_gene_run(tmp_path, self._payload_with_composite_score(), with_svgs=False)
        html = build_run_viewer(run_dir).read_text()
        assert "Composite evidence ranking" in html
        assert 'id="composite-score-table"' in html
        assert 'table class="data-table sortable"' in html
        assert 'data-sort-key="Composite_score"' in html
        assert 'data-sort="0.91"' in html
        assert "AGK" in html and "KIAA1549" in html

    def test_composite_score_section_omitted_when_algorithm_absent(self, tmp_path):
        run_dir = _write_gene_run(tmp_path, _minimal_gene_payload(), with_svgs=False)
        html = build_run_viewer(run_dir).read_text()
        assert "Composite evidence ranking" not in html
        assert "composite-score-table" not in html

    def test_composite_score_section_omitted_when_algorithm_failed(self, tmp_path):
        payload = self._payload_with_composite_score()
        payload["algorithm_results"][0]["Warnings"] = ["Algorithm failed: boom"]
        run_dir = _write_gene_run(tmp_path, payload, with_svgs=False)
        html = build_run_viewer(run_dir).read_text()
        assert "composite-score-table" not in html


class TestGenePairPage:
    def _gene_pair_payload(self) -> dict:
        return {
            "gene_symbol": "EML4-ALK",
            "study_id": "some_study",
            "summary": {
                "gene_pair": ["EML4", "ALK"],
                "eligible_event_count": 10,
                "observed_count": 8,
                "expected_count": 5.5,
                "fisher_p_value": 0.002,
                "is_enriched": True,
            },
            "events": [
                {
                    "event_id": "E1",
                    "sample_id": "S1",
                    "tumor_type": "Lung Cancer",
                    "oncotree_code": "LUAD",
                    "partner_gene": "EML4",
                    "target_role": "three_prime",
                    "breakpoint_protein_position": 1156,
                    "domain_status": "retained",
                }
            ],
            "algorithm_results": [
                {
                    "Algorithm": "joint_partner",
                    "Algorithm_version": "0.1.0",
                    "Parameters": {},
                    "Summary": {},
                    "Tables": {
                        "pair_results": [
                            {
                                "gene5": "EML4",
                                "gene3": "ALK",
                                "eligible_event_count": 10,
                                "observed_count": 8,
                                "expected_count": 5.5,
                                "p_value": 0.002,
                                "odds_ratio": None,
                            }
                        ]
                    },
                    "Warnings": [],
                }
            ],
        }

    def test_gene_pair_payload_uses_the_dedicated_template(self, tmp_path):
        run_dir = _write_gene_run(tmp_path, self._gene_pair_payload(), with_svgs=False)
        html = build_run_viewer(run_dir).read_text()
        assert "gene-pair fusion-hotspot report" in html
        assert "Gene-pair enrichment test" in html
        assert "EML4" in html and "ALK" in html
        assert "0.002" in html
        # Not force-fit into the single-gene "Domain retention" heading.
        assert "Domain retention / outlier lollipop track" not in html

    def test_gene_pair_page_still_has_tumor_type_filter_and_events(self, tmp_path):
        run_dir = _write_gene_run(tmp_path, self._gene_pair_payload(), with_svgs=False)
        html = build_run_viewer(run_dir).read_text()
        assert 'id="tumor-type-filter"' in html
        assert "Lung Cancer" in html

    def test_gene_pair_with_no_joint_partner_table_still_renders(self, tmp_path):
        payload = self._gene_pair_payload()
        payload["algorithm_results"] = []
        run_dir = _write_gene_run(tmp_path, payload, with_svgs=False)
        result = build_run_viewer(run_dir)
        assert result is not None
        assert "Gene-pair enrichment test" not in result.read_text()


class TestGenomicClusteringSection:
    def _payload_with_genomic_recurrence(self, *, determinable: bool) -> dict:
        payload = _minimal_gene_payload()
        summary = {
            "determinable": determinable,
            "reference_build": "GRCh37",
            "bin_size_bp": 1000,
            "genomic_vs_protein_clustering_note": "All 2 events share the exact same breakpoint.",
        }
        tables = (
            {
                "genomic_bin_recurrence": [
                    {"chromosome": "7", "bin_start": 140000, "bin_end": 141000, "n_events": 2}
                ],
                "exact_position_recurrence": [],
                "protein_position_genomic_spread": [
                    {
                        "protein_position_aa": 381,
                        "n_events": 2,
                        "n_distinct_genomic_positions": 1,
                        "genomic_span_bp": 0,
                        "spans_multiple_chromosomes": False,
                    }
                ],
            }
            if determinable
            else {
                "genomic_bin_recurrence": [],
                "exact_position_recurrence": [],
                "protein_position_genomic_spread": [],
            }
        )
        payload["algorithm_results"] = [
            {
                "Algorithm": "genomic_position_recurrence",
                "Algorithm_version": "0.1.0",
                "Parameters": {},
                "Summary": summary,
                "Tables": tables,
                "Warnings": [],
            }
        ]
        return payload

    def test_genomic_clustering_section_renders_when_determinable(self, tmp_path):
        payload = self._payload_with_genomic_recurrence(determinable=True)
        run_dir = _write_gene_run(tmp_path, payload, with_svgs=False)
        html = build_run_viewer(run_dir).read_text()
        assert "Genomic-coordinate clustering" in html
        assert "GRCh37" in html
        assert "381" in html
        assert "All 2 events share the exact same breakpoint." in html

    def test_genomic_clustering_section_omitted_when_not_determinable(self, tmp_path):
        payload = self._payload_with_genomic_recurrence(determinable=False)
        run_dir = _write_gene_run(tmp_path, payload, with_svgs=False)
        html = build_run_viewer(run_dir).read_text()
        assert "Genomic-coordinate clustering" not in html

    def test_genomic_clustering_section_omitted_when_algorithm_absent(self, tmp_path):
        run_dir = _write_gene_run(tmp_path, _minimal_gene_payload(), with_svgs=False)
        html = build_run_viewer(run_dir).read_text()
        assert "Genomic-coordinate clustering" not in html


class TestCrossCohortPanel:
    def _cross_cohort_report(self) -> dict:
        return {
            "gene_symbol": "BRAF",
            "cmh": {
                "statistic": 21.36,
                "p_value": 3.8e-06,
                "common_odds_ratio": 7.46,
                "informative_strata": 2,
                "total_strata": 3,
            },
            "strata": [
                {
                    "study_id": "msk_impact_50k_2026",
                    "association_direction": "positive",
                    "informative": True,
                },
                {
                    "study_id": "thca_tcga_pan_can_atlas_2018",
                    "association_direction": "neutral",
                    "informative": False,
                },
            ],
            "nominal_pooled_significant": True,
            "warnings": ["Shared sample IDs violate between-cohort independence."],
        }

    def test_cmh_panel_renders_from_sibling_file(self, tmp_path):
        run_dir = _write_gene_run(tmp_path, _minimal_gene_payload(), with_svgs=False)
        (run_dir / "cross_cohort_concordance.json").write_text(
            json.dumps(self._cross_cohort_report())
        )
        html = build_run_viewer(run_dir).read_text()
        assert "Cross-cohort concordance" in html
        assert "21.36" in html
        assert "msk_impact_50k_2026" in html
        assert "Shared sample IDs violate" in html

    def test_cmh_panel_omitted_when_sibling_file_absent(self, tmp_path):
        run_dir = _write_gene_run(tmp_path, _minimal_gene_payload(), with_svgs=False)
        html = build_run_viewer(run_dir).read_text()
        assert "Cross-cohort concordance" not in html

    def test_cmh_panel_omitted_when_sibling_file_is_malformed(self, tmp_path):
        run_dir = _write_gene_run(tmp_path, _minimal_gene_payload(), with_svgs=False)
        (run_dir / "cross_cohort_concordance.json").write_text("{not json")
        result = build_run_viewer(run_dir)
        assert result is not None
        assert "Cross-cohort concordance" not in result.read_text()


class TestConfigSourceBadge:
    def test_curated_gene_gets_curated_badge(self, tmp_path):
        run_dir = _write_gene_run(
            tmp_path, _minimal_gene_payload(gene_symbol="BRAF"), with_svgs=False
        )
        html = build_run_viewer(run_dir).read_text()
        assert "badge-curated" in html
        assert "Curated config" in html

    def test_unrecognized_gene_gets_auto_badge(self, tmp_path):
        run_dir = _write_gene_run(
            tmp_path, _minimal_gene_payload(gene_symbol="NOTAREALGENE"), with_svgs=False
        )
        html = build_run_viewer(run_dir).read_text()
        assert "badge-auto" in html
        assert "Auto-configured" in html

    def test_cohort_scan_gene_page_uses_recorded_config_source_not_a_recheck(self, tmp_path):
        run_dir = tmp_path / "cohort_run"
        cohort_dir = run_dir / "cohort_scan"
        cohort_dir.mkdir(parents=True)
        summary = {
            "study_id": "s",
            "genes": [{"gene_symbol": "NOTAREALGENE", "config_source": "curated"}],
            "significant_genes": [],
        }
        (cohort_dir / "summary.json").write_text(json.dumps(summary))
        (cohort_dir / "manhattan.svg").write_text('<svg xmlns="http://www.w3.org/2000/svg"></svg>')
        gene_reports_dir = cohort_dir / "gene_reports"
        (gene_reports_dir / "notarealgene").mkdir(parents=True)
        (gene_reports_dir / "notarealgene" / "results.json").write_text(
            json.dumps(_minimal_gene_payload(gene_symbol="NOTAREALGENE"))
        )
        build_run_viewer(run_dir)
        gene_html = (run_dir / "viewer" / "genes" / "notarealgene.html").read_text()
        # The scan already determined this (unusual/synthetic) case as curated;
        # the viewer must trust that recorded value rather than re-deriving it.
        assert "badge-curated" in gene_html


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


def test_real_braf_run_viewer_renders_the_real_composite_score_and_genomic_sections(tmp_path):
    """BRAF's real committed run has both a non-trivial ``composite_score``
    evidence ranking and a determinable ``genomic_position_recurrence``
    result -- both optional sections must actually show up with real data,
    not just degrade gracefully."""
    scratch = tmp_path / BRAF_RUN_DIR.name
    shutil.copytree(BRAF_RUN_DIR, scratch)
    payload = json.loads((BRAF_RUN_DIR / "results.json").read_text())
    algorithm_results = {item["Algorithm"]: item for item in payload["algorithm_results"]}

    html = build_run_viewer(scratch).read_text()
    assert 'id="composite-score-table"' in html
    top_partner = algorithm_results["composite_score"]["Tables"]["composite_evidence_ranking"][0][
        "Partner_gene"
    ]
    assert top_partner in html

    genomic_result = algorithm_results.get("genomic_position_recurrence")
    if genomic_result and genomic_result["Summary"].get("determinable"):
        assert "Genomic-coordinate clustering" in html
    else:
        assert "Genomic-coordinate clustering" not in html


def _format_number(value: float) -> str:
    return f"{value:.4g}"


def test_real_braf_run_viewer_renders_the_real_cross_cohort_cmh_panel_when_present(tmp_path):
    scratch = tmp_path / BRAF_RUN_DIR.name
    shutil.copytree(BRAF_RUN_DIR, scratch)
    sibling = BRAF_RUN_DIR / "cross_cohort_concordance.json"
    if not sibling.exists():
        pytest.skip("no committed cross_cohort_concordance.json next to the real BRAF run")
    report = json.loads(sibling.read_text())

    html = build_run_viewer(scratch).read_text()
    assert "Cross-cohort concordance" in html
    assert _format_number(report["cmh"]["statistic"]) in html


def test_real_braf_run_viewer_shows_the_curated_badge(tmp_path):
    scratch = tmp_path / BRAF_RUN_DIR.name
    shutil.copytree(BRAF_RUN_DIR, scratch)
    html = build_run_viewer(scratch).read_text()
    assert "badge-curated" in html


def test_real_eml4_alk_run_uses_the_gene_pair_template(tmp_path):
    scratch = tmp_path / EML4_ALK_RUN_DIR.name
    shutil.copytree(EML4_ALK_RUN_DIR, scratch)
    payload = json.loads((EML4_ALK_RUN_DIR / "results.json").read_text())

    html = build_run_viewer(scratch).read_text()
    assert "gene-pair fusion-hotspot report" in html
    assert "Gene-pair enrichment test" in html
    gene5, gene3 = payload["summary"]["gene_pair"]
    assert gene5 in html and gene3 in html
    assert "fetch(" not in html
