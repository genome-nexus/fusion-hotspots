"""Unit tests for the "honorable mentions" highly ranked non-FDR-significant tier
(:func:`cfh.cohort.outputs.build_honorable_mentions`) and the curated-gene
full-report selection (:func:`cfh.cohort.scan.genes_needing_full_report`),
plus a real-data check against the already-committed
``msk_impact_50k_2026`` cohort-scan run.
"""

from __future__ import annotations

import json
import re
from unittest.mock import MagicMock

from pypdf import PdfReader

from cfh.cohort.outputs import (
    DEFAULT_HONORABLE_MENTION_COUNT,
    build_honorable_mentions,
    write_cohort_scan_outputs,
)
from cfh.cohort.recurrence import GeneRecurrence, RecurrenceGateResult
from cfh.cohort.scan import CohortScanResult, GeneScanOutcome, genes_needing_full_report
from conftest import latest_run_dir

REAL_COHORT_SCAN_SUMMARY_JSON = (
    latest_run_dir("cohort-scan_msk_impact_50k_2026") / "cohort_scan" / "summary.json"
)


def _row(
    gene_symbol: str,
    *,
    fisher_p: float | None,
    q_value: float | None,
    significant: bool,
    n_events: int = 10,
    in_frame: float = 50.0,
    domain_retained: float = 50.0,
) -> dict:
    return {
        "gene_symbol": gene_symbol,
        "fisher_p_value": fisher_p,
        "min_fdr_adjusted_q_value": q_value,
        "fdr_significant": significant,
        "n_events_analyzed": n_events,
        "in_frame_percent": in_frame,
        "domain_retention_percent": domain_retained,
    }


class TestBuildHonorableMentions:
    def test_excludes_fdr_significant_genes(self):
        rows = [
            _row("SIG", fisher_p=1e-6, q_value=0.001, significant=True),
            _row("NOTSIG", fisher_p=0.01, q_value=0.2, significant=False),
        ]
        mentions = build_honorable_mentions(rows)
        assert [m["gene_symbol"] for m in mentions] == ["NOTSIG"]

    def test_ranks_by_raw_fisher_p_value_ascending(self):
        rows = [
            _row("WORST", fisher_p=0.04, q_value=0.3, significant=False),
            _row("BEST", fisher_p=0.001, q_value=0.15, significant=False),
            _row("MIDDLE", fisher_p=0.01, q_value=0.2, significant=False),
        ]
        mentions = build_honorable_mentions(rows)
        assert [m["gene_symbol"] for m in mentions] == ["BEST", "MIDDLE", "WORST"]
        assert [m["rank"] for m in mentions] == [1, 2, 3]

    def test_respects_limit(self):
        rows = [
            _row(f"G{i}", fisher_p=0.001 * (i + 1), q_value=0.2, significant=False)
            for i in range(20)
        ]
        mentions = build_honorable_mentions(rows, limit=5)
        assert len(mentions) == 5
        assert [m["gene_symbol"] for m in mentions] == [f"G{i}" for i in range(5)]

    def test_default_limit_matches_module_constant(self):
        rows = [
            _row(f"G{i}", fisher_p=0.001 * (i + 1), q_value=0.2, significant=False)
            for i in range(DEFAULT_HONORABLE_MENTION_COUNT + 5)
        ]
        assert len(build_honorable_mentions(rows)) == DEFAULT_HONORABLE_MENTION_COUNT

    def test_excludes_genes_with_no_p_value(self):
        rows = [
            _row("NOPVALUE", fisher_p=None, q_value=None, significant=False),
            _row("HASPVALUE", fisher_p=0.02, q_value=0.3, significant=False),
        ]
        mentions = build_honorable_mentions(rows)
        assert [m["gene_symbol"] for m in mentions] == ["HASPVALUE"]

    def test_note_is_precise_about_not_being_significant(self):
        rows = [_row("GENE", fisher_p=0.01, q_value=0.2, significant=False)]
        [mention] = build_honorable_mentions(rows)
        note = mention["note"].lower()
        assert "did not survive" in note
        assert "not a claim of statistical significance" in note

    def test_empty_when_every_gene_is_significant_or_has_no_p_value(self):
        rows = [
            _row("SIG", fisher_p=1e-6, q_value=0.001, significant=True),
            _row("NOPVALUE", fisher_p=None, q_value=None, significant=False),
        ]
        assert build_honorable_mentions(rows) == []

    def test_includes_a_gene_with_a_real_raw_p_value_but_no_computed_q_value(self):
        """A gene can have a real Fisher p-value (domain_retention ran) yet
        never contribute a hypothesis that produced a q-value (e.g. every
        p-value it did produce was excluded from the BH correction for some
        reason) -- ``min_fdr_adjusted_q_value`` is ``None`` and
        ``fdr_significant`` is the tri-state ``None`` (see
        ``build_summary_rows``), not ``False``. Such a gene is still a valid
        honorable-mentions candidate by raw p-value ranking."""
        rows = [
            _row("NOQVALUE", fisher_p=0.001, q_value=None, significant=None),
            _row("HASQVALUE", fisher_p=0.02, q_value=0.3, significant=False),
        ]
        mentions = build_honorable_mentions(rows)
        assert {m["gene_symbol"] for m in mentions} == {"NOQVALUE", "HASQVALUE"}
        # Ranked by raw p-value regardless of whether q was ever computed.
        assert [m["gene_symbol"] for m in mentions] == ["NOQVALUE", "HASQVALUE"]

    def test_note_never_claims_a_q_verdict_for_a_gene_with_no_computed_q_value(self):
        """Regression test: the single, unconditional
        ``_HONORABLE_MENTION_NOTE`` used to claim, for every candidate
        regardless of whether it actually had a q-value, that its
        "FDR-adjusted q-value [was] at or above the significance threshold".
        For a gene that never had a q-value computed at all, that is an
        invented verdict about a number that was never derived -- the note
        must instead say the q-value/FDR status is unknown, never assert a
        specific "at or above threshold" claim."""
        rows = [
            _row("NOQVALUE", fisher_p=0.001, q_value=None, significant=None),
            _row("HASQVALUE", fisher_p=0.02, q_value=0.3, significant=False),
        ]
        mentions = {m["gene_symbol"]: m for m in build_honorable_mentions(rows)}

        no_q_note = mentions["NOQVALUE"]["note"].lower()
        assert "at or above the significance threshold" not in no_q_note
        assert "no fdr-adjusted q-value was ever computed" in no_q_note
        assert "unknown" in no_q_note
        assert "not a claim of statistical significance" in no_q_note

        # The gene that DOES have a real q-value keeps the original,
        # precise wording -- this fix must not blur the two cases together.
        has_q_note = mentions["HASQVALUE"]["note"].lower()
        assert "at or above the significance threshold" in has_q_note
        assert "no fdr-adjusted q-value was ever computed" not in has_q_note


def _outcome(gene_symbol: str, config_source: str) -> GeneScanOutcome:
    return GeneScanOutcome(
        gene_symbol=gene_symbol,
        entrez_gene_id=1,
        distinct_patient_count=10,
        total_sv_count=10,
        config_source=config_source,
        status="ok",
        run=MagicMock(name=f"run-for-{gene_symbol}"),
    )


def _result(outcomes: list[GeneScanOutcome], significant_genes: list[str]) -> CohortScanResult:
    return CohortScanResult(
        study_id="test_study",
        min_distinct_patients=5,
        recurrence_gate=RecurrenceGateResult(
            study_id="test_study",
            min_distinct_patients=5,
            total_genes=len(outcomes),
            passing_genes=[
                GeneRecurrence(
                    hugo_gene_symbol=o.gene_symbol,
                    entrez_gene_id=1,
                    distinct_patient_count=10,
                    total_sv_count=10,
                )
                for o in outcomes
            ],
        ),
        curated_gene_count=sum(1 for o in outcomes if o.config_source == "curated"),
        auto_config_gene_count=sum(1 for o in outcomes if o.config_source == "auto"),
        unresolved_gene_count=0,
        gene_outcomes=outcomes,
        fdr_rows=[],
        significant_genes=significant_genes,
    )


class TestGenesNeedingFullReport:
    def test_covers_every_curated_gene_not_just_a_hardcoded_pair(self):
        """Regression test for the real bug found in this codebase: a
        hardcoded ``ALWAYS_FULL_REPORT_GENES = ("BRAF", "RET")`` tuple was
        never updated after ALK and NTRK1 curated configs were added, so
        those two hand-curated genes silently got no gene_reports/ output.
        Full-report selection must instead derive "curated" from each
        outcome's own ``config_source``."""
        outcomes = [
            _outcome("BRAF", "curated"),
            _outcome("RET", "curated"),
            _outcome("ALK", "curated"),
            _outcome("NTRK1", "curated"),
            _outcome("AUTOGENE", "auto"),
        ]
        result = _result(outcomes, significant_genes=[])
        assert set(genes_needing_full_report(result)) == {"BRAF", "RET", "ALK", "NTRK1"}

    def test_includes_fdr_significant_and_honorable_mention_genes(self):
        outcomes = [
            _outcome("CURATEDGENE", "curated"),
            _outcome("SIGGENE", "auto"),
            _outcome("HONORABLEGENE", "auto"),
            _outcome("PLAINGENE", "auto"),
        ]
        result = _result(outcomes, significant_genes=["SIGGENE"])
        full_report_genes = genes_needing_full_report(
            result, honorable_mention_genes={"HONORABLEGENE"}
        )
        assert set(full_report_genes) == {"CURATEDGENE", "SIGGENE", "HONORABLEGENE"}
        assert "PLAINGENE" not in full_report_genes

    def test_never_includes_a_gene_that_was_not_actually_scanned(self):
        outcomes = [_outcome("SCANNED", "curated")]
        result = _result(outcomes, significant_genes=[])
        full_report_genes = genes_needing_full_report(
            result, honorable_mention_genes={"NEVERSCANNED"}
        )
        assert full_report_genes == ["SCANNED"]


def _stub_out_full_gene_report_writing(monkeypatch):
    """These tests exercise the consolidated summary artifacts only; full
    per-gene report writing (:func:`cfh.real_benchmark.write_outputs`) needs
    a real ``RealBenchmarkRun``, not the ``MagicMock`` stand-in used here,
    and is already covered by ``test_cohort_scan_pipeline.py`` and the real
    live cohort-scan run."""
    monkeypatch.setattr("cfh.cohort.outputs.write_outputs", lambda *a, **k: {})


def test_summary_json_carries_an_additive_machine_readable_honorable_mentions_field(
    tmp_path, monkeypatch
):
    _stub_out_full_gene_report_writing(monkeypatch)
    outcomes = [
        _outcome("SIGGENE", "auto"),
        _outcome("HONORABLEGENE", "auto"),
        _outcome("PLAINGENE", "auto"),
    ]
    for outcome, (p_value, q_value) in zip(
        outcomes, [(1e-8, 0.001), (0.01, 0.3), (0.5, 0.9)], strict=True
    ):
        outcome.run.summary = {
            "total_fusions": 10,
            "in_frame_percent": 60.0,
            "kinase_retained_percent": 70.0,
            "fisher_p_value": p_value,
            "permutation_p_value": p_value,
        }
        outcome.run.results = []
    result = _result(outcomes, significant_genes=["SIGGENE"])
    result.fdr_rows = [
        {"gene": "SIGGENE", "bh_adjusted_q": 0.001},
        {"gene": "HONORABLEGENE", "bh_adjusted_q": 0.3},
        {"gene": "PLAINGENE", "bh_adjusted_q": 0.9},
    ]

    paths = write_cohort_scan_outputs(result, tmp_path / "runs", pdf=False)
    payload = json.loads(paths["summary_json"].read_text())

    assert "honorable_mentions" in payload  # new, additive field
    assert "genes" in payload  # existing field untouched
    mention_genes = {m["gene_symbol"] for m in payload["honorable_mentions"]}
    assert mention_genes == {"HONORABLEGENE", "PLAINGENE"}
    assert "SIGGENE" not in mention_genes
    for mention in payload["honorable_mentions"]:
        assert mention["gene_symbol"] != "SIGGENE"


def test_summary_markdown_labels_the_tier_precisely_without_claiming_significance(
    tmp_path, monkeypatch
):
    _stub_out_full_gene_report_writing(monkeypatch)
    outcomes = [_outcome("SIGGENE", "auto"), _outcome("HONORABLEGENE", "auto")]
    for outcome, (p_value, _q) in zip(outcomes, [(1e-8, 0.001), (0.01, 0.3)], strict=True):
        outcome.run.summary = {
            "total_fusions": 10,
            "in_frame_percent": 60.0,
            "kinase_retained_percent": 70.0,
            "fisher_p_value": p_value,
            "permutation_p_value": p_value,
        }
        outcome.run.results = []
    result = _result(outcomes, significant_genes=["SIGGENE"])
    result.fdr_rows = [
        {"gene": "SIGGENE", "bh_adjusted_q": 0.001},
        {"gene": "HONORABLEGENE", "bh_adjusted_q": 0.3},
    ]

    paths = write_cohort_scan_outputs(result, tmp_path / "runs", pdf=False)
    markdown = paths["summary_markdown"].read_text()

    assert "Honorable mentions" in markdown
    assert "HONORABLEGENE" in markdown
    heading_index = markdown.index("Honorable mentions")
    next_heading_index = markdown.index("\n## ", heading_index)
    section = markdown[heading_index:next_heading_index]
    assert "did not survive genome-wide" in section
    assert "not** a claim of statistical significance" in section
    # SIGGENE (the actually-significant gene) must not appear inside the
    # honorable-mentions table rows, only in the full scanned-genes table.
    mentions_table_start = section.index("| rank |")
    assert "SIGGENE" not in section[mentions_table_start:]


def test_summary_pdf_renders_the_honorable_mentions_section_as_real_text(tmp_path, monkeypatch):
    _stub_out_full_gene_report_writing(monkeypatch)
    outcomes = [_outcome("SIGGENE", "auto"), _outcome("HONORABLEGENE", "auto")]
    for outcome, (p_value, _q) in zip(outcomes, [(1e-8, 0.001), (0.01, 0.3)], strict=True):
        outcome.run.summary = {
            "total_fusions": 10,
            "in_frame_percent": 60.0,
            "kinase_retained_percent": 70.0,
            "fisher_p_value": p_value,
            "permutation_p_value": p_value,
        }
        outcome.run.results = []
    result = _result(outcomes, significant_genes=["SIGGENE"])
    result.fdr_rows = [
        {"gene": "SIGGENE", "bh_adjusted_q": 0.001},
        {"gene": "HONORABLEGENE", "bh_adjusted_q": 0.3},
    ]

    paths = write_cohort_scan_outputs(result, tmp_path / "runs", pdf=True)
    reader = PdfReader(str(paths["summary_pdf"]))
    text = "".join(page.extract_text() or "" for page in reader.pages)

    assert "Honorable mentions" in text
    assert "HONORABLEGENE" in text


def test_honorable_mention_count_is_configurable(tmp_path, monkeypatch):
    _stub_out_full_gene_report_writing(monkeypatch)
    outcomes = [_outcome(f"GENE{i}", "auto") for i in range(10)]
    for i, outcome in enumerate(outcomes):
        outcome.run.summary = {
            "total_fusions": 10,
            "in_frame_percent": 60.0,
            "kinase_retained_percent": 70.0,
            "fisher_p_value": 0.001 * (i + 1),
            "permutation_p_value": 0.001 * (i + 1),
        }
        outcome.run.results = []
    result = _result(outcomes, significant_genes=[])
    result.fdr_rows = [{"gene": o.gene_symbol, "bh_adjusted_q": 0.5} for o in outcomes]

    paths = write_cohort_scan_outputs(
        result, tmp_path / "runs", pdf=False, honorable_mention_count=3
    )
    payload = json.loads(paths["summary_json"].read_text())
    assert len(payload["honorable_mentions"]) == 3


def test_write_cohort_scan_outputs_wires_up_the_manuscript_report(tmp_path, monkeypatch):
    """Regression test for the ``cfh cohort-scan`` CLI path: every call to
    :func:`write_cohort_scan_outputs` (which the CLI always calls) must also
    produce the cross-gene manuscript report (``paper.md``/``paper.pdf``),
    not just as a one-off manual regeneration."""
    _stub_out_full_gene_report_writing(monkeypatch)
    outcomes = [_outcome("SIGGENE", "auto"), _outcome("HONORABLEGENE", "auto")]
    for outcome, (p_value, _q) in zip(outcomes, [(1e-8, 0.001), (0.01, 0.3)], strict=True):
        outcome.run.summary = {
            "total_fusions": 10,
            "in_frame_percent": 60.0,
            "kinase_retained_percent": 70.0,
            "fisher_p_value": p_value,
            "permutation_p_value": p_value,
        }
        outcome.run.results = []
    result = _result(outcomes, significant_genes=["SIGGENE"])
    result.fdr_rows = [
        {"gene": "SIGGENE", "bh_adjusted_q": 0.001},
        {"gene": "HONORABLEGENE", "bh_adjusted_q": 0.3},
    ]

    paths = write_cohort_scan_outputs(result, tmp_path / "runs", pdf=True)

    assert "manuscript_markdown" in paths
    assert "manuscript_pdf" in paths
    markdown = paths["manuscript_markdown"].read_text()
    assert markdown.startswith("# Genome-wide fusion-hotspot analysis of test_study")
    assert "## Abstract" in markdown
    assert "## Methods" in markdown
    assert "## Results" in markdown
    assert "## Discussion" in markdown
    assert "## Appendix: per-gene report index" in markdown
    assert "SIGGENE" in markdown
    assert "HONORABLEGENE" in markdown

    reader = PdfReader(str(paths["manuscript_pdf"]))
    text = "".join(page.extract_text() or "" for page in reader.pages)
    assert "Abstract" in text
    assert "SIGGENE" in text
    assert "HONORABLEGENE" in text


def test_write_cohort_scan_outputs_records_algorithms_actually_run_not_the_live_registry(
    tmp_path, monkeypatch
):
    """Regression test: the Methods section must describe the algorithm
    suite THIS run's own results actually recorded, not a live import of
    ``cfh.algorithms.registry.list_algorithms()`` -- using algorithm names
    that are not real registered algorithms proves the text is sourced from
    ``outcome.run.results``, not the registry."""
    _stub_out_full_gene_report_writing(monkeypatch)
    outcomes = [_outcome("SIGGENE", "auto")]
    fake_result_1 = MagicMock()
    fake_result_1.Algorithm = "totally_fake_algorithm"
    fake_result_2 = MagicMock()
    fake_result_2.Algorithm = "another_fake_algorithm"
    outcomes[0].run.summary = {
        "total_fusions": 10,
        "in_frame_percent": 60.0,
        "kinase_retained_percent": 70.0,
        "fisher_p_value": 1e-8,
        "permutation_p_value": 1e-8,
    }
    outcomes[0].run.results = [fake_result_1, fake_result_2]
    result = _result(outcomes, significant_genes=["SIGGENE"])
    result.fdr_rows = [{"gene": "SIGGENE", "bh_adjusted_q": 0.001}]

    paths = write_cohort_scan_outputs(result, tmp_path / "runs", pdf=False)

    payload = json.loads(paths["summary_json"].read_text())
    assert payload["algorithms_run"] == ["another_fake_algorithm", "totally_fake_algorithm"]

    markdown = paths["manuscript_markdown"].read_text()
    match = re.search(r"algorithm suite recorded for this scan \(([^)]*)\)", markdown)
    assert match is not None
    assert {name.strip() for name in match.group(1).split(",")} == {
        "another_fake_algorithm",
        "totally_fake_algorithm",
    }


def test_real_committed_run_honorable_mentions_are_precise_and_ranked_by_p_value():
    """Validates against the real, already-committed genome-wide
    ``msk_impact_50k_2026`` run: only ETV6 is FDR-significant, but RET,
    FGFR2, ALK, EGFR, BRAF, FGFR3, and NTRK1 form a real biologically
    sensible second tier of well-known fusion-driver genes that should be
    visible as honorable mentions, not silently folded into 543
    undifferentiated non-significant rows.
    """
    payload = json.loads(REAL_COHORT_SCAN_SUMMARY_JSON.read_text())
    rows = payload["genes"]
    mentions = build_honorable_mentions(rows)

    assert mentions, "expected a non-empty honorable-mentions tier"
    mention_genes = [m["gene_symbol"] for m in mentions]
    assert "ETV6" not in mention_genes  # the one real FDR-significant gene

    expected_second_tier = {"RET", "FGFR2", "ALK", "EGFR", "BRAF", "FGFR3", "NTRK1"}
    assert expected_second_tier <= set(mention_genes[: len(expected_second_tier) + 3])

    # Strictly ascending raw p-value order.
    p_values = [m["fisher_p_value"] for m in mentions]
    assert p_values == sorted(p_values)

    for mention in mentions:
        assert mention["gene_symbol"] not in payload["significant_genes"]
