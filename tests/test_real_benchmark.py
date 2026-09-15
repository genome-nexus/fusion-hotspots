import csv
import json
import re
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

import pandas as pd
import pytest
import requests
from click.testing import CliRunner

from cfh import cli
from cfh import real_benchmark as benchmark_module
from cfh.algorithms import domain_disruption as disruption_module
from cfh.algorithms import domain_retention as retention_module
from cfh.algorithms.mutation_cooccurrence import comparator_target_key
from cfh.genes.registry import load_gene_config
from cfh.ingestion import cbioportal_api
from cfh.mapping.genome_nexus_source import GenomeNexusClient
from cfh.real_benchmark import (
    RealBenchmarkNetworkError,
    _target_breakpoint,
    analyze_structural_variant_calls,
    analyze_structural_variant_calls_with_config,
    run_real_benchmark,
    write_outputs,
)
from cfh.reporting.text import render_abstract


def test_target_breakpoint_uses_genome_nexus_locus_when_cbioportal_site_labels_are_swapped():
    """Regression for P-0053901-T01-IM6's malformed STRN-ALK source row.

    The live source labels site2 as ALK but pairs it with STRN's position
    (37145859). The actual ALK coordinate is in site1 (29446375), inside
    Genome Nexus's canonical ALK exon-spanned locus (29415640-30144432).
    """
    row = {
        "Site1_Hugo_Symbol": "STRN",
        "Site1_Position": 29446375,
        "Site2_Hugo_Symbol": "ALK",
        "Site2_Position": 37145859,
    }

    assert _target_breakpoint(row, "ALK", (29415640, 30144432)) == 29446375


def test_no_exon_coordinates_skip_target_locus_records_without_crashing(
    genome_nexus_canonical_transcript_fixture_path,
):
    """Genome Nexus returned this shape for NCOA4 in a cohort-wide scan."""
    client = _genome_nexus_client(genome_nexus_canonical_transcript_fixture_path)
    transcript = json.loads(genome_nexus_canonical_transcript_fixture_path.read_text())
    transcript["exons"] = []
    client.fetch_canonical_transcript.return_value = transcript

    run = analyze_structural_variant_calls(
        [_call("NCOA4-NO-EXONS")],
        "BRAF",
        "msk_impact_50k_2026",
        genome_nexus_client=client,
        n_permutations=5,
    )

    assert run.summary["total_fusions"] == 1
    assert run.summary["mapped_fusions"] == 0
    assert run.summary["skipped_fusions"] == 1
    assert any(
        "Genome Nexus returned no exon coordinates for target-locus validation" in warning
        for warning in run.warnings
    )


@pytest.mark.parametrize(
    ("row", "message"),
    [
        (
            {"Site1_Position": 10_000_000, "Site2_Position": 20_000_000},
            "no site breakpoint within target locus",
        ),
        (
            {"Site1_Position": 150, "Site2_Position": 175},
            "ambiguous genomic breakpoints within target locus",
        ),
    ],
)
def test_target_breakpoint_rejects_unmappable_or_ambiguous_source_positions(row, message):
    with pytest.raises(ValueError, match=message):
        _target_breakpoint(row, "TARGET", (100, 200))


def _genome_nexus_client(fixture_path):
    client = MagicMock(spec=GenomeNexusClient)
    client.fetch_canonical_transcript.return_value = json.loads(fixture_path.read_text())
    return client


def _call(
    sample_id,
    *,
    event_info="Protein Fusion: in frame  {KIAA1549:BRAF}",
    breakpoint=140493152,
    connection_type="3to3",
):
    return {
        "sampleId": sample_id,
        "site1HugoSymbol": "KIAA1549",
        "site2HugoSymbol": "BRAF",
        "site2Position": breakpoint,
        "site2EffectOnFrame": "NA",
        "connectionType": connection_type,
        "eventInfo": event_info,
    }


@pytest.mark.parametrize("override_budget", [None, 11])
def test_disruption_uses_genomic_null_and_requested_budget(
    genome_nexus_canonical_transcript_fixture_path, monkeypatch, override_budget
):
    client = _genome_nexus_client(genome_nexus_canonical_transcript_fixture_path)
    retention_spy = MagicMock(wraps=retention_module.permutation_null_test)
    disruption_spy = MagicMock(wraps=disruption_module.permutation_null_test)
    monkeypatch.setattr(retention_module, "permutation_null_test", retention_spy)
    monkeypatch.setattr(disruption_module, "permutation_null_test", disruption_spy)
    overrides = (
        {"domain_disruption": {"n_permutations": override_budget}}
        if override_budget is not None
        else {}
    )
    run = analyze_structural_variant_calls(
        [_call("S1")],
        "BRAF",
        "test",
        genome_nexus_client=client,
        n_permutations=7,
        algorithm_names=["domain_retention", "domain_disruption"],
        algorithm_params=overrides,
    )
    assert all(not result.Warnings for result in run.results)
    assert retention_spy.call_args.kwargs["genome_nexus_client"] is client
    assert retention_spy.call_args.kwargs["n_permutations"] == 7
    assert disruption_spy.call_args.kwargs["genome_nexus_client"] is client
    assert disruption_spy.call_args.kwargs["n_permutations"] == (override_budget or 7)


def test_real_benchmark_pipeline_writes_tsv_json_and_markdown(
    tmp_path,
    genome_nexus_canonical_transcript_fixture_path,
):
    client = _genome_nexus_client(genome_nexus_canonical_transcript_fixture_path)
    calls = [
        {
            "sampleId": "SAMPLE-1",
            "site1HugoSymbol": "KIAA1549",
            "site2HugoSymbol": "BRAF",
            "site2Position": 140493152,
            "site2EffectOnFrame": "NA",
            "connectionType": "3to3",
            "eventInfo": "Protein Fusion: in frame  {KIAA1549:BRAF}",
        },
        {
            "sampleId": "SAMPLE-2",
            "site1HugoSymbol": "BRAF",
            "site2HugoSymbol": "AGK",
            "site1Position": 140493152,
            "site2EffectOnFrame": "out-of-frame",
            "connectionType": "5to5",
            "eventInfo": "Protein Fusion: in frame  {BRAF:AGK}",
        },
    ]

    run = analyze_structural_variant_calls(
        calls,
        "BRAF",
        "msk_impact_50k_2026",
        genome_nexus_client=client,
        n_permutations=5,
    )
    paths = write_outputs(run, tmp_path, output_stem="benchmark")

    assert run.summary["total_fusions"] == 2
    assert run.summary["mapped_fusions"] == 2
    assert run.summary["in_frame_count"] == 1
    # Both records mapped successfully (no skips); any warning present is the
    # fusion-annotation QA cross-check flagging its own known
    # directional-intronic-breakpoint-snapping difference from this
    # project's nearest-exon-by-distance approximation (see
    # cfh.mapping.fusion_annotation_crosscheck), never a processing failure.
    assert all("Fusion-annotation cross-check" in warning for warning in run.warnings)
    assert paths["tsv"].read_text().splitlines()[0].startswith("event_id\tsample_id")
    payload = json.loads(paths["json"].read_text())
    assert payload["summary"]["domain_accession"] == "PF07714"
    assert len(payload["events"]) == 2
    assert "fusion_annotation_crosscheck" in payload["summary"]
    for event_row in payload["events"]:
        assert "fusion_annotation_crosscheck_agrees" in event_row
    assert paths["tsv"].name == "results.tsv"
    assert paths["json"].name == "results.json"
    assert paths["markdown"].name == "report.md"
    assert paths["manifest"].name == "manifest.json"
    assert paths["tsv"].parent.name.startswith("braf_msk-impact-50k-2026_")
    manifest = json.loads(paths["manifest"].read_text())
    assert manifest["gene"] == "BRAF"
    assert manifest["study_id"] == "msk_impact_50k_2026"
    assert len(manifest["endpoints_used"]) == 2
    assert paths["domain_svg"].parent.name == "visualizations"
    domain_svg = paths["domain_svg"].read_text()
    assert "#d62728" in domain_svg
    assert "#2878b5" in domain_svg
    assert "#f2a93b" in domain_svg
    assert "#777777" in domain_svg
    assert "fully retained" in domain_svg
    assert "truncated" in domain_svg
    assert "fully lost" in domain_svg
    assert "reference 100.0%" in paths["comparison_svg"].read_text()
    outliers = paths["outliers"].read_text()
    assert "reference_discrepancy" in outliers
    assert "source_vs_derived_qa_mismatch" in outliers
    assert "Protein Fusion: in frame" in outliers
    outlier_rows = list(csv.DictReader(paths["outliers"].open(), delimiter="\t"))
    qa_event_ids = {
        row["event_id"]
        for row in outlier_rows
        if row["discrepancy_type"] == "source_vs_derived_qa_mismatch"
    }
    assert qa_event_ids == {"EVT-SAMPLE-2-2"}
    # Every discrepancy/outlier row also states the mapped breakpoint's
    # protein position AND the exon it falls in (or was clamped to, for an
    # intronic breakpoint), not amino-acid position alone.
    assert outlier_rows
    for row in outlier_rows:
        assert row["breakpoint_protein_position"] == "380"
        assert row["breakpoint_exon"] == "8"
        assert row["is_intronic_breakpoint"] == "True"
    report = paths["markdown"].read_text()
    assert "does **not** reproduce" in report
    assert "Protein kinase domain (458-712 aa)" in report
    svg_paths = sorted(paths["run_directory"].glob("visualizations/*.svg"))
    assert svg_paths
    for svg_path in svg_paths:
        relative_path = svg_path.relative_to(paths["run_directory"]).as_posix()
        assert re.search(rf"!\[[^]]+\]\({re.escape(relative_path)}\)", report)
    assert "![Domain retention diagram](visualizations/domain_retention_outliers.svg)" in report
    assert "![Reference comparison](visualizations/reference_comparison.svg)" in report


def test_outputs_without_benchmark_reference_omit_reference_svg_and_section(
    tmp_path,
    genome_nexus_canonical_transcript_fixture_path,
):
    client = _genome_nexus_client(genome_nexus_canonical_transcript_fixture_path)
    config_without_reference = load_gene_config("BRAF").model_copy(
        update={"benchmark_reference": None}
    )
    run = analyze_structural_variant_calls_with_config(
        [_call("NO-REFERENCE")],
        config_without_reference,
        "msk_impact_50k_2026",
        genome_nexus_client=client,
        n_permutations=5,
    )

    paths = write_outputs(run, tmp_path, pdf=False)

    report = paths["markdown"].read_text()
    assert "## Reference comparison" not in report
    assert "reference_comparison.svg" not in report
    assert "comparison_svg" not in paths
    assert not (paths["run_directory"] / "visualizations" / "reference_comparison.svg").exists()


def test_markdown_cutpoint_line_states_the_exon_the_inferred_position_falls_in(
    tmp_path,
    genome_nexus_canonical_transcript_fixture_path,
):
    """The ``Cutpoint detection`` line in report.md is a specific breakpoint
    amino-acid position, so it must also state the protein exon it falls
    in -- not the aa number alone.
    """
    from cfh.algorithms.registry import list_algorithms

    client = _genome_nexus_client(genome_nexus_canonical_transcript_fixture_path)

    def retained_call(sample):
        return {
            "sampleId": sample,
            "site1HugoSymbol": "KIAA1549",
            "site2HugoSymbol": "BRAF",
            "site2Position": 140493152,
            "site2EffectOnFrame": "NA",
            "connectionType": "3to3",
            "eventInfo": "Protein Fusion: in frame  {KIAA1549:BRAF}",
        }

    def lost_call(sample):
        return {
            "sampleId": sample,
            "site1HugoSymbol": "BRAF",
            "site2HugoSymbol": "SOMEPARTNER",
            "site1Position": 140534500,
            "site2EffectOnFrame": "NA",
            "connectionType": "5to5",
            "eventInfo": "Protein Fusion: in frame  {BRAF:SOMEPARTNER}",
        }

    calls = [retained_call("S1"), retained_call("S2"), lost_call("S3"), lost_call("S4")]
    run = analyze_structural_variant_calls(
        calls,
        "BRAF",
        "msk_impact_50k_2026",
        genome_nexus_client=client,
        n_permutations=50,
        algorithm_names=list_algorithms(),
    )
    cutpoint_result = next(r for r in run.results if r.Algorithm == "cutpoint_detection")
    assert cutpoint_result.Summary["determinable"] is True
    # aa 138 sits strictly inside exon 3 (81-168 aa) of this fixture's
    # canonical transcript -- not on a shared exon/intron boundary, so a
    # single, unambiguous exon number is expected here.
    assert cutpoint_result.Summary["inferred_cutpoint_aa"] == 138

    paths = write_outputs(run, tmp_path, output_stem="benchmark", pdf=False)
    report = paths["markdown"].read_text()
    cutpoint_lines = [line for line in report.splitlines() if line.startswith("- Cutpoint")]

    assert cutpoint_lines == [
        "- Cutpoint detection: inferred breakpoint 138 aa (exon 3); "
        "corrected permutation p=0.431373."
    ]


def test_write_outputs_renders_pdf_report_by_default_and_can_be_disabled(
    tmp_path,
    genome_nexus_canonical_transcript_fixture_path,
):
    from pypdf import PdfReader

    client = _genome_nexus_client(genome_nexus_canonical_transcript_fixture_path)
    calls = [
        _call("SAMPLE-1"),
        _call("SAMPLE-2", event_info="Protein Fusion: in frame  {BRAF:AGK}"),
    ]
    run = analyze_structural_variant_calls(
        calls,
        "BRAF",
        "msk_impact_50k_2026",
        genome_nexus_client=client,
        n_permutations=5,
    )

    paths = write_outputs(run, tmp_path, run_id="with-pdf")

    assert paths["pdf"].name == "report.pdf"
    assert paths["pdf"].exists()
    reader = PdfReader(str(paths["pdf"]))
    assert len(reader.pages) >= 1
    text = "".join(page.extract_text() or "" for page in reader.pages)
    assert "BRAF" in text
    assert "Abstract" in text
    assert "Results summary" in text

    paths_no_pdf = write_outputs(run, tmp_path, run_id="without-pdf", pdf=False)

    assert "pdf" not in paths_no_pdf
    assert not (paths_no_pdf["run_directory"] / "report.pdf").exists()


def test_malformed_fusion_rows_are_warned_and_skipped_without_losing_valid_rows(
    genome_nexus_canonical_transcript_fixture_path,
):
    client = _genome_nexus_client(genome_nexus_canonical_transcript_fixture_path)
    calls = [
        _call("VALID"),
        _call(
            "BAD-ROLE",
            event_info="Protein Fusion: in frame; order unavailable",
            connection_type="3to3",
        ),
        _call("NO-BREAKPOINT", breakpoint=None),
    ]

    run = analyze_structural_variant_calls(
        calls,
        "BRAF",
        "edge-study",
        genome_nexus_client=client,
        n_permutations=5,
    )

    assert run.summary["total_fusions"] == 3
    assert run.summary["mapped_fusions"] == 1
    assert run.summary["skipped_fusions"] == 2
    assert [event.Sample_id for event in run.events] == ["VALID"]
    assert any(
        "BAD-ROLE" in warning and "could not determine" in warning for warning in run.warnings
    )
    assert any(
        "NO-BREAKPOINT" in warning and "no genomic breakpoint" in warning
        for warning in run.warnings
    )


@pytest.mark.parametrize("skip_reason", ["missing", "malformed", "ambiguous"])
@pytest.mark.parametrize("include_mapped", [True, False])
def test_in_frame_and_retention_summaries_exclude_skipped_in_frame_rows(
    tmp_path, genome_nexus_canonical_transcript_fixture_path, skip_reason, include_mapped
):
    client = _genome_nexus_client(genome_nexus_canonical_transcript_fixture_path)
    mapped_calls = (
        [
            _call("IN-FRAME"),
            _call("OUT-OF-FRAME", event_info="Protein Fusion: out of frame  {KIAA1549:BRAF}"),
        ]
        if include_mapped
        else []
    )
    skipped_call = _call("SKIPPED-IN-FRAME", breakpoint=None)
    if skip_reason == "malformed":
        skipped_call["site2Position"] = "not-a-coordinate"
    elif skip_reason == "ambiguous":
        skipped_call.update(site1Position=140493152, site2Position=140493153)
    baseline = analyze_structural_variant_calls(
        mapped_calls, "BRAF", "edge-study", genome_nexus_client=client, n_permutations=5
    )
    run = analyze_structural_variant_calls(
        [*mapped_calls, skipped_call],
        "BRAF",
        "edge-study",
        genome_nexus_client=client,
        n_permutations=5,
    )

    assert run.summary["total_fusions"] == len(mapped_calls) + 1
    assert run.summary["mapped_fusions"] == len(mapped_calls)
    assert run.summary["skipped_fusions"] == 1
    assert any("SKIPPED-IN-FRAME" in warning for warning in run.warnings)
    assert run.rows == baseline.rows
    assert run.summary["in_frame_count"] == (1 if include_mapped else 0)
    assert run.summary["in_frame_percent"] == (50.0 if include_mapped else 0.0)
    assert run.summary["kinase_retained_count"] == len(mapped_calls)
    assert run.summary["kinase_retained_percent"] == (100.0 if include_mapped else 0.0)
    for key in (
        "in_frame_kinase_retained_count",
        "frame_domain_contingency_table",
        "fisher_p_value",
        "fisher_odds_ratio",
        "permutation_p_value",
    ):
        assert json.dumps(run.summary[key]) == json.dumps(baseline.summary[key])

    paths = write_outputs(run, tmp_path, pdf=False)
    report = paths["markdown"].read_text()
    if include_mapped:
        assert "- In-frame: 1/2 (50.0%)" in report
        abstract = render_abstract(json.loads(paths["json"].read_text()))
        assert "1/2 fusions (50.0%) were in-frame" in abstract
        assert "2/2 fusions (100.0%) retained" in abstract


def test_tcga_fusion_annotation_uses_shared_benchmark_pipeline(
    genome_nexus_canonical_transcript_fixture_path,
):
    client = _genome_nexus_client(genome_nexus_canonical_transcript_fixture_path)
    calls = [
        {
            "sampleId": "TCGA-DE-A0Y2-01",
            "site1HugoSymbol": "MACF1",
            "site1Chromosome": "1",
            "site1Position": 39430908,
            "site2HugoSymbol": "BRAF",
            "site2Chromosome": "7",
            "site2Position": 140787584,
            "site2EffectOnFrame": "in-frame",
            "eventInfo": "MACF1-BRAF Fusion",
        }
    ]

    run = analyze_structural_variant_calls(
        calls,
        "BRAF",
        "thca_tcga_pan_can_atlas_2018",
        genome_nexus_client=client,
        n_permutations=5,
    )

    assert run.summary["total_fusions"] == 1
    assert run.summary["mapped_fusions"] == 1
    assert run.events[0].Three_prime_gene == "BRAF"
    assert run.events[0].Frame_status == "in-frame"


def test_real_pipeline_derives_target_exon_and_measures_its_retention(
    genome_nexus_canonical_transcript_fixture_path,
):
    client = _genome_nexus_client(genome_nexus_canonical_transcript_fixture_path)

    run = analyze_structural_variant_calls(
        [_call("SAMPLE-1")],
        "BRAF",
        "msk_impact_50k_2026",
        genome_nexus_client=client,
        n_permutations=5,
        algorithm_names=["exon_retention"],
    )

    result = run.results[0]
    assert result.Algorithm == "exon_retention"
    assert result.Summary["target_exon"] == 11
    assert result.Summary["retained_event_count"] == 1
    assert result.Summary["retained_fraction"] == 1.0


def test_zero_in_frame_records_emit_outputs_with_unavailable_statistics(
    tmp_path,
    genome_nexus_canonical_transcript_fixture_path,
):
    client = _genome_nexus_client(genome_nexus_canonical_transcript_fixture_path)
    run = analyze_structural_variant_calls(
        [_call("OUT", event_info="Protein Fusion: out of frame  {KIAA1549:BRAF}")],
        "BRAF",
        "edge-study",
        genome_nexus_client=client,
        n_permutations=5,
    )

    assert run.summary["total_fusions"] == 1
    assert run.summary["in_frame_count"] == 0
    assert run.summary["fisher_p_value"] is None
    assert any("no mapped in-frame" in warning for warning in run.warnings)
    paths = write_outputs(run, tmp_path)
    assert "p=unavailable" in paths["markdown"].read_text()
    assert json.loads(paths["json"].read_text())["summary"]["fisher_p_value"] is None


def test_empty_study_result_is_actionable_and_does_not_query_genome_nexus(tmp_path):
    client = MagicMock(spec=GenomeNexusClient)

    run = analyze_structural_variant_calls(
        [],
        "BRAF",
        "typo-study",
        genome_nexus_client=client,
        n_permutations=5,
    )

    assert run.summary["total_fusions"] == 0
    assert run.summary["fisher_p_value"] is None
    assert any("Verify the gene and study ID" in warning for warning in run.warnings)
    client.fetch_canonical_transcript.assert_not_called()
    paths = write_outputs(run, tmp_path)
    assert all(path.exists() for path in paths.values())
    assert "comparison is not possible" in paths["markdown"].read_text()


def test_cbioportal_network_failure_is_wrapped_with_actionable_context(monkeypatch):
    def fail_fetch(*_args, **_kwargs):
        raise requests.Timeout("timed out")

    monkeypatch.setattr(cbioportal_api, "fetch_structural_variants", fail_fetch)
    # BRAF's real config opts into mutual_exclusivity_targets (see
    # genes/configs/braf.yaml), so run_real_benchmark also live-fetches that
    # comparator's cohort sample list before the structural-variant fetch
    # under test below -- mock it out so this stays a real-network-free test.
    monkeypatch.setattr(
        cbioportal_api, "fetch_sample_list_ids", MagicMock(return_value=["SAMPLE-1"])
    )
    monkeypatch.setattr(cbioportal_api, "fetch_mutations", MagicMock(return_value=[]))

    with pytest.raises(RealBenchmarkNetworkError) as caught:
        run_real_benchmark("BRAF", "some-study")

    message = str(caught.value)
    assert "cBioPortal request failed" in message
    assert "some-study_structural_variants" in message
    assert "Check the study ID" in message


def test_genome_nexus_network_failure_is_wrapped_with_actionable_context():
    client = MagicMock(spec=GenomeNexusClient)
    client.fetch_canonical_transcript.side_effect = requests.ConnectionError("DNS failed")

    with pytest.raises(RealBenchmarkNetworkError) as caught:
        analyze_structural_variant_calls(
            [_call("SAMPLE")],
            "BRAF",
            "some-study",
            genome_nexus_client=client,
            n_permutations=5,
        )

    message = str(caught.value)
    assert "Genome Nexus/domain lookup failed" in message
    assert "Check network access" in message


def test_real_benchmark_click_command_wires_arguments_and_echoes_summary(
    monkeypatch,
    tmp_path,
):
    fake_run = SimpleNamespace(
        gene_symbol="BRAF",
        summary={
            "total_fusions": 7,
            "mapped_fusions": 6,
            "in_frame_count": 5,
            "kinase_retained_count": 4,
            "fisher_p_value": 0.012345,
        },
        warnings=["Skipped one malformed row"],
        is_gene_pair=False,
    )
    run_mock = MagicMock(return_value=fake_run)
    output_paths = {
        "tsv": Path("out.tsv"),
        "json": Path("out.json"),
        "markdown": Path("out.md"),
    }
    write_mock = MagicMock(return_value=output_paths)
    monkeypatch.setattr(cli, "run_real_benchmark", run_mock)
    monkeypatch.setattr(cli, "write_outputs", write_mock)

    result = CliRunner().invoke(
        cli.main,
        [
            "real-benchmark",
            "BRAF",
            "study-id",
            "--output-dir",
            str(tmp_path),
            "--output-stem",
            "custom",
            "--n-permutations",
            "25",
        ],
    )

    assert result.exit_code == 0, result.output
    run_mock.assert_called_once_with("BRAF", "study-id", n_permutations=25)
    write_mock.assert_called_once_with(
        fake_run,
        tmp_path,
        output_stem="custom",
        pdf=True,
        cli_args=[
            "real-benchmark",
            "BRAF",
            "study-id",
            "--output-dir",
            str(tmp_path),
            "--n-permutations",
            "25",
        ],
    )
    assert "Analyzed 7 BRAF fusions; mapped=6, in-frame=5" in result.output
    assert "Fisher p=0.012345" in result.output
    assert "Warning: Skipped one malformed row" in result.output
    assert "markdown: out.md" in result.output


def test_real_benchmark_click_command_no_pdf_flag_disables_pdf_rendering(monkeypatch, tmp_path):
    fake_run = SimpleNamespace(
        gene_symbol="BRAF",
        summary={
            "total_fusions": 1,
            "mapped_fusions": 1,
            "in_frame_count": 1,
            "kinase_retained_count": 1,
            "fisher_p_value": 0.5,
        },
        warnings=[],
        is_gene_pair=False,
    )
    run_mock = MagicMock(return_value=fake_run)
    write_mock = MagicMock(return_value={"markdown": Path("out.md")})
    monkeypatch.setattr(cli, "run_real_benchmark", run_mock)
    monkeypatch.setattr(cli, "write_outputs", write_mock)

    result = CliRunner().invoke(
        cli.main,
        ["real-benchmark", "BRAF", "study-id", "--output-dir", str(tmp_path), "--no-pdf"],
    )

    assert result.exit_code == 0, result.output
    assert write_mock.call_args.kwargs["pdf"] is False


def test_real_benchmark_click_command_renders_network_error_without_traceback(monkeypatch):
    monkeypatch.setattr(
        cli,
        "run_real_benchmark",
        MagicMock(side_effect=RealBenchmarkNetworkError("cBioPortal timed out; retry later")),
    )

    result = CliRunner().invoke(cli.main, ["real-benchmark", "BRAF", "study-id"])

    assert result.exit_code == 1
    assert "Error: cBioPortal timed out; retry later" in result.output
    assert "Traceback" not in result.output


def test_analyze_click_command_runs_registered_orchestrator_path(monkeypatch, tmp_path):
    fake_run = SimpleNamespace(
        gene_symbol="BRAF",
        summary={"total_fusions": 3},
        results=[object(), object(), object()],
        warnings=[],
        is_gene_pair=False,
    )
    analyze_mock = MagicMock(return_value=fake_run)
    write_mock = MagicMock(return_value={"run_directory": Path("runs/example")})
    monkeypatch.setattr(cli, "run_analysis", analyze_mock)
    monkeypatch.setattr(cli, "write_outputs", write_mock)

    result = CliRunner().invoke(
        cli.main,
        ["analyze", "BRAF", "study-id", "--output-dir", str(tmp_path), "--n-permutations", "9"],
    )

    assert result.exit_code == 0, result.output
    analyze_mock.assert_called_once_with("BRAF", "study-id", n_permutations=9)
    assert "with 3 registered algorithms" in result.output
    assert write_mock.call_args.kwargs["cli_args"][0] == "analyze"


def test_run_analysis_requests_every_registered_algorithm(monkeypatch):
    run_mock = MagicMock(return_value=object())
    monkeypatch.setattr(benchmark_module, "list_algorithms", lambda: ["alpha", "beta"])
    monkeypatch.setattr(benchmark_module, "run_real_benchmark", run_mock)

    result = benchmark_module.run_analysis("GENE", "study", n_permutations=7)

    assert result is run_mock.return_value
    run_mock.assert_called_once_with(
        "GENE",
        "study",
        n_permutations=7,
        algorithm_names=["alpha", "beta"],
    )


def test_tcga_study_config_selects_profile_and_grch38_genome_nexus(monkeypatch):
    """BRAF's real curated config now opts into ``mutual_exclusivity_targets``
    (see genes/configs/braf.yaml), so ``run_real_benchmark`` also live-fetches
    that comparator's data for ANY study, not just msk_impact_50k_2026 --
    the cohort-sample-list/mutation fetches used for that must be mocked
    here too, or this test would otherwise attempt a real network call.
    """
    fetched_calls = [{"sampleId": "TCGA-SAMPLE"}]
    monkeypatch.setattr(
        cbioportal_api,
        "fetch_structural_variants",
        MagicMock(return_value=fetched_calls),
    )
    clinical = pd.DataFrame()
    monkeypatch.setattr(
        cbioportal_api, "fetch_sample_tumor_types", MagicMock(return_value=clinical)
    )
    monkeypatch.setattr(
        cbioportal_api,
        "fetch_sample_list_ids",
        MagicMock(return_value=["TCGA-SAMPLE", "TCGA-OTHER"]),
    )
    monkeypatch.setattr(cbioportal_api, "fetch_mutations", MagicMock(return_value=[]))
    fetch_molecular_data_mock = MagicMock(return_value=[])
    monkeypatch.setattr(cbioportal_api, "fetch_molecular_data", fetch_molecular_data_mock)
    client = MagicMock(spec=GenomeNexusClient)
    client_factory = MagicMock(return_value=client)
    analyze_mock = MagicMock(return_value=SimpleNamespace(endpoints=[], warnings=[]))
    monkeypatch.setattr(benchmark_module, "GenomeNexusClient", client_factory)
    monkeypatch.setattr(benchmark_module, "analyze_structural_variant_calls", analyze_mock)

    result = run_real_benchmark("BRAF", "thca_tcga_pan_can_atlas_2018", n_permutations=7)

    assert result is analyze_mock.return_value
    profile_id = "thca_tcga_pan_can_atlas_2018_structural_variants"
    cbioportal_api.fetch_structural_variants.assert_called_once_with([673], [profile_id])
    client_factory.assert_called_once_with(base_url="https://grch38.genomenexus.org")
    cbioportal_api.fetch_sample_list_ids.assert_called_once_with(
        "thca_tcga_pan_can_atlas_2018_all", base_url=cbioportal_api.DEFAULT_BASE_URL, session=None
    )
    cbioportal_api.fetch_mutations.assert_called_once_with(
        [673],
        ["thca_tcga_pan_can_atlas_2018_mutations"],
        base_url=cbioportal_api.DEFAULT_BASE_URL,
        session=None,
    )
    # The expression-association fetch queries the study's configured mRNA
    # z-score profile (distinct from the structural-variant profile above).
    fetch_molecular_data_mock.assert_called_once_with(
        [673],
        "thca_tcga_pan_can_atlas_2018_rna_seq_v2_mrna_median_Zscores",
        sample_list_id="thca_tcga_pan_can_atlas_2018_all",
        base_url=cbioportal_api.DEFAULT_BASE_URL,
        session=None,
    )
    analyze_mock.assert_called_once_with(
        fetched_calls,
        "BRAF",
        "thca_tcga_pan_can_atlas_2018",
        molecular_profile_id=profile_id,
        clinical_df=clinical,
        genome_nexus_client=client,
        n_permutations=7,
        algorithm_names=None,
        algorithm_params={
            "mutation_cooccurrence": {
                "cohort_sample_ids": ["TCGA-SAMPLE", "TCGA-OTHER"],
                "comparator_alterations": [],
                "comparator_availability": {
                    comparator_target_key(
                        load_gene_config("BRAF").mutual_exclusivity_targets[0]
                    ): True
                },
            }
        },
        extra_warnings=[],
    )
    # No expression records were returned, so a warning was appended rather
    # than the run failing -- expression evidence is optional, corroborating.
    assert result.warnings == [
        "No mRNA-expression records were returned for BRAF from "
        "thca_tcga_pan_can_atlas_2018_rna_seq_v2_mrna_median_Zscores; "
        "expression_association analysis was skipped."
    ]


def test_real_benchmark_click_command_explains_unknown_gene_without_traceback():
    result = CliRunner().invoke(cli.main, ["real-benchmark", "SOMEGENE", "study-id"])

    assert result.exit_code == 1
    assert "Error: Unknown gene 'SOMEGENE'" in result.output
    assert "cfh list-genes" in result.output
    assert "Traceback" not in result.output


@pytest.mark.parametrize("clinical_value", ["Glioma", None, float("nan")])
def test_serialized_clinical_annotations_preserve_statistics(
    tmp_path, genome_nexus_canonical_transcript_fixture_path, clinical_value
):
    client = _genome_nexus_client(genome_nexus_canonical_transcript_fixture_path)
    calls = [_call("SAMPLE-1"), _call("MISSING")]
    baseline = analyze_structural_variant_calls(
        calls, "BRAF", "study", genome_nexus_client=client, n_permutations=5
    )
    enriched = analyze_structural_variant_calls(
        calls,
        "BRAF",
        "study",
        genome_nexus_client=client,
        n_permutations=5,
        clinical_df=pd.DataFrame(
            [
                {
                    "Sample_id": "SAMPLE-1",
                    "Tumor_type": clinical_value,
                    "Oncotree_code": "PA",
                    "Patient_id": "MUST-NOT-CHANGE-GROUPING",
                }
            ]
        ),
    )
    before = json.loads(
        write_outputs(baseline, tmp_path, run_id="before", pdf=False)["json"].read_text()
    )
    after = json.loads(
        write_outputs(enriched, tmp_path, run_id="after", pdf=False)["json"].read_text()
    )
    assert after["events"][0]["tumor_type"] == ("Glioma" if clinical_value == "Glioma" else "")
    assert after["events"][0]["oncotree_code"] == "PA"
    assert after["events"][1]["tumor_type"] is None
    assert after["events"][1]["oncotree_code"] is None
    assert enriched.events[0].Tumor_type == after["events"][0]["tumor_type"]
    for payload in (before, after):
        payload.pop("retrieved_at")
        for result in payload["algorithm_results"]:
            result.pop("Created_at")
            result.pop("Input_fingerprint")  # includes the newly populated annotations
            result["Summary"].pop("Runtime_seconds", None)
        for row in payload["events"]:
            row.pop("tumor_type")
            row.pop("oncotree_code")
    assert json.dumps(before).encode() == json.dumps(after).encode()
