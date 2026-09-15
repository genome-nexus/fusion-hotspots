"""Tests for the gene_pair (joint-partner) CLI/benchmark path.

Before this, a ``gene_pair`` config like ``eml4-alk.yaml`` could not run
through ``cfh analyze``/``cfh real-benchmark`` at all: the single-gene
``_load_benchmark_config`` path requires ``gene_symbol``/``entrez_gene_id``/
``key_domains``, none of which a ``gene_pair`` config declares by design
(see PR #62, which had to invoke ``JointPartnerMode`` directly against
already-fetched ALK events as a manual workaround). These tests exercise the
new dedicated gene-pair path end-to-end without live network access.
"""

from __future__ import annotations

import json
from unittest.mock import MagicMock

import pandas as pd
import pytest
from click.testing import CliRunner

from cfh import cli
from cfh import real_benchmark as benchmark_module
from cfh.genes.registry import GeneConfig, load_gene_config
from cfh.ingestion import cbioportal_api
from cfh.mapping.genome_nexus_source import GenomeNexusClient
from cfh.real_benchmark import (
    RealBenchmarkInputError,
    _maybe_load_gene_pair_config,
    _partner_component_configs,
    run_gene_pair_benchmark,
    run_real_benchmark,
)


def _genome_nexus_client(fixture_path):
    client = MagicMock(spec=GenomeNexusClient)
    client.fetch_canonical_transcript.return_value = json.loads(fixture_path.read_text())
    return client


def _eml4_alk_call(sample_id, *, partner="EML4", breakpoint=140493152):
    """A protein-fusion SV record with ``partner`` as the 5' gene and ALK
    as the target/3' gene, mirroring the ``{FIVE:THREE}`` transcript-order
    annotation real cBioPortal protein-fusion records carry.
    """
    return {
        "sampleId": sample_id,
        "site1HugoSymbol": partner,
        "site2HugoSymbol": "ALK",
        "site2Position": breakpoint,
        "site2EffectOnFrame": "NA",
        "connectionType": "3to3",
        "eventInfo": f"Protein Fusion: in frame  {{{partner}:ALK}}",
    }


def _tmprss2_erg_call(sample_id, *, partner="TMPRSS2", breakpoint=140493152):
    """A protein-fusion SV record with ``partner`` as the 5' gene and ERG
    as the target/3' gene, mirroring the ``{FIVE:THREE}`` transcript-order
    annotation real cBioPortal protein-fusion records carry.

    ``breakpoint`` reuses the same value as :func:`_eml4_alk_call`: both
    mock ``GenomeNexusClient.fetch_canonical_transcript`` to always return
    the shared BRAF canonical-transcript fixture regardless of which gene
    is requested (see ``_genome_nexus_client``), so the breakpoint must
    fall inside *that* fixture's exon locus for the target-locus validation
    to accept it, whichever gene the test is nominally about.
    """
    return {
        "sampleId": sample_id,
        "site1HugoSymbol": partner,
        "site2HugoSymbol": "ERG",
        "site2Position": breakpoint,
        "site2EffectOnFrame": "NA",
        "connectionType": "3to5",
        "eventInfo": f"Protein Fusion: in frame  {{{partner}:ERG}}",
    }


def test_maybe_load_gene_pair_config_detects_pair_configs_only():
    assert _maybe_load_gene_pair_config("eml4-alk") is not None
    assert _maybe_load_gene_pair_config("EML4-ALK") is not None
    assert _maybe_load_gene_pair_config("tmprss2-erg") is not None
    assert _maybe_load_gene_pair_config("TMPRSS2-ERG") is not None
    assert _maybe_load_gene_pair_config("ALK") is None
    assert _maybe_load_gene_pair_config("ERG") is None
    assert _maybe_load_gene_pair_config("NOT-A-REAL-GENE") is None


def test_single_gene_configs_are_unaffected_by_gene_pair_detection():
    """BRAF/RET/ALK/NTRK1 must still be treated as ordinary single-gene
    configs -- never routed through the new gene_pair path -- so their
    existing ``real-benchmark``/``analyze`` behavior is untouched.
    """
    for symbol in ("BRAF", "RET", "ALK", "NTRK1"):
        config = load_gene_config(symbol)
        assert config.gene_pair is None
        assert config.gene_symbol == symbol
        assert config.entrez_gene_id is not None
        assert config.canonical_transcript_id is not None
        assert config.protein_id is not None
        assert _maybe_load_gene_pair_config(symbol) is None


def test_partner_component_configs_uses_only_the_curated_partner():
    pair_config = load_gene_config("eml4-alk")

    configs = _partner_component_configs(pair_config)

    assert [config.gene_symbol for config in configs] == ["ALK"]


def test_partner_component_configs_uses_only_the_curated_partner_for_tmprss2_erg():
    """TMPRSS2-ERG mirrors EML4-ALK's curation shape: only the 3' partner
    (ERG, the recurrently-hub gene fused to several 5' promoter donors, the
    same role ALK plays for EML4-ALK) has a curated single-gene config."""
    pair_config = load_gene_config("tmprss2-erg")

    configs = _partner_component_configs(pair_config)

    assert [config.gene_symbol for config in configs] == ["ERG"]


def test_partner_component_configs_raises_actionable_error_with_no_curated_partner():
    pair_config = GeneConfig(gene_pair=("FAKE5PRIME", "FAKE3PRIME"))

    with pytest.raises(RealBenchmarkInputError, match="no partner with a curated"):
        _partner_component_configs(pair_config)


def test_gene_pair_benchmark_pools_component_events_and_computes_enrichment(
    genome_nexus_canonical_transcript_fixture_path, monkeypatch
):
    client = _genome_nexus_client(genome_nexus_canonical_transcript_fixture_path)
    monkeypatch.setattr(benchmark_module, "GenomeNexusClient", MagicMock(return_value=client))
    calls = [_eml4_alk_call(f"EML4-{i}") for i in range(8)] + [
        _eml4_alk_call(f"OTHER-{i}", partner="OTHERGENE") for i in range(2)
    ]
    monkeypatch.setattr(cbioportal_api, "fetch_structural_variants", MagicMock(return_value=calls))
    monkeypatch.setattr(
        cbioportal_api, "fetch_sample_tumor_types", MagicMock(return_value=pd.DataFrame())
    )

    pair_config = load_gene_config("eml4-alk")
    run = run_gene_pair_benchmark(pair_config, "msk_impact_50k_2026", n_permutations=5)

    assert run.is_gene_pair is True
    assert run.gene_symbol == "EML4-ALK"
    assert run.summary["gene_pair"] == ["EML4", "ALK"]
    assert run.summary["component_genes"] == ["ALK"]
    assert run.summary["eligible_event_count"] == 10
    assert run.summary["observed_count"] == 8
    assert run.summary["fisher_p_value"] is not None
    assert len(run.events) == 10
    # Events are deduplicated by Event_id when pooling component runs.
    assert len({event.Event_id for event in run.events}) == 10
    assert run.results[0].Algorithm == "joint_partner"
    assert run.summary["mechanism_note"] == pair_config.mechanism_note
    assert "dimerization" in run.summary["mechanism_note"].lower()


def test_gene_pair_benchmark_pools_component_events_for_tmprss2_erg_and_surfaces_mechanism_note(
    genome_nexus_canonical_transcript_fixture_path, monkeypatch
):
    """TMPRSS2-ERG mirrors the EML4-ALK pooling/enrichment mechanics exactly
    -- JointPartnerMode's co-occurrence machinery is mechanism-agnostic --
    but its curated ``mechanism_note`` describes a genuinely different
    mechanism (promoter-swap/expression-driven, not domain-retention),
    since the report must state whichever mechanism this pair's own config
    actually curates rather than defaulting to domain-retention language."""
    client = _genome_nexus_client(genome_nexus_canonical_transcript_fixture_path)
    monkeypatch.setattr(benchmark_module, "GenomeNexusClient", MagicMock(return_value=client))
    calls = [_tmprss2_erg_call(f"TMPRSS2-{i}") for i in range(8)] + [
        _tmprss2_erg_call(f"OTHER-{i}", partner="OTHERGENE") for i in range(2)
    ]
    monkeypatch.setattr(cbioportal_api, "fetch_structural_variants", MagicMock(return_value=calls))
    monkeypatch.setattr(
        cbioportal_api, "fetch_sample_tumor_types", MagicMock(return_value=pd.DataFrame())
    )

    pair_config = load_gene_config("tmprss2-erg")
    run = run_gene_pair_benchmark(pair_config, "msk_impact_50k_2026", n_permutations=5)

    assert run.is_gene_pair is True
    assert run.gene_symbol == "TMPRSS2-ERG"
    assert run.summary["gene_pair"] == ["TMPRSS2", "ERG"]
    assert run.summary["component_genes"] == ["ERG"]
    assert run.summary["eligible_event_count"] == 10
    assert run.summary["observed_count"] == 8
    assert run.summary["fisher_p_value"] is not None
    assert run.results[0].Algorithm == "joint_partner"

    assert run.summary["mechanism_note"] == pair_config.mechanism_note
    assert "promoter-swap" in run.summary["mechanism_note"].lower()

    report = benchmark_module._gene_pair_markdown_summary(run)
    assert "## Mechanism" in report
    assert pair_config.mechanism_note in report


def test_run_real_benchmark_routes_a_gene_pair_symbol_automatically(
    genome_nexus_canonical_transcript_fixture_path, monkeypatch
):
    """This is the exact CLI/library gap PR #62 flagged: calling the
    single-gene entry point with a gene_pair-named argument must now
    succeed via JointPartnerMode instead of raising.
    """
    client = _genome_nexus_client(genome_nexus_canonical_transcript_fixture_path)
    monkeypatch.setattr(benchmark_module, "GenomeNexusClient", MagicMock(return_value=client))
    calls = [_eml4_alk_call(f"S{i}") for i in range(5)]
    monkeypatch.setattr(cbioportal_api, "fetch_structural_variants", MagicMock(return_value=calls))
    monkeypatch.setattr(
        cbioportal_api, "fetch_sample_tumor_types", MagicMock(return_value=pd.DataFrame())
    )

    run = run_real_benchmark("EML4-ALK", "msk_impact_50k_2026", n_permutations=5)

    assert run.is_gene_pair is True
    assert run.summary["observed_count"] == 5


def test_gene_pair_cli_analyze_writes_a_run_directory_instead_of_erroring(
    tmp_path, genome_nexus_canonical_transcript_fixture_path, monkeypatch
):
    client = _genome_nexus_client(genome_nexus_canonical_transcript_fixture_path)
    monkeypatch.setattr(benchmark_module, "GenomeNexusClient", MagicMock(return_value=client))
    calls = [_eml4_alk_call(f"S{i}") for i in range(6)] + [
        _eml4_alk_call("O1", partner="OTHERGENE")
    ]
    monkeypatch.setattr(cbioportal_api, "fetch_structural_variants", MagicMock(return_value=calls))
    monkeypatch.setattr(
        cbioportal_api, "fetch_sample_tumor_types", MagicMock(return_value=pd.DataFrame())
    )

    result = CliRunner().invoke(
        cli.main,
        [
            "analyze",
            "EML4-ALK",
            "msk_impact_50k_2026",
            "--output-dir",
            str(tmp_path),
            "--n-permutations",
            "5",
            "--no-pdf",
        ],
    )

    assert result.exit_code == 0, result.output
    assert "observed EML4->ALK=6" in result.output
    run_dirs = sorted(tmp_path.glob("eml4-alk_msk-impact-50k-2026_*"))
    assert len(run_dirs) == 1
    run_dir = run_dirs[0]
    payload = json.loads((run_dir / "results.json").read_text())
    assert payload["summary"]["gene_pair"] == ["EML4", "ALK"]
    assert payload["summary"]["observed_count"] == 6
    assert payload["summary"]["eligible_event_count"] == 7
    assert (run_dir / "report.md").exists()
    assert (run_dir / "results.tsv").exists()
    assert (run_dir / "manifest.json").exists()
    assert not (run_dir / "report.pdf").exists()
    manifest = json.loads((run_dir / "manifest.json").read_text())
    assert manifest["gene"] == "EML4-ALK"
    assert manifest["study_id"] == "msk_impact_50k_2026"


def test_gene_pair_cli_real_benchmark_also_routes_through_joint_partner(
    tmp_path, genome_nexus_canonical_transcript_fixture_path, monkeypatch
):
    client = _genome_nexus_client(genome_nexus_canonical_transcript_fixture_path)
    monkeypatch.setattr(benchmark_module, "GenomeNexusClient", MagicMock(return_value=client))
    calls = [_eml4_alk_call(f"S{i}") for i in range(4)]
    monkeypatch.setattr(cbioportal_api, "fetch_structural_variants", MagicMock(return_value=calls))
    monkeypatch.setattr(
        cbioportal_api, "fetch_sample_tumor_types", MagicMock(return_value=pd.DataFrame())
    )

    result = CliRunner().invoke(
        cli.main,
        [
            "real-benchmark",
            "EML4-ALK",
            "msk_impact_50k_2026",
            "--output-dir",
            str(tmp_path),
            "--n-permutations",
            "5",
            "--no-pdf",
        ],
    )

    assert result.exit_code == 0, result.output
    assert "pair-eligible" in result.output
    run_dirs = list(tmp_path.glob("eml4-alk_msk-impact-50k-2026_*"))
    assert len(run_dirs) == 1
