"""Environment-gated live checks for expression_association/mutation_cooccurrence
coverage across the curated genes beyond BRAF/RET (see PR description for the
full gene-by-gene live-data investigation this locks in).

Mirrors ``test_ret_tcga_live_generalization.py``'s pattern: skipped unless
``CFH_RUN_NETWORK_TESTS=1`` (the same scheduled ``live-data-check`` workflow
that already exercises the BRAF/RET live paths), and asserts on structural
facts (mapped-fusion counts, which comparisons ran) rather than exact
p-values, since those can shift slightly as cBioPortal's underlying data is
revised.
"""

from __future__ import annotations

import os

import pytest

from cfh.ingestion import cbioportal_api
from cfh.mapping.genome_nexus_source import GenomeNexusGeneNotFound
from cfh.real_benchmark import run_analysis
from cfh.studies.registry import load_study_config

pytestmark = [
    pytest.mark.network,
    pytest.mark.skipif(
        os.getenv("CFH_RUN_NETWORK_TESTS") != "1",
        reason="set CFH_RUN_NETWORK_TESTS=1 to run live expression/mutation-cooccurrence "
        "benchmarks",
    ),
]


def _algorithm_result(run, name: str):
    return next(result for result in run.results if result.Algorithm == name)


def test_alk_expression_association_against_real_luad_tcga_sv():
    run = run_analysis("ALK", "luad_tcga_pan_can_atlas_2018", n_permutations=25)
    assert run.summary["mapped_fusions"] == 5
    expression = _algorithm_result(run, "expression_association")
    comparison = expression.Summary["fusion_positive_vs_negative"]
    assert comparison["n_fusion_positive"] == 5
    assert comparison["p_value"] < 0.01


def test_ntrk1_expression_association_against_real_thca_tcga_sv():
    run = run_analysis("NTRK1", "thca_tcga_pan_can_atlas_2018", n_permutations=25)
    assert run.summary["mapped_fusions"] == 6
    expression = _algorithm_result(run, "expression_association")
    comparison = expression.Summary["fusion_positive_vs_negative"]
    assert comparison["n_fusion_positive"] == 6
    assert comparison["p_value"] < 0.01


def test_etv6_expression_association_against_real_thca_tcga_sv():
    run = run_analysis("ETV6", "thca_tcga_pan_can_atlas_2018", n_permutations=25)
    assert run.summary["mapped_fusions"] == 6
    expression = _algorithm_result(run, "expression_association")
    comparison = expression.Summary["fusion_positive_vs_negative"]
    assert comparison["n_fusion_positive"] == 5


def test_ros1_expression_association_against_real_luad_tcga_sv():
    run = run_analysis("ROS1", "luad_tcga_pan_can_atlas_2018", n_permutations=25)
    assert run.summary["mapped_fusions"] == 10
    expression = _algorithm_result(run, "expression_association")
    comparison = expression.Summary["fusion_positive_vs_negative"]
    assert comparison["n_fusion_positive"] == 7
    assert comparison["p_value"] < 0.01


def test_fgfr2_mutation_cooccurrence_against_real_brca_tcga_data():
    """FGFR2 S252W (Pollock et al. 2007, PMID 17297457) is a real,
    literature-established recurrent hotspot -- but live cBioPortal data
    confirms it does not co-occur with any FGFR2 fusion in this cohort;
    this is a real reported null result, not a skipped comparison."""
    run = run_analysis("FGFR2", "brca_tcga_pan_can_atlas_2018", n_permutations=25)
    assert run.summary["mapped_fusions"] == 5
    cooccurrence = _algorithm_result(run, "mutation_cooccurrence")
    target = cooccurrence.Summary["targets"][0]
    assert target["protein_change"] == "S252W"
    assert target["comparator_altered_sample_count"] == 0
    expression = _algorithm_result(run, "expression_association")
    comparison = expression.Summary["fusion_positive_vs_negative"]
    assert comparison["n_fusion_positive"] == 5


def test_ntrk3_has_no_genome_nexus_grch38_canonical_transcript():
    """Real, external data-availability gap: NTRK3 fusion events exist in
    TCGA PanCancer Atlas (thca_tcga_pan_can_atlas_2018 carries real
    ETV6-NTRK3/RBPMS-NTRK3 SV records), but the GRCh38 Genome Nexus
    instance every TCGA PanCancer Atlas study is configured against has no
    canonical-transcript mapping for NTRK3 at all, so no domain-retention,
    expression-association, or mutation-cooccurrence result can be produced
    for it against this cohort family. This locks in that this is a live,
    external-service limitation, not a regression, if Genome Nexus later
    adds the mapping this test starts failing and should be revisited."""
    with pytest.raises(GenomeNexusGeneNotFound):
        run_analysis("NTRK3", "thca_tcga_pan_can_atlas_2018", n_permutations=25)


def test_fli1_has_no_fusion_records_anywhere_in_tcga_pan_cancer_atlas():
    """Real gap: EWSR1-FLI1 is an Ewing-sarcoma-defining fusion, and Ewing
    sarcoma is not one of the TCGA PanCancer Atlas's 32 adult-cancer
    cohorts, so FLI1 has zero structural-variant records across the whole
    atlas -- there is no cohort to run expression_association or
    mutation_cooccurrence against."""
    study_ids = load_study_config("thca_tcga_pan_can_atlas_2018").study_ids
    profile_ids = [f"{study_id}_structural_variants" for study_id in study_ids]
    calls = cbioportal_api.fetch_structural_variants([2313], profile_ids)
    assert calls == []
