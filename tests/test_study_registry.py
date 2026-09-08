import pytest

from cfh.studies.registry import StudyConfig, load_study_config


def test_tcga_pan_cancer_atlas_studies_use_grch38_genome_nexus():
    config = load_study_config("thca_tcga_pan_can_atlas_2018")

    assert isinstance(config, StudyConfig)
    assert len(config.study_ids) == 32
    assert config.genome_nexus_base_url == "https://grch38.genomenexus.org"
    assert (
        config.molecular_profile_id("thca_tcga_pan_can_atlas_2018")
        == "thca_tcga_pan_can_atlas_2018_structural_variants"
    )


def test_unconfigured_study_uses_pipeline_defaults():
    assert load_study_config("msk_impact_50k_2026") is None


def test_study_config_rejects_profile_request_for_unlisted_study():
    config = load_study_config("thca_tcga_pan_can_atlas_2018")

    with pytest.raises(ValueError, match="not covered"):
        config.molecular_profile_id("not_in_config")


def test_tcga_pan_cancer_atlas_config_resolves_mrna_expression_profile():
    config = load_study_config("thca_tcga_pan_can_atlas_2018")

    assert (
        config.mrna_expression_profile_id("thca_tcga_pan_can_atlas_2018")
        == "thca_tcga_pan_can_atlas_2018_rna_seq_v2_mrna_median_Zscores"
    )


def test_mrna_expression_profile_id_is_none_without_a_configured_template():
    config = StudyConfig(study_ids=["some_study"])

    assert config.mrna_expression_profile_id("some_study") is None


def test_mrna_expression_profile_id_rejects_an_uncovered_study():
    config = load_study_config("thca_tcga_pan_can_atlas_2018")

    with pytest.raises(ValueError, match="not covered"):
        config.mrna_expression_profile_id("msk_impact_50k_2026")


def test_msk_impact_has_no_mrna_expression_profile_available():
    """msk_impact_50k_2026 is a targeted DNA panel with no mRNA-expression
    molecular profile at all, so it deliberately has no StudyConfig entry --
    load_study_config falls back to None, the same graceful-skip signal the
    expression_association wiring in real_benchmark.py relies on."""
    assert load_study_config("msk_impact_50k_2026") is None
