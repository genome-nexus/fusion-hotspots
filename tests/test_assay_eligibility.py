"""Assay coverage is evidence for the denominator, never inferred from zero calls."""

from unittest.mock import MagicMock

import pytest
import requests

from cfh.algorithms.mutation_cooccurrence import (
    MutationCooccurrenceAlgorithm,
    comparator_target_key,
)
from cfh.genes.registry import load_gene_config
from cfh.ingestion import cbioportal_api
from cfh.model.fusion_event import FusionEvent
from cfh.real_benchmark import _fetch_mutual_exclusivity_params


def test_panel_membership_profile_participation_and_unknown_coverage():
    session = MagicMock()
    session.post.return_value.json.return_value = [
        {"sampleId": "covered", "profiled": True, "genePanelId": "A"},
        {"sampleId": "replicate", "profiled": True, "genePanelId": "A"},
        {"sampleId": "off_panel", "profiled": True, "genePanelId": "B"},
        {"sampleId": "not_profiled", "profiled": False, "genePanelId": "A"},
        {"sampleId": "unknown", "profiled": True},
    ]
    session.get.side_effect = [
        MagicMock(json=MagicMock(return_value={"genes": [{"entrezGeneId": 673}]})),
        MagicMock(json=MagicMock(return_value={"genes": [{"entrezGeneId": 1}]})),
    ]
    assert cbioportal_api.fetch_gene_panel_eligibility(
        "study_sv", "study_all", 673, base_url="https://portal.example/api", session=session
    ) == {
        "covered": True,
        "replicate": True,
        "off_panel": False,
        "not_profiled": False,
        "unknown": None,
    }
    assert session.get.call_count == 2
    session.post.assert_called_once_with(
        "https://portal.example/api/molecular-profiles/study_sv/gene-panel-data/fetch",
        json={"sampleListId": "study_all"},
        timeout=30,
    )


@pytest.mark.parametrize(
    "payload", [None, [{}], [{"sampleId": "S", "profiled": True, "genePanelId": "A"}]]
)
def test_malformed_coverage_cannot_become_confirmed_negative(payload):
    session = MagicMock()
    session.post.return_value.json.return_value = payload
    session.get.return_value.json.return_value = {"genes": None}
    with pytest.raises(ValueError):
        cbioportal_api.fetch_gene_panel_eligibility("sv", "all", 673, session=session)


def test_joint_coverage_excludes_unmeasured_calls_and_preserves_denominator(monkeypatch):
    config = load_gene_config("BRAF")
    target_key = comparator_target_key(config.mutual_exclusivity_targets[0])
    monkeypatch.setattr(
        cbioportal_api, "fetch_sample_list_ids", lambda *a, **k: ["A", "B", "C", "D"]
    )
    monkeypatch.setattr(
        cbioportal_api,
        "fetch_mutations",
        lambda *a, **k: [
            {"sampleId": "B", "proteinChange": "V600E"},
            {"sampleId": "D", "proteinChange": "V600E"},
            {"sampleId": "OUTSIDE", "proteinChange": "V600E"},
        ],
    )
    monkeypatch.setattr(
        cbioportal_api,
        "fetch_gene_panel_eligibility",
        MagicMock(
            side_effect=[
                {"A": True, "B": True, "C": True, "D": False},
                {"A": True, "B": True, "C": None, "D": True},
            ]
        ),
    )
    params, warnings = _fetch_mutual_exclusivity_params(config, "study", None)
    values = params["mutation_cooccurrence"]
    assert values["eligible_sample_ids_by_target"][target_key] == ["A", "B"]
    result = MutationCooccurrenceAlgorithm().run(
        [FusionEvent(Event_id=s, Sample_id=s, Cohort="study") for s in ["A", "C", "OUTSIDE"]],
        [],
        config,
        values,
    )
    row = result.Summary["targets"][0]
    assert row["contingency_table"] == [[0, 1], [1, 0]]
    assert row["cohort_sample_count"] == sum(map(sum, row["contingency_table"])) == 2
    assert row["excluded_fusion_positive_sample_count"] == 2
    assert row["excluded_comparator_sample_count"] == 2
    assert warnings


def test_failed_eligibility_skips_available_zero_call_comparator(monkeypatch):
    config = load_gene_config("BRAF")
    monkeypatch.setattr(cbioportal_api, "fetch_sample_list_ids", lambda *a, **k: ["A"])
    monkeypatch.setattr(cbioportal_api, "fetch_mutations", lambda *a, **k: [])
    monkeypatch.setattr(
        cbioportal_api,
        "fetch_gene_panel_eligibility",
        MagicMock(side_effect=requests.Timeout("offline")),
    )
    params, warnings = _fetch_mutual_exclusivity_params(config, "study", None)
    result = MutationCooccurrenceAlgorithm().run([], [], config, params["mutation_cooccurrence"])
    assert not result.Summary["targets"]
    assert warnings


def test_patient_identity_preserved_in_annotations():
    session = MagicMock()
    session.post.return_value.json.return_value = [
        {"sampleId": "S", "patientId": "P", "clinicalAttributeId": "CANCER_TYPE", "value": "Glioma"}
    ]
    assert (
        cbioportal_api.fetch_sample_tumor_types("study", ["S"], session=session).iloc[0][
            "Patient_id"
        ]
        == "P"
    )


def test_expression_background_requires_sv_coverage(monkeypatch):
    from cfh.real_benchmark import _fetch_expression_association_params
    from cfh.studies.registry import StudyConfig

    config = load_gene_config("BRAF")
    study = StudyConfig(study_ids=["study"], mrna_expression_profile_template="{study_id}_rna")
    monkeypatch.setattr(
        cbioportal_api,
        "fetch_molecular_data",
        lambda *a, **k: [{"sampleId": "A", "value": 1.0}, {"sampleId": "B", "value": 2.0}],
    )
    monkeypatch.setattr(
        cbioportal_api, "fetch_gene_panel_eligibility", lambda *a, **k: {"A": True, "B": None}
    )
    params, warning = _fetch_expression_association_params(config, study, "study")
    assert params["cohort_sample_ids"] == ["A"]
    assert set(params["expression_by_sample"]) == {"A", "B"}
    assert "1 measured samples excluded" in warning
