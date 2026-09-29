"""Live comparator-fetch availability regressions.

These tests deliberately exercise only the additive comparator-fetch helper,
so fetch failure behavior stays isolated from the larger benchmark suite.
"""

from unittest.mock import MagicMock

import pytest
import requests

from cfh.algorithms.mutation_cooccurrence import comparator_target_key
from cfh.genes.registry import GeneConfig, MutualExclusivityTarget
from cfh.ingestion import cbioportal_api
from cfh.real_benchmark import _fetch_mutual_exclusivity_params


@pytest.fixture(autouse=True)
def _profile_eligibility(monkeypatch):
    monkeypatch.setattr(
        cbioportal_api,
        "fetch_gene_panel_eligibility",
        MagicMock(return_value={"S1": True, "S2": True}),
    )


def _config() -> GeneConfig:
    return GeneConfig(
        gene_symbol="FUSION",
        entrez_gene_id=100,
        canonical_transcript_id="NM_000001",
        protein_id="P00001",
        mutual_exclusivity_targets=[
            MutualExclusivityTarget(
                gene="SUCCESS", alteration_type="point_mutation", entrez_gene_id=1
            ),
            MutualExclusivityTarget(
                gene="FAILED", alteration_type="point_mutation", entrez_gene_id=2
            ),
            MutualExclusivityTarget(gene="NO_ID", alteration_type="point_mutation"),
            MutualExclusivityTarget(
                gene="UNKNOWN", alteration_type="unsupported", entrez_gene_id=4
            ),
        ],
    )


def _single_target_config() -> GeneConfig:
    config = _config()
    return config.model_copy(
        update={"mutual_exclusivity_targets": config.mutual_exclusivity_targets[:1]}
    )


@pytest.mark.parametrize("unsupported_type", ["unsupported", "cna_unrecognized"])
def test_live_fetch_marks_each_unavailable_target_and_retains_successful_zero_result(
    monkeypatch, unsupported_type
):
    config = _config()
    config.mutual_exclusivity_targets[-1].alteration_type = unsupported_type
    sample_fetch = MagicMock(return_value=["S1", "S2"])
    mutation_fetch = MagicMock(side_effect=[[], requests.Timeout("timed out")])
    client_session = MagicMock()
    monkeypatch.setattr(cbioportal_api, "fetch_sample_list_ids", sample_fetch)
    monkeypatch.setattr(cbioportal_api, "fetch_mutations", mutation_fetch)

    params, warnings = _fetch_mutual_exclusivity_params(
        config, "study", None, base_url="https://portal.example/api", session=client_session
    )

    assert params is not None
    comparator_params = params["mutation_cooccurrence"]
    assert comparator_params["comparator_alterations"] == []
    assert comparator_params["comparator_availability"] == {
        comparator_target_key(config.mutual_exclusivity_targets[0]): True,
        comparator_target_key(config.mutual_exclusivity_targets[1]): False,
        comparator_target_key(config.mutual_exclusivity_targets[2]): False,
        comparator_target_key(config.mutual_exclusivity_targets[3]): False,
    }
    assert any("Could not fetch or validate point_mutation data for FAILED" in w for w in warnings)
    assert any("NO_ID has no entrez_gene_id" in warning for warning in warnings)
    assert any("Unrecognized mutual_exclusivity_targets" in warning for warning in warnings)
    sample_fetch.assert_called_once_with(
        "study_all", base_url="https://portal.example/api", session=client_session
    )
    assert mutation_fetch.call_count == 2
    assert all(
        call.kwargs["base_url"] == "https://portal.example/api"
        for call in mutation_fetch.call_args_list
    )


def test_cohort_universe_failure_skips_entire_comparator_layer(monkeypatch):
    sample_fetch = MagicMock(side_effect=requests.ConnectionError("offline"))
    mutation_fetch = MagicMock()
    monkeypatch.setattr(cbioportal_api, "fetch_sample_list_ids", sample_fetch)
    monkeypatch.setattr(cbioportal_api, "fetch_mutations", mutation_fetch)

    params, warnings = _fetch_mutual_exclusivity_params(_config(), "study", None)

    assert params is None
    assert any("Could not fetch cohort sample universe" in warning for warning in warnings)
    mutation_fetch.assert_not_called()


def test_none_comparator_response_marks_target_unavailable(monkeypatch):
    config = _single_target_config()
    monkeypatch.setattr(cbioportal_api, "fetch_sample_list_ids", MagicMock(return_value=["S1"]))
    monkeypatch.setattr(cbioportal_api, "fetch_mutations", MagicMock(return_value=None))

    params, warnings = _fetch_mutual_exclusivity_params(config, "study", None)

    assert params is not None
    assert params["mutation_cooccurrence"]["comparator_alterations"] == []
    assert params["mutation_cooccurrence"]["comparator_availability"] == {
        comparator_target_key(config.mutual_exclusivity_targets[0]): False
    }
    assert any("Could not fetch or validate point_mutation data" in warning for warning in warnings)


@pytest.mark.parametrize("calls", [[{}], [{"sampleId": "S1", "proteinChange": "V600E"}, {}]])
def test_malformed_comparator_rows_mark_target_unavailable(monkeypatch, calls):
    config = _single_target_config()
    monkeypatch.setattr(cbioportal_api, "fetch_sample_list_ids", MagicMock(return_value=["S1"]))
    monkeypatch.setattr(cbioportal_api, "fetch_mutations", MagicMock(return_value=calls))

    params, warnings = _fetch_mutual_exclusivity_params(config, "study", None)

    assert params is not None
    assert params["mutation_cooccurrence"]["comparator_alterations"] == []
    assert params["mutation_cooccurrence"]["comparator_availability"] == {
        comparator_target_key(config.mutual_exclusivity_targets[0]): False
    }
    assert any("contained malformed record(s)" in warning for warning in warnings)


@pytest.mark.parametrize("failed", [False, True])
def test_cna_fetch_preserves_availability_and_connection_context(monkeypatch, failed):
    config = _single_target_config()
    target = config.mutual_exclusivity_targets[0].model_copy(update={"alteration_type": "cna_amp"})
    config = config.model_copy(update={"mutual_exclusivity_targets": [target]})
    monkeypatch.setattr(cbioportal_api, "fetch_sample_list_ids", MagicMock(return_value=["S1"]))
    fetch = MagicMock(
        return_value=[{"sampleId": "S1", "alteration": 2}],
        side_effect=requests.Timeout("unavailable") if failed else None,
    )
    session = MagicMock()
    monkeypatch.setattr(cbioportal_api, "fetch_discrete_copy_number", fetch)
    params, warnings = _fetch_mutual_exclusivity_params(
        config, "study", None, base_url="https://custom.example/api", session=session
    )
    values = params["mutation_cooccurrence"]
    assert values["comparator_availability"][comparator_target_key(target)] is not failed
    if failed:
        assert values["comparator_alterations"] == []
        assert warnings
    else:
        assert values["comparator_alterations"][0]["Alteration_type"] == "cna_amp"
        assert warnings == []
    fetch.assert_called_once_with(
        [1], "study_cna", "study_all", base_url="https://custom.example/api", session=session
    )
