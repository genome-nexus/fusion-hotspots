"""Held-out ranking comparisons require actual labels and disjoint patients."""

import copy

import pytest

from cfh.stats.independent_validation import compare_heldout_rankings


@pytest.fixture
def inputs():
    discovery = {"patient_id_namespace": "registry-v1", "patient_ids": ["D1", "D2"]}
    validation = {
        "patient_id_namespace": "registry-v1",
        "patient_ids": ["V1", "V2"],
        "candidates": [
            {"candidate_id": "A", "label": 1, "composite_score": 0.8, "recurrence_score": 0.5},
            {"candidate_id": "B", "label": 0, "composite_score": 0.8, "recurrence_score": 0.5},
            {"candidate_id": "C", "label": 0, "composite_score": 0.1, "recurrence_score": 0.2},
        ],
    }
    return discovery, validation


def test_ties_use_group_average_and_fractional_precision_at_k(inputs):
    discovery, validation = inputs
    report = compare_heldout_rankings(discovery, validation, k=1)
    assert report["composite"]["average_precision"] == pytest.approx(0.5)
    assert report["composite"]["precision_at_k"] == pytest.approx(0.5)
    assert report["composite"]["tie_groups"] == 1
    assert report["patient_overlap_count"] == 0
    assert report["delta_average_precision"] == pytest.approx(0)


def test_rejects_overlap_or_unknown_namespace(inputs):
    discovery, validation = inputs
    validation["patient_ids"] = ["D1", "V2"]
    with pytest.raises(ValueError, match="share 1 patient"):
        compare_heldout_rankings(discovery, validation)
    validation["patient_ids"] = ["V1", "V2"]
    validation.pop("patient_id_namespace")
    with pytest.raises(ValueError, match="shared, nonempty"):
        compare_heldout_rankings(discovery, validation)
    validation["patient_id_namespace"] = "other-namespace"
    with pytest.raises(ValueError, match="shared, nonempty"):
        compare_heldout_rankings(discovery, validation)


def test_rejects_missing_scores_or_labels(inputs):
    discovery, validation = inputs
    validation = copy.deepcopy(validation)
    validation["candidates"][0].pop("composite_score")
    with pytest.raises(ValueError, match="finite composite_score"):
        compare_heldout_rankings(discovery, validation)
    validation["candidates"][0]["composite_score"] = 0.8
    validation["candidates"][0].pop("label")
    with pytest.raises(ValueError, match="external binary label"):
        compare_heldout_rankings(discovery, validation)


def test_composite_and_recurrence_compared_on_identical_candidates(inputs):
    discovery, validation = inputs
    validation["candidates"][0]["composite_score"] = 1.0
    validation["candidates"][0]["recurrence_score"] = 0.1
    validation["candidates"][1]["recurrence_score"] = 0.9
    validation["candidates"][2]["recurrence_score"] = 0.8
    report = compare_heldout_rankings(discovery, validation, k=1)
    assert report["composite"]["average_precision"] == 1.0
    assert report["recurrence"]["average_precision"] == pytest.approx(1 / 3)
    assert report["delta_average_precision"] == pytest.approx(2 / 3)
    assert report["delta_precision_at_k"] == 1.0
    validation["candidates"].reverse()
    assert report == compare_heldout_rankings(discovery, validation, k=1)


@pytest.mark.parametrize("value", [None, float("nan"), float("inf"), True])
def test_invalid_score_cannot_silently_drop_a_candidate(inputs, value):
    discovery, validation = inputs
    validation["candidates"][0]["recurrence_score"] = value
    with pytest.raises(ValueError, match="finite recurrence_score"):
        compare_heldout_rankings(discovery, validation)


def test_missing_patient_ids_and_malformed_candidate_fail_closed(inputs):
    discovery, validation = inputs
    discovery["patient_ids"] = []
    with pytest.raises(ValueError, match="complete patient_ids"):
        compare_heldout_rankings(discovery, validation)
    discovery["patient_ids"] = ["D1"]
    validation["candidates"] = [None]
    with pytest.raises(ValueError, match="must be an object"):
        compare_heldout_rankings(discovery, validation)
