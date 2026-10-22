"""Repeated observations of one fusion in one patient are counted once."""

from __future__ import annotations

import pytest

from cfh.model.fusion_event import FusionEvent
from cfh.model.fusion_feature import FusionFeature
from cfh.stats.observation_units import collapse_repeated_observations


def _pair(event_id, *, patient="P1", sample=None, partner="KIAA1549", aa=438, frame="in-frame"):
    event = FusionEvent(
        Event_id=event_id,
        Cohort="c",
        Patient_id=patient,
        Sample_id=sample or f"{patient}-{event_id}",
        Site1_gene=partner,
        Site2_gene="BRAF",
        Frame_status=frame,
    )
    feature = FusionFeature(
        Event_id=event_id, Gene="BRAF", Role="three_prime", Junction_position_aa=aa
    )
    return event, feature


def _collapse(pairs, **kwargs):
    events, features = map(list, zip(*pairs, strict=True))
    return collapse_repeated_observations(events, features, "BRAF", **kwargs)


_REPEATS_AND_DISTINCT_FUSIONS = [
    ("E1", {}),
    ("E2", {}),  # same patient, partner, junction: repeat observation
    ("E3", {"partner": "AGAP3"}),  # different partner in same patient
    ("E4", {"aa": 500}),  # different junction in same patient
    ("E5", {"patient": "P2"}),
]


def test_patient_unit_collapses_repeats_but_keeps_distinct_fusions():
    events, features, report = _collapse(
        [_pair(event_id, **kwargs) for event_id, kwargs in _REPEATS_AND_DISTINCT_FUSIONS],
        observation_unit="patient",
    )
    assert [event.Event_id for event in events] == ["E1", "E3", "E4", "E5"]
    assert [feature.Event_id for feature in features] == ["E1", "E3", "E4", "E5"]
    assert report.collapsed_event_count == 1
    assert report.patients_with_repeats == 1
    assert report.as_summary()["analyzed_observation_count"] == 4
    assert report.merged_distinct_fusion_count == 0


def test_default_patient_gene_unit_counts_each_patient_once():
    events, features, report = _collapse(
        [_pair(event_id, **kwargs) for event_id, kwargs in _REPEATS_AND_DISTINCT_FUSIONS]
    )
    assert [event.Event_id for event in events] == ["E1", "E5"]
    assert [feature.Event_id for feature in features] == ["E1", "E5"]
    summary = report.as_summary()
    assert summary["observation_unit"] == "patient_gene"
    assert summary["collapsed_repeat_observation_count"] == 1
    assert summary["merged_distinct_fusion_count"] == 2
    assert summary["patients_with_multiple_distinct_fusions"] == 1
    assert summary["analyzed_observation_count"] == 2


def test_patient_gene_representative_prefers_a_determinate_frame():
    events, _, _ = _collapse(
        [_pair("E1", frame="unknown"), _pair("E2", partner="AGAP3", frame="out-of-frame")]
    )
    assert [event.Event_id for event in events] == ["E2"]


def test_missing_patient_falls_back_to_sample_identity():
    events, _, report = _collapse(
        [
            _pair("E1", patient=None, sample="S1"),
            _pair("E2", patient=None, sample="S1"),
            _pair("E3", patient=None, sample="S2"),
        ]
    )
    assert [event.Event_id for event in events] == ["E1", "E3"]
    assert report.identity_fallback_counts == {"sample": 3}


def test_representative_prefers_determinate_frame_and_counts_conflicts():
    events, _, report = _collapse([_pair("E1", frame="unknown"), _pair("E2", frame="in-frame")])
    assert [event.Event_id for event in events] == ["E2"]
    assert report.frame_conflict_groups == 0
    _, _, conflicted = _collapse([_pair("E1", frame="out-of-frame"), _pair("E2", frame="in-frame")])
    assert conflicted.frame_conflict_groups == 1


def test_event_unit_is_identity_and_bad_inputs_raise():
    pairs = [_pair("E1"), _pair("E2")]
    events, features, report = _collapse(pairs, observation_unit="event")
    assert len(events) == len(features) == 2
    assert report.collapsed_event_count == 0
    with pytest.raises(ValueError):
        _collapse(pairs, observation_unit="sample")
    with pytest.raises(ValueError):
        collapse_repeated_observations([pairs[0][0]], [pairs[1][1]], "BRAF")
