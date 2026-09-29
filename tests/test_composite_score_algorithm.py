"""Unit tests for the composite_score algorithm's aggregation math.

All fixtures here are synthetic and fabricated (a made-up ``FAKE1`` gene,
made-up partner genes) specifically to prove the aggregation formula itself,
independent of any real gene's biology. The real-data BRAF/RET proof lives
in ``tests/benchmark/test_composite_score_real_msk_impact.py``.
"""

from __future__ import annotations

import pytest

from cfh.algorithms import registry
from cfh.algorithms.composite_score import (
    DEFAULT_WEIGHTS,
    CompositeScoreAlgorithm,
    _clipped_neg_log10,
)
from cfh.genes.registry import GeneConfig, KeyDomain
from cfh.model.algorithm_result import AlgorithmResult
from cfh.model.fusion_event import FusionEvent
from cfh.model.fusion_feature import FusionFeature

_FAKE_GENE = GeneConfig(
    gene_symbol="FAKE1",
    canonical_transcript_id="NM_000001",
    protein_id="P00001",
    key_domains=[KeyDomain(name="Made-up domain", source="test", key="made_up")],
    disruption_required_domains=[KeyDomain(name="Made-up autoinhib", source="test", key="autoinh")],
)

_FAKE_GENE_NO_DISRUPTION = GeneConfig(
    gene_symbol="FAKE1",
    canonical_transcript_id="NM_000001",
    protein_id="P00001",
    key_domains=[KeyDomain(name="Made-up domain", source="test", key="made_up")],
)


def _frequency_result(counts: dict[str, int]) -> AlgorithmResult:
    table = [{"Partner_gene": partner, "Event_count": count} for partner, count in counts.items()]
    return AlgorithmResult(
        Algorithm="frequency",
        Summary={
            "input_event_count": sum(counts.values()),
            "analyzed_event_count": sum(counts.values()),
            "unique_partner_gene_count": len(counts),
        },
        Tables={"Partner_gene_counts": table},
        Warnings=[],
    )


def _domain_result(
    algorithm: str, *, fisher_p=None, permutation_p=None, warnings=None
) -> AlgorithmResult:
    return AlgorithmResult(
        Algorithm=algorithm,
        Summary={
            "fisher_odds_ratio": None if fisher_p is None else 3.0,
            "fisher_p_value": fisher_p,
            "permutation_empirical_p_value": permutation_p,
            "observed_in_frame_retention_rate": None,
        },
        Tables={},
        Warnings=warnings or [],
    )


def _cutpoint_result(
    *, determinable: bool, inferred_cutpoint_aa=None, corrected_p_value=0.01
) -> AlgorithmResult:
    return AlgorithmResult(
        Algorithm="cutpoint_detection",
        Summary={
            "determinable": determinable,
            "reason": None if determinable else "not enough data",
            "n_events_analyzed": 10,
            "inferred_cutpoint_aa": inferred_cutpoint_aa,
            "observed_statistic_neg_log10_p": None,
            "observed_p_value": None,
            "observed_odds_ratio": None,
            "corrected_p_value": corrected_p_value if determinable else None,
            "known_domain_boundary_comparison": None,
        },
        Tables={},
        Warnings=[],
    )


def _confidence_result(groups: dict[str, dict]) -> AlgorithmResult:
    return AlgorithmResult(
        Algorithm="confidence_stats",
        Summary={
            "group_field": "Frame_status",
            "group_a_label": "True",
            "group_b_label": "False",
            "n_group_a": groups.get("True", {}).get("n", 0),
            "n_group_b": groups.get("False", {}).get("n", 0),
            "mle": {"outcome_field": "made_up", "groups": groups},
        },
        Tables={},
        Warnings=[],
    )


def _events_and_features(partner_positions: dict[str, list[int]]) -> tuple[list, list]:
    """One event+feature pair per breakpoint position, grouped by partner."""
    events = []
    features = []
    counter = 0
    for partner, positions in partner_positions.items():
        for position in positions:
            counter += 1
            event_id = f"evt{counter}"
            events.append(
                FusionEvent(
                    Event_id=event_id,
                    Cohort="synthetic",
                    Site1_gene="FAKE1",
                    Site2_gene=partner,
                )
            )
            features.append(
                FusionFeature(Event_id=event_id, Gene="FAKE1", Junction_position_aa=position)
            )
    return events, features


def test_algorithm_registered():
    assert registry.get("composite_score") is CompositeScoreAlgorithm
    assert "composite_score" in registry.list_algorithms()


@pytest.mark.parametrize("corrected_p", [None, 0.05, 0.8])
def test_cutpoint_component_requires_significant_corrected_p(corrected_p):
    events, features = _events_and_features({"PARTNER_A": [500]})
    result = CompositeScoreAlgorithm().run(
        events,
        features,
        _FAKE_GENE,
        {
            "algorithm_results": [
                _frequency_result({"PARTNER_A": 1}),
                _cutpoint_result(
                    determinable=True,
                    inferred_cutpoint_aa=500,
                    corrected_p_value=corrected_p,
                ),
            ]
        },
    )
    row = result.Tables["composite_evidence_ranking"][0]
    assert row["Cutpoint_proximity_score"] is None
    assert "cutpoint_proximity" not in row["Components_applicable"]
    assert result.Summary["score_interpretation"] == (
        "prioritization_heuristic_not_probability_or_significance"
    )
    assert any("corrected cutpoint p-value" in warning for warning in result.Warnings)


def test_clipped_neg_log10_saturates_and_bounds_to_unit_interval():
    assert _clipped_neg_log10(1.0, cap=10.0) == 0.0
    assert _clipped_neg_log10(0.1, cap=10.0) == pytest.approx(0.1)
    assert _clipped_neg_log10(1e-10, cap=10.0) == pytest.approx(1.0)
    assert _clipped_neg_log10(1e-20, cap=10.0) == 1.0  # saturates, never exceeds 1.0
    assert _clipped_neg_log10(0.0, cap=10.0) == 1.0


def test_requires_frequency_result():
    with pytest.raises(ValueError, match="frequency"):
        CompositeScoreAlgorithm().run([], [], _FAKE_GENE, {"algorithm_results": []})


def test_negative_weight_override_is_rejected():
    """A negative weight breaks the documented [0, 1] composite-score
    guarantee (the weighted average is only bounded when every weight and
    every sub-score is non-negative), so it must be rejected up front.
    """
    events, features = _events_and_features({"PARTNER_A": [500]})
    algorithm_results = [_frequency_result({"PARTNER_A": 1})]

    with pytest.raises(ValueError, match="recurrence"):
        CompositeScoreAlgorithm().run(
            events,
            features,
            _FAKE_GENE_NO_DISRUPTION,
            {"algorithm_results": algorithm_results, "weights": {"recurrence": -5.0}},
        )


def test_non_finite_weight_override_is_rejected():
    events, features = _events_and_features({"PARTNER_A": [500]})
    algorithm_results = [_frequency_result({"PARTNER_A": 1})]

    for bad_weight in (float("nan"), float("inf"), float("-inf")):
        with pytest.raises(ValueError, match="recurrence"):
            CompositeScoreAlgorithm().run(
                events,
                features,
                _FAKE_GENE_NO_DISRUPTION,
                {"algorithm_results": algorithm_results, "weights": {"recurrence": bad_weight}},
            )


def test_non_finite_or_non_positive_neg_log10_p_cap_is_rejected():
    events, features = _events_and_features({"PARTNER_A": [500]})
    algorithm_results = [_frequency_result({"PARTNER_A": 1})]

    for bad_cap in (0.0, -1.0, float("nan"), float("inf")):
        with pytest.raises(ValueError, match="neg_log10_p_cap"):
            CompositeScoreAlgorithm().run(
                events,
                features,
                _FAKE_GENE_NO_DISRUPTION,
                {"algorithm_results": algorithm_results, "neg_log10_p_cap": bad_cap},
            )


def test_all_zero_weights_are_rejected_outright():
    """An all-zero weight override would otherwise make every row's weight
    denominator zero, producing a fabricated 0.0 composite score
    indistinguishable from a genuine "no evidence" score. Reject it as a
    configuration error instead.
    """
    events, features = _events_and_features({"PARTNER_A": [500]})
    algorithm_results = [_frequency_result({"PARTNER_A": 1})]

    with pytest.raises(ValueError, match="zero"):
        CompositeScoreAlgorithm().run(
            events,
            features,
            _FAKE_GENE_NO_DISRUPTION,
            {
                "algorithm_results": algorithm_results,
                "weights": {name: 0.0 for name in DEFAULT_WEIGHTS},
            },
        )


def test_zero_weight_for_a_rows_only_applicable_subscore_yields_none_not_zero():
    """A *valid* (not all-zero) weight set can still leave a specific row's
    own applicable sub-scores all weighted 0.0 -- here, recurrence is the
    row's only applicable component and its weight is overridden to 0.0
    while other weights stay positive. That row's Composite_score must be
    None (with a warning), not a fabricated 0.0, and must sort after every
    row with a real score.
    """
    events, features = _events_and_features({"PARTNER_A": [500], "PARTNER_B": [500]})
    algorithm_results = [_frequency_result({"PARTNER_A": 1, "PARTNER_B": 1})]

    result = CompositeScoreAlgorithm().run(
        events,
        features,
        _FAKE_GENE_NO_DISRUPTION,
        {"algorithm_results": algorithm_results, "weights": {"recurrence": 0.0}},
    )

    ranking = result.Tables["composite_evidence_ranking"]
    assert len(ranking) == 2
    for row in ranking:
        assert row["Composite_score"] is None
    assert any("zero total weight" in warning for warning in result.Warnings)


def test_unrecognized_weight_override_key_is_rejected():
    events, features = _events_and_features({"PARTNER_A": [500]})
    algorithm_results = [_frequency_result({"PARTNER_A": 1})]

    with pytest.raises(ValueError, match="unrecognized"):
        CompositeScoreAlgorithm().run(
            events,
            features,
            _FAKE_GENE_NO_DISRUPTION,
            {"algorithm_results": algorithm_results, "weights": {"recurence": 0.5}},
        )


def test_non_numeric_weight_value_raises_value_error_not_type_error():
    events, features = _events_and_features({"PARTNER_A": [500]})
    algorithm_results = [_frequency_result({"PARTNER_A": 1})]

    with pytest.raises(ValueError, match="recurrence"):
        CompositeScoreAlgorithm().run(
            events,
            features,
            _FAKE_GENE_NO_DISRUPTION,
            {"algorithm_results": algorithm_results, "weights": {"recurrence": "high"}},
        )


def test_composite_score_matches_hand_computed_weighted_average():
    events, features = _events_and_features(
        {
            "PARTNER_A": [499, 500, 501],  # mean distance to cutpoint(500) = 0.6667
            "PARTNER_B": [100, 110],  # mean distance = 395
            "PARTNER_C": [1000],  # distance = 500 (the max -> proximity 0.0)
        }
    )
    algorithm_results = [
        _frequency_result({"PARTNER_A": 10, "PARTNER_B": 5, "PARTNER_C": 1}),
        _domain_result("domain_retention", fisher_p=0.005, permutation_p=0.001),
        _domain_result("domain_disruption", fisher_p=0.02, permutation_p=0.01),
        _cutpoint_result(determinable=True, inferred_cutpoint_aa=500),
        _confidence_result(
            {
                "True": {"n": 20, "successes": 18, "ci_low": 0.65, "ci_high": 0.95},
                "False": {"n": 5, "successes": 1, "ci_low": 0.05, "ci_high": 0.60},
            }
        ),
    ]

    result = CompositeScoreAlgorithm().run(
        events, features, _FAKE_GENE, {"algorithm_results": algorithm_results}
    )

    assert isinstance(result, AlgorithmResult)
    assert result.Algorithm == "composite_score"
    rows = {row["Partner_gene"]: row for row in result.Tables["composite_evidence_ranking"]}
    assert set(rows) == {"PARTNER_A", "PARTNER_B", "PARTNER_C"}

    # domain_retention score = min(-log10(0.001), 10)/10 = 0.3
    # domain_disruption score = min(-log10(0.01), 10)/10 = 0.2
    # confidence certainty = 1 - mean(0.30, 0.55) = 0.575
    for row in rows.values():
        assert row["Domain_retention_score"] == pytest.approx(0.3)
        assert row["Domain_disruption_score"] == pytest.approx(0.2)
        assert row["Confidence_certainty_score"] == pytest.approx(0.575)
        assert set(row["Components_applicable"]) == {
            "recurrence",
            "domain_retention",
            "domain_disruption",
            "cutpoint_proximity",
            "confidence_certainty",
        }

    assert rows["PARTNER_A"]["Recurrence_score"] == pytest.approx(10 / 16)
    assert rows["PARTNER_A"]["Cutpoint_proximity_score"] == pytest.approx(1 - (2 / 3) / 500)
    assert rows["PARTNER_B"]["Cutpoint_proximity_score"] == pytest.approx(1 - 395 / 500)
    assert rows["PARTNER_C"]["Cutpoint_proximity_score"] == pytest.approx(0.0)

    def _expected_composite(row: dict) -> float:
        w = DEFAULT_WEIGHTS
        return (
            w["recurrence"] * row["Recurrence_score"]
            + w["domain_retention"] * row["Domain_retention_score"]
            + w["domain_disruption"] * row["Domain_disruption_score"]
            + w["cutpoint_proximity"] * row["Cutpoint_proximity_score"]
            + w["confidence_certainty"] * row["Confidence_certainty_score"]
        )

    for row in rows.values():
        assert row["Composite_score"] == pytest.approx(_expected_composite(row))

    # Ranked descending by composite score, with a dense 1..N Rank column.
    ranking = result.Tables["composite_evidence_ranking"]
    assert [row["Rank"] for row in ranking] == [1, 2, 3]
    assert [row["Composite_score"] for row in ranking] == sorted(
        (row["Composite_score"] for row in ranking), reverse=True
    )
    assert result.Warnings == []


def test_cutpoint_proximity_can_outrank_higher_recurrence_partner():
    """The composite score is a genuine aggregation, not just a recurrence
    sort: a low-recurrence partner whose breakpoints sit exactly on the
    inferred cutpoint can outrank a higher-recurrence partner whose
    breakpoints sit far from it, once the sub-scores are combined.
    """
    events, features = _events_and_features(
        {
            "HIGH_RECURRENCE_FAR": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # far from cutpoint
            "LOW_RECURRENCE_ON_CUTPOINT": [500],  # exactly at the inferred cutpoint
        }
    )
    algorithm_results = [
        _frequency_result({"HIGH_RECURRENCE_FAR": 10, "LOW_RECURRENCE_ON_CUTPOINT": 1}),
        _cutpoint_result(determinable=True, inferred_cutpoint_aa=500),
    ]

    result = CompositeScoreAlgorithm().run(
        events,
        features,
        _FAKE_GENE,
        {
            "algorithm_results": algorithm_results,
            "weights": {"recurrence": 0.2, "cutpoint_proximity": 0.8},
        },
    )

    ranking = result.Tables["composite_evidence_ranking"]
    assert ranking[0]["Partner_gene"] == "LOW_RECURRENCE_ON_CUTPOINT"
    assert ranking[0]["Rank"] == 1


def test_gene_with_only_domain_retention_configured_excludes_the_rest_gracefully():
    """Mirrors the real RET case: only frequency + domain_retention are
    applicable. The other sub-scores must be excluded from the weighted
    average entirely, not silently scored as zero evidence.
    """
    events, features = _events_and_features({"KIF5B": [200, 210], "CCDC6": [220]})
    algorithm_results = [
        _frequency_result({"KIF5B": 2, "CCDC6": 1}),
        _domain_result("domain_retention", fisher_p=0.02, permutation_p=0.01),
        # No domain_disruption, cutpoint_detection, or confidence_stats result at all.
    ]

    result = CompositeScoreAlgorithm().run(
        events, features, _FAKE_GENE_NO_DISRUPTION, {"algorithm_results": algorithm_results}
    )

    ranking = result.Tables["composite_evidence_ranking"]
    for row in ranking:
        assert set(row["Components_applicable"]) == {"recurrence", "domain_retention"}
        assert row["Domain_disruption_score"] is None
        assert row["Cutpoint_proximity_score"] is None
        assert row["Confidence_certainty_score"] is None

        w = DEFAULT_WEIGHTS
        weight_sum = w["recurrence"] + w["domain_retention"]
        expected = (
            w["recurrence"] * row["Recurrence_score"] + w["domain_retention"] * 0.2
        ) / weight_sum
        assert row["Composite_score"] == pytest.approx(expected)

    assert result.Summary["components_applicable"] == {
        "recurrence": True,
        "domain_retention": True,
        "domain_disruption": False,
        "cutpoint_proximity": False,
        "confidence_certainty": False,
    }
    warning_text = " ".join(result.Warnings)
    assert "domain_disruption" in warning_text
    assert "cutpoint_proximity" in warning_text
    assert "confidence_certainty" in warning_text


def test_domain_disruption_no_op_result_is_excluded_not_zeroed():
    """When domain_disruption ran but the gene has no
    disruption_required_domains configured, it returns a graceful no-op
    result (None statistics) -- composite_score must exclude it from the
    weighted average exactly like a missing result, never average in a 0.0.
    """
    events, features = _events_and_features({"PARTNER_A": [500]})
    no_op = _domain_result(
        "domain_disruption",
        fisher_p=None,
        permutation_p=None,
        warnings=["FAKE1 has no disruption_required_domains configured; skipped."],
    )
    with_disruption = CompositeScoreAlgorithm().run(
        events,
        features,
        _FAKE_GENE_NO_DISRUPTION,
        {"algorithm_results": [_frequency_result({"PARTNER_A": 1}), no_op]},
    )
    without_disruption_result_at_all = CompositeScoreAlgorithm().run(
        events,
        features,
        _FAKE_GENE_NO_DISRUPTION,
        {"algorithm_results": [_frequency_result({"PARTNER_A": 1})]},
    )

    row_with = with_disruption.Tables["composite_evidence_ranking"][0]
    row_without = without_disruption_result_at_all.Tables["composite_evidence_ranking"][0]
    assert row_with["Domain_disruption_score"] is None
    assert row_with["Composite_score"] == pytest.approx(row_without["Composite_score"])


def test_orchestrator_failed_result_is_treated_as_unavailable():
    """A result the orchestrator wrapped after an exception (Warnings
    starting with 'Algorithm failed') must be excluded like a missing one.
    """
    events, features = _events_and_features({"PARTNER_A": [500]})
    failed_confidence = AlgorithmResult(
        Algorithm="confidence_stats",
        Summary={},
        Tables={},
        Warnings=["Algorithm failed: ValueError: params['group_field'] is required"],
    )

    result = CompositeScoreAlgorithm().run(
        events,
        features,
        _FAKE_GENE_NO_DISRUPTION,
        {"algorithm_results": [_frequency_result({"PARTNER_A": 1}), failed_confidence]},
    )

    row = result.Tables["composite_evidence_ranking"][0]
    assert row["Confidence_certainty_score"] is None
    assert "confidence_certainty" not in row["Components_applicable"]


def test_accepts_json_round_tripped_algorithm_results():
    """The real-data path loads AlgorithmResult objects back out of a
    committed runs/*/results.json file, i.e. as plain dicts, not live
    AlgorithmResult instances. composite_score must accept both.
    """
    events, features = _events_and_features({"PARTNER_A": [500]})
    frequency_dict = _frequency_result({"PARTNER_A": 1}).model_dump(mode="json")
    domain_retention_dict = _domain_result(
        "domain_retention", fisher_p=0.03, permutation_p=0.02
    ).model_dump(mode="json")

    result = CompositeScoreAlgorithm().run(
        events,
        features,
        _FAKE_GENE_NO_DISRUPTION,
        {"algorithm_results": [frequency_dict, domain_retention_dict]},
    )

    row = result.Tables["composite_evidence_ranking"][0]
    assert row["Domain_retention_score"] == pytest.approx(_clipped_neg_log10(0.02, cap=10.0))


def test_result_schema_matches_canonical_algorithm_result_fields():
    events, features = _events_and_features({"PARTNER_A": [500]})
    result = CompositeScoreAlgorithm().run(
        events,
        features,
        _FAKE_GENE_NO_DISRUPTION,
        {"algorithm_results": [_frequency_result({"PARTNER_A": 1})]},
    )
    assert set(type(result).model_fields) == set(AlgorithmResult.model_fields)
