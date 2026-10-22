"""Tests for adaptive permutation-count budgeting: the shared helper module
and its opt-in integration into domain_retention/domain_disruption/
cutpoint_detection, without changing any non-adaptive (default) behavior.
"""

from __future__ import annotations

import pytest

from cfh.algorithms.cutpoint_detection import CutpointDetectionAlgorithm
from cfh.algorithms.domain_retention import DomainRetentionAlgorithm
from cfh.genes.registry import load_gene_config
from cfh.model.fusion_event import FusionEvent
from cfh.model.fusion_feature import FusionFeature
from cfh.stats.adaptive_permutation import (
    DEFAULT_BORDERLINE_FACTOR,
    DEFAULT_SIGNIFICANCE_THRESHOLD,
    DEFAULT_SMALL_N_PERMUTATIONS,
    at_permutation_floor,
    is_borderline,
    resolve_permutation_budget,
    run_permutation_budget,
    should_escalate,
)


def test_is_borderline_true_near_threshold():
    assert is_borderline(0.05, threshold=0.05, factor=2.0) is True
    assert is_borderline(0.03, threshold=0.05, factor=2.0) is True
    assert is_borderline(0.09, threshold=0.05, factor=2.0) is True


def test_is_borderline_false_far_from_threshold():
    assert is_borderline(0.9, threshold=0.05, factor=2.0) is False
    assert is_borderline(0.0001, threshold=0.05, factor=2.0) is False


def test_is_borderline_none_p_value_never_borderline():
    assert is_borderline(None) is False


def test_is_borderline_rejects_non_positive_factor():
    with pytest.raises(ValueError):
        is_borderline(0.05, factor=0)


def test_resolve_permutation_budget_defaults_are_non_adaptive():
    budget = resolve_permutation_budget({}, default_full_n=10_000)
    assert budget["adaptive"] is False
    assert budget["full_n"] == 10_000
    assert budget["small_n"] == DEFAULT_SMALL_N_PERMUTATIONS
    assert budget["threshold"] == DEFAULT_SIGNIFICANCE_THRESHOLD
    assert budget["factor"] == DEFAULT_BORDERLINE_FACTOR


def test_resolve_permutation_budget_reads_existing_n_permutations_key():
    budget = resolve_permutation_budget({"n_permutations": 500}, default_full_n=10_000)
    assert budget["full_n"] == 500


def test_resolve_permutation_budget_honors_adaptive_overrides():
    budget = resolve_permutation_budget(
        {
            "adaptive": True,
            "n_permutations_small": 25,
            "n_permutations": 2_000,
            "significance_threshold": 0.01,
            "borderline_factor": 3.0,
        },
        default_full_n=10_000,
    )
    assert budget == {
        "adaptive": True,
        "small_n": 25,
        "full_n": 2_000,
        "threshold": 0.01,
        "factor": 3.0,
        "max_n": 2_000,
    }


# --- Integration: domain_retention -----------------------------------------


def _clearly_significant_events_and_features() -> tuple[list[FusionEvent], list[FusionFeature]]:
    """6 in-frame retained vs 8 out-of-frame lost -- the same clean-separation
    fixture already used in tests/benchmark/test_braf_kinase_retention_msk_impact_50k.py.
    """
    events, features = [], []
    for i in range(6):
        events.append(
            FusionEvent(
                Event_id=f"IN-{i}",
                Cohort="test",
                Frame_status="in-frame",
                Is_protein_fusion=True,
                Three_prime_gene="BRAF",
            )
        )
        features.append(
            FusionFeature(
                Event_id=f"IN-{i}",
                Gene="BRAF",
                Role="three_prime",
                Junction_position_aa=600 + i,
                Domain_retention_flags={"kinase": "retained"},
            )
        )
    for i in range(8):
        events.append(
            FusionEvent(
                Event_id=f"OUT-{i}",
                Cohort="test",
                Frame_status="out-of-frame",
                Is_protein_fusion=True,
                Three_prime_gene="BRAF",
            )
        )
        features.append(
            FusionFeature(
                Event_id=f"OUT-{i}",
                Gene="BRAF",
                Role="three_prime",
                Junction_position_aa=100 + i,
                Domain_retention_flags={"kinase": "lost"},
            )
        )
    return events, features


def test_domain_retention_non_adaptive_default_matches_existing_behavior():
    """Regression: omitting every adaptive key must reproduce the exact
    pre-adaptive-feature Parameters/Summary shape and value."""
    events, features = _clearly_significant_events_and_features()
    config = load_gene_config("braf")

    result = DomainRetentionAlgorithm().run(
        events, features, config, {"seed": 42, "n_permutations": 500}
    )

    assert result.Parameters == {"seed": 42, "n_permutations": 500, "adaptive": False}
    assert "adaptive_permutations" not in (result.Summary or {})


def test_domain_retention_adaptive_escalates_a_clearly_significant_result():
    """A clean 6/0 vs 0/8 separation sits at the small budget's floor, which is
    too coarse for a genome-wide family, so it must get the full budget."""
    events, features = _clearly_significant_events_and_features()
    config = load_gene_config("braf")

    result = DomainRetentionAlgorithm().run(
        events,
        features,
        config,
        {
            "seed": 42,
            "n_permutations": 5_000,
            "adaptive": True,
            "n_permutations_small": 50,
        },
    )

    adaptive_info = result.Summary["adaptive_permutations"]
    assert adaptive_info["enabled"] is True
    assert adaptive_info["n_permutations_small"] == 50
    assert adaptive_info["n_permutations_full"] == 5_000
    assert adaptive_info["escalated_to_full"] is True
    assert adaptive_info["n_permutations_used"] == 5_000
    assert result.Parameters["n_permutations"] == 5_000


def test_domain_retention_adaptive_keeps_small_budget_when_clearly_null(monkeypatch):
    import cfh.algorithms.domain_retention as module

    events, features = _clearly_significant_events_and_features()
    config = load_gene_config("braf")
    calls: list[int] = []

    def _fake_permutation_null_test(*args, **kwargs):
        calls.append(kwargs["n_permutations"])
        return 0.9, 0.5, tuple()

    monkeypatch.setattr(module, "permutation_null_test", _fake_permutation_null_test)

    result = DomainRetentionAlgorithm().run(
        events,
        features,
        config,
        {"seed": 1, "n_permutations": 3_000, "adaptive": True, "n_permutations_small": 100},
    )

    assert calls == [100]
    assert result.Summary["adaptive_permutations"]["escalated_to_full"] is False


def test_domain_retention_adaptive_is_deterministic_for_a_fixed_seed():
    events, features = _clearly_significant_events_and_features()
    config = load_gene_config("braf")
    params = {"seed": 7, "n_permutations": 1_000, "adaptive": True, "n_permutations_small": 50}

    first = DomainRetentionAlgorithm().run(events, features, config, dict(params))
    second = DomainRetentionAlgorithm().run(events, features, config, dict(params))

    assert first.Summary == second.Summary


def test_domain_retention_adaptive_escalates_when_borderline(monkeypatch):
    """Force a borderline small-N p-value and confirm escalation actually
    reruns with the full permutation count (a different, larger n)."""
    import cfh.algorithms.domain_retention as module

    events, features = _clearly_significant_events_and_features()
    config = load_gene_config("braf")

    calls: list[int] = []
    original = module.permutation_null_test

    def _fake_permutation_null_test(*args, **kwargs):
        n_permutations = kwargs["n_permutations"]
        calls.append(n_permutations)
        table = module.build_frame_domain_contingency_table(events, features, config)
        _, _ = module.fishers_frame_domain_test(table)
        # First (small) call reports a borderline p-value; the escalated
        # (full) call reports a clearly non-borderline one.
        p_value = 0.05 if len(calls) == 1 else 0.9
        return p_value, 0.5, tuple()

    monkeypatch.setattr(module, "permutation_null_test", _fake_permutation_null_test)

    result = DomainRetentionAlgorithm().run(
        events,
        features,
        config,
        {"seed": 1, "n_permutations": 3_000, "adaptive": True, "n_permutations_small": 100},
    )

    assert calls == [100, 3_000]
    assert result.Summary["adaptive_permutations"]["escalated_to_full"] is True
    assert result.Summary["permutation_empirical_p_value"] == 0.9
    monkeypatch.setattr(module, "permutation_null_test", original)


# --- Integration: cutpoint_detection ----------------------------------------


def _cutpoint_events_and_features() -> tuple[list[FusionEvent], list[FusionFeature]]:
    events, features = [], []
    positions_and_statuses = [
        (100, "retained"),
        (110, "retained"),
        (120, "retained"),
        (600, "lost"),
        (610, "lost"),
        (620, "lost"),
    ]
    for i, (position, status) in enumerate(positions_and_statuses):
        events.append(FusionEvent(Event_id=f"E-{i}", Cohort="test"))
        features.append(
            FusionFeature(
                Event_id=f"E-{i}",
                Gene="BRAF",
                Junction_position_aa=position,
                Domain_retention_flags={"kinase": status},
            )
        )
    return events, features


def test_cutpoint_detection_non_adaptive_default_matches_existing_behavior():
    events, features = _cutpoint_events_and_features()
    config = load_gene_config("braf")

    result = CutpointDetectionAlgorithm().run(
        events, features, config, {"seed": 42, "n_permutations": 200}
    )

    assert result.Parameters == {"seed": 42, "n_permutations": 200, "adaptive": False}
    assert "adaptive_permutations" not in (result.Summary or {})


def test_cutpoint_detection_adaptive_skips_escalation_when_indeterminable():
    """Too few events to even scan: adaptive mode must not attempt to
    escalate (there is no p-value to judge as borderline)."""
    config = load_gene_config("braf")
    events = [FusionEvent(Event_id="E-0", Cohort="test")]
    features = [
        FusionFeature(
            Event_id="E-0",
            Gene="BRAF",
            Junction_position_aa=100,
            Domain_retention_flags={"kinase": "retained"},
        )
    ]

    result = CutpointDetectionAlgorithm().run(
        events,
        features,
        config,
        {"seed": 42, "n_permutations": 500, "adaptive": True, "n_permutations_small": 10},
    )

    assert result.Summary["determinable"] is False
    assert result.Summary["adaptive_permutations"]["escalated_to_full"] is False
    assert result.Summary["adaptive_permutations"]["n_permutations_used"] == 10


def test_cutpoint_detection_adaptive_is_deterministic_for_a_fixed_seed():
    events, features = _cutpoint_events_and_features()
    config = load_gene_config("braf")
    params = {"seed": 3, "n_permutations": 500, "adaptive": True, "n_permutations_small": 20}

    first = CutpointDetectionAlgorithm().run(events, features, config, dict(params))
    second = CutpointDetectionAlgorithm().run(events, features, config, dict(params))

    assert first.Summary == second.Summary


# --- Escalation rule and floor-driven escalation -----------------------------


def test_should_escalate_includes_results_far_below_threshold():
    assert should_escalate(0.0001, threshold=0.05, factor=2.0) is True
    assert should_escalate(0.05, threshold=0.05, factor=2.0) is True
    assert should_escalate(0.1, threshold=0.05, factor=2.0) is True
    assert should_escalate(0.11, threshold=0.05, factor=2.0) is False
    assert should_escalate(None) is False


def test_at_permutation_floor():
    assert at_permutation_floor(1 / 101, 100) is True
    assert at_permutation_floor(2 / 101, 100) is False
    assert at_permutation_floor(None, 100) is False


def _floor_until(n_with_exceedance: int):
    """A fake test whose p-value sits at its floor below ``n_with_exceedance``."""
    calls: list[int] = []

    def run(n: int) -> float:
        calls.append(n)
        return 1 / (n + 1) if n < n_with_exceedance else 5 / (n + 1)

    return run, calls


def test_run_permutation_budget_escalates_floor_results_up_to_max():
    run, calls = _floor_until(10**9)
    budget = resolve_permutation_budget(
        {"adaptive": True, "n_permutations": 1_000, "n_permutations_max": 50_000},
        default_full_n=1_000,
    )

    p_value, used, escalated = run_permutation_budget(run, lambda p: p, budget)

    assert calls == [100, 1_000, 10_000, 50_000]
    assert used == 50_000
    assert escalated is True
    assert p_value == 1 / 50_001


def test_run_permutation_budget_stops_once_off_the_floor():
    run, calls = _floor_until(5_000)
    budget = resolve_permutation_budget(
        {"n_permutations": 1_000, "n_permutations_max": 1_000_000}, default_full_n=1_000
    )

    _, used, _ = run_permutation_budget(run, lambda p: p, budget)

    assert calls == [1_000, 10_000]
    assert used == 10_000


def test_run_permutation_budget_without_max_never_exceeds_full_budget():
    run, calls = _floor_until(10**9)
    budget = resolve_permutation_budget({"n_permutations": 1_000}, default_full_n=1_000)

    run_permutation_budget(run, lambda p: p, budget)

    assert calls == [1_000]
