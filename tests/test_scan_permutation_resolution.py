"""Tests for the cohort scan's permutation budgets and floor annotation."""

from __future__ import annotations

from cfh.cohort.scan import (
    ADAPTIVE_ALGORITHMS,
    _adaptive_algorithm_params,
    _annotate_permutation_resolution,
)


def test_every_permutation_algorithm_gets_the_max_budget_even_without_adaptive():
    params = _adaptive_algorithm_params(False, 100, 50_000)

    assert set(params) == set(ADAPTIVE_ALGORITHMS)
    # The descriptive scans keep their fixed base budget.
    assert not {"cutpoint_detection", "window_detection"} & set(params)
    assert all(value == {"n_permutations_max": 50_000} for value in params.values())


def test_adaptive_params_carry_small_and_max_budgets():
    params = _adaptive_algorithm_params(True, 100, 50_000)

    assert params["junction_concentration"] == {
        "n_permutations_max": 50_000,
        "adaptive": True,
        "n_permutations_small": 100,
    }


def test_floor_results_that_miss_significance_are_resolution_limited():
    rows = [
        # At the floor of 1,000 permutations and not significant: limited.
        {
            "test": "permutation_corrected",
            "raw_p": 1 / 1001,
            "n_permutations": 1_000,
            "bh_adjusted_q": 0.2,
        },
        # At the floor but already significant: not limited.
        {
            "test": "permutation_corrected",
            "raw_p": 1 / 100_001,
            "n_permutations": 100_000,
            "bh_adjusted_q": 0.01,
        },
        # Off the floor: not limited.
        {"test": "permutation", "raw_p": 0.3, "n_permutations": 1_000, "bh_adjusted_q": 0.9},
        # Fisher p-values have no permutation floor.
        {"test": "fisher", "raw_p": 1 / 1001, "n_permutations": 1_000, "bh_adjusted_q": 0.2},
        # Algorithms without permutations.
        {"test": "enrichment", "raw_p": 0.001, "n_permutations": None, "bh_adjusted_q": 0.2},
    ]

    limited = _annotate_permutation_resolution(rows, 0.05)

    assert limited == 1
    assert [row["at_permutation_floor"] for row in rows] == [True, True, False, False, False]
    assert [row["resolution_limited"] for row in rows] == [True, False, False, False, False]
