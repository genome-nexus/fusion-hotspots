"""Small deterministic checks for the offline scan calibration runner."""

import numpy as np
import pytest

from cfh.stats.adaptive_permutation import permutation_resolution
from cfh.stats.calibration import calibrate_scans, simulated_observations, wilson_interval


def test_resolution_reports_scan_floor_and_across_family_diagnostic():
    report = permutation_resolution(99, family_size=100)
    assert report["minimum_empirical_p_value"] == pytest.approx(0.01)
    assert report["minimum_bonferroni_adjusted_p_value"] == 1.0
    assert report["permutations_for_strict_0_05_first_rank"] == 2000
    with pytest.raises(ValueError):
        permutation_resolution(0)


def test_simulations_preserve_patient_cluster_and_snapping():
    positions, statuses, patients = simulated_observations(
        "repeated_null", np.random.default_rng(7), n_patients=8
    )
    assert len(positions) == len(statuses) == len(patients) == 24
    for i in range(0, 24, 3):
        assert len(set(positions[i : i + 3])) == 1
        assert len(set(statuses[i : i + 3])) == 1
        assert len(set(patients[i : i + 3])) == 1
    positions, _, _ = simulated_observations(
        "mapped_pile_null", np.random.default_rng(7), n_patients=24
    )
    assert positions.count(500) > 1


def test_wilson_and_calibration_are_deterministic():
    assert wilson_interval(0, 20)[0] == 0
    assert wilson_interval(20, 20)[1] == 1
    kwargs = dict(replicates=2, n_patients=12, n_permutations=19, seed=13)
    first = calibrate_scans(**kwargs)
    assert first == calibrate_scans(**kwargs)
    assert first["permutation_resolution"]["minimum_empirical_p_value"] == 0.05
    assert set(first["scenarios"]) == {
        "iid_null",
        "uneven_null",
        "planted_cutpoint",
        "planted_window",
        "mapped_pile_null",
        "repeated_null",
    }


def test_family_budget_requires_fixed_attainable_resolution():
    from cfh.stats.adaptive_permutation import resolve_permutation_budget

    with pytest.raises(ValueError, match="cannot resolve"):
        resolve_permutation_budget(
            {"n_permutations": 100, "correction_family_size": 100}, default_full_n=100
        )
    with pytest.raises(ValueError, match="non-adaptive"):
        resolve_permutation_budget(
            {"n_permutations": 2000, "correction_family_size": 100, "adaptive": True},
            default_full_n=100,
        )
    with pytest.raises(ValueError, match="positive integer"):
        resolve_permutation_budget(
            {"n_permutations": 2000, "correction_family_size": 1.5}, default_full_n=100
        )
    result = resolve_permutation_budget(
        {"n_permutations": 2000, "correction_family_size": 100}, default_full_n=100
    )
    assert result["permutation_resolution"]["minimum_empirical_p_value"] < 0.0005
