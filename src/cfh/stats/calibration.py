"""Deterministic offline stress tests for breakpoint label-separation scans.

These simulations exercise the cutpoint/window hypothesis tests. They do not
test genomic breakpoint density or establish real
cohort replication. Repeated observations are deliberately kept as separate
rows in one scenario to expose possible anti-conservative behavior.
"""

from __future__ import annotations

import argparse
import json
import math
from collections.abc import Sequence
from typing import Any

import numpy as np

from cfh.stats.adaptive_permutation import permutation_resolution
from cfh.stats.cutpoint_scan import detect_cutpoint
from cfh.stats.window_scan import detect_window


def wilson_interval(successes: int, trials: int, z: float = 1.959963984540054) -> list[float]:
    """Two-sided approximate 95% binomial interval by default."""
    if trials <= 0 or not 0 <= successes <= trials:
        raise ValueError("Require 0 <= successes <= trials and trials > 0")
    p = successes / trials
    denominator = 1 + z * z / trials
    center = (p + z * z / (2 * trials)) / denominator
    half = z * math.sqrt(p * (1 - p) / trials + z * z / (4 * trials * trials))
    return [max(0.0, center - half / denominator), min(1.0, center + half / denominator)]


def simulated_observations(
    scenario: str, rng: np.random.Generator, *, n_patients: int = 24
) -> tuple[list[int], list[str], list[str]]:
    """Generate mapped positions, binary labels and stable patient IDs.

    In ``planted_window`` positive-label probability is 0.85 inside the
    100-aa interval and 0.15 outside. ``mapped_pile_null`` maps half of
    independent genomic positions to one intronic boundary. ``repeated_null``
    replicates each patient's position and label three times, violating the
    event-level exchangeability assumption of the current scans.
    """
    if n_patients < 8:
        raise ValueError("n_patients must be at least 8")
    if scenario not in {
        "iid_null",
        "uneven_null",
        "planted_window",
        "planted_cutpoint",
        "mapped_pile_null",
        "repeated_null",
    }:
        raise ValueError(f"Unknown scenario: {scenario}")
    raw_positions = rng.integers(1, 1001, size=n_patients)
    if scenario == "uneven_null":
        raw_positions = 1 + np.floor(999 * rng.beta(2, 5, size=n_patients)).astype(int)
    if scenario == "planted_window":
        inside = (raw_positions >= 400) & (raw_positions <= 500)
        probabilities = np.where(inside, 0.85, 0.15)
        labels = rng.random(n_patients) < probabilities
    elif scenario == "planted_cutpoint":
        labels = rng.random(n_patients) < np.where(raw_positions <= 500, 0.9, 0.1)
    else:
        labels = rng.random(n_patients) < 0.5
    if scenario == "mapped_pile_null":
        # A simplified exon-boundary snapping artifact, independent of label.
        raw_positions = np.where(raw_positions < 500, 500, raw_positions)
    repetitions = 3 if scenario == "repeated_null" else 1
    positions = np.repeat(raw_positions, repetitions).astype(int).tolist()
    statuses = ["retained" if label else "lost" for label in np.repeat(labels, repetitions)]
    patient_ids = [f"P{i}" for i in np.repeat(np.arange(n_patients), repetitions)]
    return positions, statuses, patient_ids


def calibrate_scans(
    *,
    replicates: int = 20,
    n_patients: int = 24,
    n_permutations: int = 99,
    seed: int = 42,
    alpha: float = 0.05,
    family_size: int = 1,
    widths: Sequence[int] = (100,),
) -> dict[str, Any]:
    """Return empirical rejection rates, Wilson intervals and diagnostics.

    A null rejection rate estimates the false-positive rate for that
    simulation's data-generating process; the planted rate estimates power
    for its specified effect. All replicates use independent seeded streams.
    """
    if replicates <= 0 or not 0 < alpha < 1:
        raise ValueError("replicates must be positive and alpha between 0 and 1")
    if not widths or any(width <= 0 for width in widths):
        raise ValueError("widths must contain positive values")
    resolution = permutation_resolution(n_permutations, family_size=family_size)
    rng = np.random.default_rng(seed)
    scenarios: dict[str, Any] = {}
    for scenario in (
        "iid_null",
        "uneven_null",
        "planted_cutpoint",
        "planted_window",
        "mapped_pile_null",
        "repeated_null",
    ):
        outcomes: dict[str, list[bool]] = {"cutpoint": [], "window": []}
        determinable = {"cutpoint": 0, "window": 0}
        for _ in range(replicates):
            positions, statuses, patient_ids = simulated_observations(
                scenario, rng, n_patients=n_patients
            )
            run_seed = int(rng.integers(0, 2**32))
            event_ids = [f"E{i}" for i in range(len(positions))]
            results = {
                "cutpoint": detect_cutpoint(
                    positions, statuses, seed=run_seed, n_permutations=n_permutations
                ),
                "window": detect_window(
                    positions,
                    statuses,
                    event_ids,
                    seed=run_seed,
                    n_permutations=n_permutations,
                    widths=widths,
                ),
            }
            for method, result in results.items():
                determinable[method] += bool(result["determinable"])
                p_value = result["corrected_p_value"]
                outcomes[method].append(p_value is not None and p_value < alpha)
            assert len(patient_ids) == len(positions)
        scenarios[scenario] = {
            method: {
                "rejections": sum(calls),
                "replicates": replicates,
                "determinable": determinable[method],
                "rate": sum(calls) / replicates,
                "wilson_95_interval": wilson_interval(sum(calls), replicates),
            }
            for method, calls in outcomes.items()
        }
    return {
        "scope": "cutpoint/window label-separation tests; no breakpoint-density test exists here",
        "seed": seed,
        "alpha": alpha,
        "n_patients": n_patients,
        "n_permutations": n_permutations,
        "widths": list(widths),
        "permutation_resolution": resolution,
        "scenarios": scenarios,
        "notes": [
            "Planted-window rejection rate is sensitivity only for the specified synthetic effect.",
            "Repeated-patient scenario retains duplicate events intentionally; rates do not "
            "represent patient-level permutation validity.",
            "Intervals describe simulation uncertainty, not uncertainty in a real cohort.",
            "Across-gene FDR is not established by these scan-level simulations.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=20)
    parser.add_argument("--n-patients", type=int, default=24)
    parser.add_argument("--n-permutations", type=int, default=99)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--family-size", type=int, default=1)
    args = parser.parse_args()
    report = calibrate_scans(
        replicates=args.replicates,
        n_patients=args.n_patients,
        n_permutations=args.n_permutations,
        seed=args.seed,
        family_size=args.family_size,
    )
    print(json.dumps(report, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
