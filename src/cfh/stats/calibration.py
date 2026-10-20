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
from cfh.stats.window_scan import DEFAULT_WIDTHS, detect_window

DEFAULT_POWER_SIZES: tuple[int, ...] = (10, 25, 50, 100, 250)
DEFAULT_POWER_EFFECTS: tuple[float, ...] = (0.4, 0.6, 0.8)
POWER_WINDOW_AA: tuple[int, int] = (400, 500)
POWER_CUTPOINT_AA = 500


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


def collapse_patient_observations(
    positions: Sequence[int], statuses: Sequence[str], patient_ids: Sequence[str]
) -> tuple[list[int], list[str], list[str]]:
    """Keep one observation per (patient, position), as the live pipeline does.

    Mirrors :func:`cfh.stats.observation_units.collapse_repeated_observations`
    for simulated data, where every observation shares one partner and role.
    """
    seen: set[tuple[str, int]] = set()
    kept_positions: list[int] = []
    kept_statuses: list[str] = []
    kept_patients: list[str] = []
    for position, status, patient in zip(positions, statuses, patient_ids, strict=True):
        if (patient, position) in seen:
            continue
        seen.add((patient, position))
        kept_positions.append(position)
        kept_statuses.append(status)
        kept_patients.append(patient)
    return kept_positions, kept_statuses, kept_patients


def calibrate_scans(
    *,
    replicates: int = 20,
    n_patients: int = 24,
    n_permutations: int = 99,
    seed: int = 42,
    alpha: float = 0.05,
    family_size: int = 1,
    widths: Sequence[int] = (100,),
    observation_unit: str = "event",
) -> dict[str, Any]:
    """Return empirical rejection rates, Wilson intervals and diagnostics.

    ``observation_unit="patient"`` collapses repeated observations of one
    patient before scanning, matching the live pipeline's default.

    A null rejection rate estimates the false-positive rate for that
    simulation's data-generating process; the planted rate estimates power
    for its specified effect. All replicates use independent seeded streams.
    """
    if replicates <= 0 or not 0 < alpha < 1:
        raise ValueError("replicates must be positive and alpha between 0 and 1")
    if not widths or any(width <= 0 for width in widths):
        raise ValueError("widths must contain positive values")
    if observation_unit not in {"event", "patient"}:
        raise ValueError("observation_unit must be 'event' or 'patient'")
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
            if observation_unit == "patient":
                positions, statuses, patient_ids = collapse_patient_observations(
                    positions, statuses, patient_ids
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
        "observation_unit": observation_unit,
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


def planted_observations(
    kind: str, n_events: int, effect: float, rng: np.random.Generator
) -> tuple[list[int], list[str]]:
    """Independent events with a planted retained-probability difference.

    Positions are uniform on 1..1000 aa. Inside the signal region (positions
    <= 500 for ``cutpoint``; 400..500 for ``window``) an event is retained
    with probability ``0.5 + effect / 2``; outside, ``0.5 - effect / 2``.
    ``effect`` is therefore the absolute difference in retained probability.
    """
    if kind not in {"cutpoint", "window"}:
        raise ValueError(f"Unknown planted signal: {kind}")
    if not 0 <= effect <= 1:
        raise ValueError("effect must be between 0 and 1")
    if n_events < 1:
        raise ValueError("n_events must be positive")
    positions = rng.integers(1, 1001, size=n_events)
    if kind == "cutpoint":
        inside = positions <= POWER_CUTPOINT_AA
    else:
        inside = (positions >= POWER_WINDOW_AA[0]) & (positions <= POWER_WINDOW_AA[1])
    labels = rng.random(n_events) < np.where(inside, 0.5 + effect / 2, 0.5 - effect / 2)
    return positions.astype(int).tolist(), ["retained" if label else "lost" for label in labels]


def power_grid(
    *,
    sizes: Sequence[int] = DEFAULT_POWER_SIZES,
    effects: Sequence[float] = DEFAULT_POWER_EFFECTS,
    replicates: int = 50,
    n_permutations: int = 99,
    seed: int = 42,
    alpha: float = 0.05,
    widths: Sequence[int] = DEFAULT_WIDTHS,
    signals: Sequence[str] = ("cutpoint", "window"),
) -> dict[str, Any]:
    """Rejection rate of each scan against its own planted signal.

    Each (signal, events, effect) cell runs ``replicates`` independent
    simulations. The cutpoint scan is tested on a planted cutpoint and the
    window scan on a planted 100-aa internal window, using the production
    window widths by default. ``n_permutations`` matches the cohort scan's
    first adaptive stage (100 draws); borderline results there escalate to a
    larger budget, which this grid does not model.
    """
    if replicates <= 0 or not 0 < alpha < 1:
        raise ValueError("replicates must be positive and alpha between 0 and 1")
    if not sizes or not effects:
        raise ValueError("sizes and effects must be nonempty")
    if not signals or set(signals) - {"cutpoint", "window"}:
        raise ValueError("signals must be a nonempty subset of {'cutpoint', 'window'}")
    rng = np.random.default_rng(seed)
    rows: list[dict[str, Any]] = []
    for kind in signals:
        for n_events in sizes:
            for effect in effects:
                rejections = 0
                determinable = 0
                for _ in range(replicates):
                    positions, statuses = planted_observations(kind, n_events, effect, rng)
                    run_seed = int(rng.integers(0, 2**32))
                    if kind == "cutpoint":
                        result = detect_cutpoint(
                            positions, statuses, seed=run_seed, n_permutations=n_permutations
                        )
                    else:
                        result = detect_window(
                            positions,
                            statuses,
                            [f"E{i}" for i in range(n_events)],
                            seed=run_seed,
                            n_permutations=n_permutations,
                            widths=widths,
                        )
                    determinable += bool(result["determinable"])
                    p_value = result["corrected_p_value"]
                    rejections += p_value is not None and p_value < alpha
                rows.append(
                    {
                        "signal": kind,
                        "n_events": n_events,
                        "effect": effect,
                        "rejections": rejections,
                        "replicates": replicates,
                        "determinable": determinable,
                        "power": rejections / replicates,
                        "wilson_95_interval": wilson_interval(rejections, replicates),
                    }
                )
    return {
        "scope": "power of cutpoint/window label-separation scans for planted signals",
        "seed": seed,
        "alpha": alpha,
        "n_permutations": n_permutations,
        "widths": list(widths),
        "cutpoint_aa": POWER_CUTPOINT_AA,
        "window_aa": list(POWER_WINDOW_AA),
        "effect_definition": "absolute difference in retained probability inside vs outside",
        "rows": rows,
        "notes": [
            "Events are independent; repeated-patient dependence is not modeled.",
            "Positions are uniform on 1..1000 aa; real genes have clustered breakpoints.",
            "Adaptive escalation beyond the first permutation stage is not modeled.",
        ],
    }


def power_markdown(report: dict[str, Any]) -> str:
    """Render a power grid as one Markdown table per signal."""
    effects = sorted({row["effect"] for row in report["rows"]})
    lines = [
        f"Power at alpha={report['alpha']} with {report['n_permutations']} permutations, "
        f"seed {report['seed']}. Cells: rejections/replicates (Wilson 95% interval).",
        "",
    ]
    for kind in ("cutpoint", "window"):
        rows = [row for row in report["rows"] if row["signal"] == kind]
        if not rows:
            continue
        lines += [
            f"**Planted {kind}**",
            "",
            "| Events | " + " | ".join(f"effect {effect:g}" for effect in effects) + " |",
            "| ---: |" + " ---: |" * len(effects),
        ]
        for n_events in sorted({row["n_events"] for row in rows}):
            cells = []
            for effect in effects:
                row = next(r for r in rows if r["n_events"] == n_events and r["effect"] == effect)
                low, high = row["wilson_95_interval"]
                cells.append(
                    f"{row['rejections']}/{row['replicates']} ({100 * low:.0f}–{100 * high:.0f}%)"
                )
            lines.append(f"| {n_events} | " + " | ".join(cells) + " |")
        lines.append("")
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=20)
    parser.add_argument("--n-patients", type=int, default=24)
    parser.add_argument("--n-permutations", type=int, default=99)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--family-size", type=int, default=1)
    parser.add_argument(
        "--observation-unit",
        choices=("event", "patient"),
        default="event",
        help="Collapse repeated patient observations before scanning (null scenarios only).",
    )
    parser.add_argument(
        "--power-grid",
        action="store_true",
        help="Estimate power over --sizes x --effects instead of running the null scenarios.",
    )
    parser.add_argument("--sizes", type=int, nargs="+", default=list(DEFAULT_POWER_SIZES))
    parser.add_argument("--effects", type=float, nargs="+", default=list(DEFAULT_POWER_EFFECTS))
    parser.add_argument(
        "--markdown", action="store_true", help="With --power-grid, print a Markdown table."
    )
    args = parser.parse_args()
    if args.power_grid:
        report = power_grid(
            sizes=args.sizes,
            effects=args.effects,
            replicates=args.replicates,
            n_permutations=args.n_permutations,
            seed=args.seed,
        )
        print(power_markdown(report) if args.markdown else json.dumps(report, indent=2))
        return
    report = calibrate_scans(
        replicates=args.replicates,
        n_patients=args.n_patients,
        n_permutations=args.n_permutations,
        seed=args.seed,
        family_size=args.family_size,
        observation_unit=args.observation_unit,
    )
    print(json.dumps(report, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
