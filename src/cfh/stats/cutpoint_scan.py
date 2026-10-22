"""Gene-agnostic, frequency/recurrence-only cutpoint (boundary) detection.

Given a gene's observed breakpoint protein positions and a binary outcome
label per breakpoint (e.g. domain-retention status), this scans every
candidate cutpoint along the protein and finds the position that best
separates the two outcome classes, using the same underlying Fisher's-exact
primitive (``scipy.stats.fisher_exact``) the domain-retention algorithm's
:func:`cfh.stats.breakpoint_tests.fishers_frame_domain_test` already uses --
no new statistics dependency is introduced. Unlike that fixed one-sided
test, a cutpoint scan makes no a-priori assumption about which side of a
candidate boundary the positive class should cluster on, so a two-sided
alternative is used here.

Scanning many candidate cutpoints and taking the best one inflates
significance (the multiple-comparisons problem), so the corrected p-value
is estimated by permuting the outcome labels across the observed breakpoint
positions and recomputing the max statistic under each permutation -- the
same permutation-test pattern already used by
:func:`cfh.stats.breakpoint_tests.permutation_null_test`.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

import numpy as np
from scipy.stats import fisher_exact

from cfh.stats.max_statistic_permutation import empirical_p_value, null_max_statistics

_POSITIVE_STATUS_DEFAULT = "retained"


def candidate_cutpoints(positions: Sequence[int]) -> list[int]:
    """Distinct observed positions except the largest, so every candidate
    cutpoint produces a non-empty "above" group.
    """
    distinct = sorted({int(position) for position in positions})
    return distinct[:-1]


def _contingency_table(
    positions: Sequence[int],
    statuses: Sequence[str],
    cutpoint: int,
    positive_status: str,
) -> list[list[int]]:
    table = [[0, 0], [0, 0]]
    for position, status in zip(positions, statuses):
        row = 0 if status == positive_status else 1
        column = 0 if position <= cutpoint else 1
        table[row][column] += 1
    return table


def _neg_log10_p(p_value: float) -> float:
    if p_value <= 0:
        return math.inf
    return -math.log10(p_value)


def fishers_cutpoint_test(contingency_table: Sequence[Sequence[int]]) -> tuple[float, float]:
    """Two-sided Fisher's exact test for one candidate cutpoint's 2x2 split."""
    table = np.asarray(contingency_table, dtype=int)
    if table.shape != (2, 2):
        raise ValueError("contingency_table must be a 2x2 table")
    if np.any(table < 0):
        raise ValueError("contingency_table counts must be non-negative")

    odds_ratio, p_value = fisher_exact(table, alternative="two-sided")
    return float(odds_ratio), float(p_value)


def scan_cutpoints(
    positions: Sequence[int],
    statuses: Sequence[str],
    positive_status: str = _POSITIVE_STATUS_DEFAULT,
) -> list[dict]:
    """Compute the separation statistic at every candidate cutpoint.

    Returns one row per candidate cutpoint (empty if fewer than two distinct
    positions are observed), suitable for a future lollipop-viz overlay.
    """
    if len(positions) != len(statuses):
        raise ValueError("positions and statuses must be the same length")

    scan: list[dict] = []
    for cutpoint in candidate_cutpoints(positions):
        table = _contingency_table(positions, statuses, cutpoint, positive_status)
        odds_ratio, p_value = fishers_cutpoint_test(table)
        scan.append(
            {
                "cutpoint": cutpoint,
                "n_positive_at_or_below": table[0][0],
                "n_positive_above": table[0][1],
                "n_negative_at_or_below": table[1][0],
                "n_negative_above": table[1][1],
                "odds_ratio": odds_ratio,
                "p_value": p_value,
                "neg_log10_p_value": _neg_log10_p(p_value),
            }
        )
    return scan


def _best_scan_row(scan: list[dict]) -> dict:
    """Pick the max-separation row, preferring the smallest cutpoint on ties."""
    return max(scan, key=lambda row: (row["neg_log10_p_value"], -row["cutpoint"]))


def detect_cutpoint(
    positions: Sequence[int],
    statuses: Sequence[str],
    *,
    positive_status: str = _POSITIVE_STATUS_DEFAULT,
    seed: int = 42,
    n_permutations: int = 10_000,
) -> dict:
    """Find the recurrence-based cutpoint that best separates ``statuses``.

    Degenerate inputs (too few events, a single distinct breakpoint
    position, or a single outcome class) never raise: they return
    ``determinable: False`` with a human-readable ``reason`` instead, so
    callers can surface a clear "not determinable" result.

    Returns a plain ``dict`` (JSON-serializable, suitable for an
    :class:`~cfh.model.algorithm_result.AlgorithmResult`) with keys:
    ``determinable``, ``reason``, ``best_cutpoint``, ``observed_statistic``
    (``-log10(p)`` at the best cutpoint), ``observed_p_value``,
    ``observed_odds_ratio``, ``corrected_p_value`` (the permutation-based
    empirical p-value), ``scan`` (the full per-cutpoint scan), and
    ``null_max_statistics`` (present only when determinable).
    """
    if n_permutations <= 0:
        raise ValueError("n_permutations must be positive")
    if len(positions) != len(statuses):
        raise ValueError("positions and statuses must be the same length")

    positions = [int(position) for position in positions]
    statuses = list(statuses)

    reason: str | None = None
    if len(positions) < 4:
        reason = "fewer than 4 mapped breakpoint events with a known outcome status"
    elif len(set(positions)) < 2:
        reason = "fewer than 2 distinct breakpoint positions to scan"
    elif len(set(statuses)) < 2:
        reason = "all events share a single outcome class; no separation is possible"

    scan = scan_cutpoints(positions, statuses, positive_status) if reason is None else []

    if reason is not None:
        return {
            "determinable": False,
            "reason": reason,
            "best_cutpoint": None,
            "observed_statistic": None,
            "observed_p_value": None,
            "observed_odds_ratio": None,
            "corrected_p_value": None,
            "scan": scan,
        }

    best = _best_scan_row(scan)
    observed_statistic = best["neg_log10_p_value"]

    position_array = np.asarray(positions)
    membership = np.stack(
        [position_array <= cutpoint for cutpoint in candidate_cutpoints(positions)]
    )
    null_statistics = null_max_statistics(
        membership,
        [status == positive_status for status in statuses],
        seed=seed,
        n_permutations=n_permutations,
    )
    corrected_p_value = empirical_p_value(observed_statistic, null_statistics)

    return {
        "determinable": True,
        "reason": None,
        "best_cutpoint": best["cutpoint"],
        "observed_statistic": observed_statistic,
        "observed_p_value": best["p_value"],
        "observed_odds_ratio": best["odds_ratio"],
        "corrected_p_value": float(corrected_p_value),
        "scan": scan,
        "null_max_statistics": tuple(float(stat) for stat in null_statistics),
    }
