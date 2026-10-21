"""Vectorized max-statistic label permutation for the cutpoint and window scans.

Both scans ask the same question of many candidate splits of a gene's
breakpoints: does a two-sided Fisher's exact test separate retained from
lost events? Candidate membership depends only on breakpoint positions, so
it is fixed across permutations; only the outcome labels move. Every
candidate's 2x2 table is therefore determined by two numbers, the events
inside it and the positives among them, with the cohort totals fixed. This
module evaluates each distinct ``(inside, positives)`` pair once with
``scipy.stats.fisher_exact`` and reuses the value, so the permutation null is
identical to calling the test per candidate per replicate, only faster.

Permutations are drawn one index vector at a time from
``np.random.default_rng(seed)``. That consumes the same random stream as
permuting the label array itself, so a given ``(seed, n_permutations)``
reproduces the scans' earlier per-replicate results exactly.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

import numpy as np
from scipy.stats import fisher_exact

_CHUNK_SIZE = 2_000


class FisherNegLog10Lookup:
    """``-log10`` two-sided Fisher p-values for one cohort's fixed totals."""

    def __init__(self, n_events: int, n_positive: int) -> None:
        self._n_events = n_events
        self._n_positive = n_positive
        self._values = np.full((n_events + 1, n_positive + 1), np.nan)

    def _table(self, inside: int, positive_inside: int) -> list[list[int]]:
        positive_outside = self._n_positive - positive_inside
        outside = self._n_events - inside
        return [
            [positive_inside, positive_outside],
            [inside - positive_inside, outside - positive_outside],
        ]

    def values(self, inside: np.ndarray, positive_inside: np.ndarray) -> np.ndarray:
        """Return ``-log10(p)`` for every ``(inside, positive_inside)`` pair."""
        found = self._values[inside, positive_inside]
        missing = np.isnan(found)
        if missing.any():
            pairs = np.unique(
                np.stack([inside[missing], positive_inside[missing]], axis=-1), axis=0
            )
            for n_inside, n_positive in pairs:
                _, p_value = fisher_exact(self._table(int(n_inside), int(n_positive)))
                self._values[n_inside, n_positive] = (
                    math.inf if p_value <= 0 else -math.log10(float(p_value))
                )
            found = self._values[inside, positive_inside]
        return found


def null_max_statistics(
    membership: np.ndarray,
    is_positive: Sequence[bool],
    *,
    seed: int,
    n_permutations: int,
) -> np.ndarray:
    """Max ``-log10`` Fisher p over candidates for each label permutation.

    ``membership`` is a ``(n_candidates, n_events)`` boolean matrix: row ``c``
    marks the events inside candidate ``c``.
    """
    if n_permutations <= 0:
        raise ValueError("n_permutations must be positive")
    labels = np.asarray(is_positive, dtype=bool)
    membership = np.asarray(membership, dtype=bool)
    lookup = FisherNegLog10Lookup(labels.size, int(labels.sum()))
    inside = membership.sum(axis=1)
    weights = membership.T.astype(np.float32)

    rng = np.random.default_rng(seed)
    maxima = np.empty(n_permutations)
    for start in range(0, n_permutations, _CHUNK_SIZE):
        size = min(_CHUNK_SIZE, n_permutations - start)
        permuted = np.stack([labels[rng.permutation(labels.size)] for _ in range(size)])
        positive_inside = np.rint(permuted.astype(np.float32) @ weights).astype(np.int64)
        statistics = lookup.values(np.broadcast_to(inside, positive_inside.shape), positive_inside)
        maxima[start : start + size] = statistics.max(axis=1)
    return maxima


def empirical_p_value(observed_statistic: float, null_statistics: np.ndarray) -> float:
    """``(1 + #{null >= observed}) / (n + 1)``, the scans' corrected p-value."""
    exceed = int(np.count_nonzero(null_statistics >= observed_statistic))
    return float((1 + exceed) / (null_statistics.size + 1))
