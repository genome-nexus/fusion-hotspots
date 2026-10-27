"""How many patients hold a significant 2x2 result up?

A removal fragility index: the fewest observations that, taken out of the
cells supporting the association, raise the p-value to ``alpha`` or above.
With one observation per patient per gene, that is a count of patients. A
result that falls apart when two or three patients are removed rests on
those patients, however small its p-value; FGFR2's tumor-type-adjusted
retention (CMH p 0.03) is the motivating case -- it went to 0.30 once six
patients with two distinct FGFR2 fusions were counted once.

Only supporting cells are emptied: in-frame with the hit status and
other-frame without it for a positive association, the off-diagonal for a
negative one. Removing an observation from a non-supporting cell can only
strengthen the association. The search is greedy -- each step removes the
one observation, from any stratum, that raises the p-value most -- so the
index is an upper bound on the true minimum: a result reported as fragile
is at least that fragile.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from typing import Any

from cfh.stats.breakpoint_tests import fishers_frame_domain_test
from cfh.stats.cochran_mantel_haenszel import cochran_mantel_haenszel

DEFAULT_ALPHA = 0.05
_SUPPORTING_CELLS = {
    "positive": ((0, 0), (1, 1)),
    "negative": ((0, 1), (1, 0)),
}
_CELL_LABELS = {
    (0, 0): "hit_in_frame",
    (0, 1): "hit_other_frame",
    (1, 0): "non_hit_in_frame",
    (1, 1): "non_hit_other_frame",
}

PValueFunction = Callable[[list[list[list[int]]]], float | None]


def fisher_p_value(tables: list[list[list[int]]]) -> float:
    """One-sided pooled Fisher p-value of a single-table list."""
    (table,) = tables
    return fishers_frame_domain_test(table)[1]


def cmh_p_value(tables: list[list[list[int]]]) -> float | None:
    """Two-sided CMH p-value across strata."""
    return cochran_mantel_haenszel(tables).p_value


def removal_fragility(
    tables: Sequence[Sequence[Sequence[int]]],
    p_value_of: PValueFunction,
    *,
    direction: str = "positive",
    alpha: float = DEFAULT_ALPHA,
    labels: Sequence[str] | None = None,
) -> dict[str, Any]:
    """Greedily remove supporting observations until ``p >= alpha``.

    ``tables`` holds one 2x2 table per stratum (a single table for the pooled
    Fisher test), oriented as in
    :func:`~cfh.stats.breakpoint_tests.build_frame_domain_contingency_table`.
    ``labels`` names the strata in the returned removals. A p-value that
    becomes undefined (every stratum uninformative) counts as non-significant.
    """
    if direction not in _SUPPORTING_CELLS:
        raise ValueError(f"direction must be one of {sorted(_SUPPORTING_CELLS)}")
    current = [[list(row) for row in table] for table in tables]
    labels = list(labels) if labels is not None else [str(i) for i in range(len(current))]
    n = sum(sum(sum(row) for row in table) for table in current)
    p_value = p_value_of(current)
    result: dict[str, Any] = {
        "alpha": alpha,
        "direction": direction,
        "n_observations": n,
        "p_value": p_value,
    }
    if p_value is None or p_value >= alpha:
        return {**result, "significant": False, "fragility_index": None}

    removals: list[tuple[int, tuple[int, int]]] = []
    while p_value is not None and p_value < alpha:
        best: tuple[float, int, tuple[int, int]] | None = None
        for index, table in enumerate(current):
            for row, column in _SUPPORTING_CELLS[direction]:
                if not table[row][column]:
                    continue
                table[row][column] -= 1
                candidate = p_value_of(current)
                table[row][column] += 1
                candidate = 1.0 if candidate is None else candidate
                if best is None or candidate > best[0]:
                    best = (candidate, index, (row, column))
        if best is None:
            # Every supporting cell is empty yet p < alpha; cannot happen for
            # a well-oriented test, but never loop forever.
            break
        p_value, index, (row, column) = best
        current[index][row][column] -= 1
        removals.append((index, (row, column)))

    removed: dict[tuple[str, str], int] = {}
    for index, cell in removals:
        key = (labels[index], _CELL_LABELS[cell])
        removed[key] = removed.get(key, 0) + 1
    return {
        **result,
        "significant": True,
        "fragility_index": len(removals),
        "fragility_quotient": len(removals) / n if n else None,
        "p_value_after_removal": p_value,
        "removed": [
            {"stratum": stratum, "cell": cell, "count": count}
            for (stratum, cell), count in removed.items()
        ],
    }


def fisher_fragility(table: Sequence[Sequence[int]], *, alpha: float = DEFAULT_ALPHA) -> dict:
    """Fragility of the pooled one-sided Fisher test."""
    return removal_fragility([table], fisher_p_value, alpha=alpha, labels=["pooled"])


def cmh_fragility(
    tables_by_stratum: Mapping[str, Sequence[Sequence[int]]],
    direction: str | None,
    *,
    alpha: float = DEFAULT_ALPHA,
) -> dict | None:
    """Fragility of the tumor-type-stratified CMH test, or ``None`` when the
    association has no direction (no odds ratio or an odds ratio of 1)."""
    if direction not in _SUPPORTING_CELLS:
        return None
    return removal_fragility(
        list(tables_by_stratum.values()),
        cmh_p_value,
        direction=direction,
        alpha=alpha,
        labels=list(tables_by_stratum),
    )
