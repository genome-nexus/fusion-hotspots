"""Cohort-wide mutation/CNA co-occurrence vs mutual-exclusivity Fisher's-exact test.

Distinct from :mod:`cfh.stats.breakpoint_tests` (which tests domain
retention among a gene's OWN in-frame fusions): this tests whether a gene's
fusion-positive sample set co-occurs with, or is mutually exclusive with, a
configured comparator alteration ACROSS THE FULL COHORT -- every profiled
sample, not just the gene's own fusion-positive ones -- the standard
two-way mutation co-occurrence design (e.g. BRAF fusions vs BRAF V600E
point mutations).

A two-sided Fisher's exact test on the resulting 2x2 table is the
pre-specified statistic here: unlike domain retention (where the tested
direction -- retention enriched among in-frame fusions -- is fixed by
biology ahead of time), co-occurrence vs mutual exclusivity is not assumed
in either direction, so the observed odds ratio's direction is reported
alongside the two-sided p-value rather than baked into a one-sided
alternative.
"""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np
from scipy.stats import fisher_exact


def build_cooccurrence_contingency_table(
    fusion_positive_sample_ids: Iterable[str],
    comparator_altered_sample_ids: Iterable[str],
    cohort_sample_ids: Iterable[str],
) -> list[list[int]]:
    """Build the 2x2 table ``[[both, fusion_only], [comparator_only, neither]]``.

    Rows are fusion status (positive, negative); columns are the configured
    comparator's alteration status (altered, not). The table's total is the
    size of the union of ``cohort_sample_ids`` with the two alteration sets
    -- in the expected case both alteration sets are subsets of the cohort
    universe and the total is simply the cohort size, but a sample outside
    the supplied cohort universe still contributes to its own row/column
    rather than being silently dropped.
    """
    cohort = set(cohort_sample_ids)
    fusion_positive = set(fusion_positive_sample_ids)
    comparator_altered = set(comparator_altered_sample_ids)
    universe = cohort | fusion_positive | comparator_altered

    both = len(fusion_positive & comparator_altered)
    fusion_only = len(fusion_positive - comparator_altered)
    comparator_only = len(comparator_altered - fusion_positive)
    neither = len(universe) - both - fusion_only - comparator_only
    return [[both, fusion_only], [comparator_only, neither]]


def fishers_cooccurrence_test(
    contingency_table: Iterable[Iterable[int]],
) -> tuple[float, float, str]:
    """Run the pre-specified two-sided Fisher exact test.

    Returns ``(odds_ratio, p_value, direction)``, where ``direction`` is
    ``"co_occurring"`` (odds ratio > 1), ``"mutually_exclusive"`` (odds
    ratio < 1), or ``"indeterminate"`` (odds ratio exactly 1, e.g. an
    all-zero table). A zero off-diagonal cell can still produce an infinite
    or zero odds ratio; SciPy still computes the exact finite p-value for
    that table.
    """
    table = np.asarray(list(contingency_table), dtype=int)
    if table.shape != (2, 2):
        raise ValueError("contingency_table must be a 2x2 table")
    if np.any(table < 0):
        raise ValueError("contingency_table counts must be non-negative")

    odds_ratio, p_value = fisher_exact(table, alternative="two-sided")
    if odds_ratio > 1:
        direction = "co_occurring"
    elif odds_ratio < 1:
        direction = "mutually_exclusive"
    else:
        direction = "indeterminate"
    return float(odds_ratio), float(p_value), direction
