"""Gene-agnostic, frequency/recurrence-only sliding-WINDOW detection.

Generalizes :mod:`cfh.stats.cutpoint_scan` from a single best-separating
*position* to a best-separating *region* ``[start, start + width]``: not
every gene's functionally critical region sits at a protein terminus (as
BRAF/RET/ALK/NTRK1's kinase domains do), so some future gene may need an
internal region bounded on both sides instead of a one-sided cutpoint.

This reuses the exact same max-statistic permutation-correction pattern
:mod:`cfh.stats.cutpoint_scan` already uses (two-sided Fisher's exact test
per candidate, then a label-permutation null over the *max* statistic
across every candidate to correct for having scanned many of them) -- just
extended from a 1-D scan over candidate cutpoints to a 2-D scan over
candidate ``(start, width)`` windows. ``width`` is drawn from a small
predeclared set (:data:`DEFAULT_WIDTHS`) rather than searched continuously,
which keeps the multiple-testing correction tractable and avoids overfitting
window width to noise.

Two real, known data-quality issues this module deliberately defends
against (see the ``directional-intronic-breakpoint-snapping`` fix this
branch is built on):

1. Intronic breakpoints commonly clamp onto a single exon-boundary
   position (historically ~87% of ALK events, before that fix). A naive
   window scan over such data reports a whole *family* of near-identical
   "best windows" around any such clamped pile -- every window that
   happens to contain the same clamped events looks equally significant,
   which looks deceptively precise. :func:`dedup_windows_by_event_mask`
   collapses windows that contain the exact same set of events into one
   candidate before they are surfaced to a human reviewer. This is a
   *reporting* fix, not a statistical one: the permutation-based
   correction below remains valid even when many correlated/duplicate
   candidates are scanned, since it corrects on the empirical distribution
   of the max statistic across the real candidate set, whatever its
   correlation structure.
2. A per-window minimum-event-count guard (``min_events_per_window``)
   excludes candidate windows so small (on either side of the split) that
   their apparent separation is likely sampling noise rather than signal.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

from cfh.stats.cutpoint_scan import fishers_cutpoint_test

DEFAULT_WIDTHS: tuple[int, ...] = (25, 50, 100, 200)
DEFAULT_MIN_EVENTS_PER_WINDOW = 4
_POSITIVE_STATUS_DEFAULT = "retained"


@dataclass(frozen=True)
class _WindowCandidate:
    start: int
    end: int
    width: int
    inside_indices: tuple[int, ...]
    outside_indices: tuple[int, ...]
    event_mask: frozenset[str]


def _neg_log10_p(p_value: float) -> float:
    if p_value <= 0:
        return math.inf
    return -math.log10(p_value)


def candidate_windows(
    positions: Sequence[int], widths: Sequence[int] = DEFAULT_WIDTHS
) -> list[tuple[int, int, int]]:
    """Distinct ``(start, end, width)`` candidate windows for ``positions``.

    For a fixed width, a window's event membership only changes when one of
    its edges crosses an observed breakpoint position, so every candidate is
    anchored so that one edge coincides with an observed position: either
    ``start == p`` (the window opens at an observed breakpoint) or
    ``end == p`` (it closes at one), for every distinct observed ``p``.
    Anchoring anywhere else can never produce a different partition of
    events, so this keeps the candidate set small (``O(n_distinct * n_widths)``)
    without missing any distinct partition.
    """
    distinct = sorted({int(position) for position in positions})
    windows: set[tuple[int, int, int]] = set()
    for width in widths:
        if width <= 0:
            raise ValueError(f"widths must be positive; got {width!r}")
        for position in distinct:
            windows.add((position, position + width, width))
            windows.add((position - width, position, width))
    return sorted(windows)


def _build_candidates(
    positions: Sequence[int],
    event_ids: Sequence[str],
    widths: Sequence[int],
    min_events_per_window: int,
) -> list[_WindowCandidate]:
    candidates = []
    for start, end, width in candidate_windows(positions, widths):
        inside = tuple(
            index for index, position in enumerate(positions) if start <= position <= end
        )
        outside = tuple(index for index in range(len(positions)) if index not in set(inside))
        if len(inside) < min_events_per_window or len(outside) < min_events_per_window:
            continue
        candidates.append(
            _WindowCandidate(
                start=start,
                end=end,
                width=width,
                inside_indices=inside,
                outside_indices=outside,
                event_mask=frozenset(event_ids[index] for index in inside),
            )
        )
    return candidates


def _evaluate_candidate(
    candidate: _WindowCandidate,
    statuses: Sequence[str],
    positive_status: str,
) -> dict:
    n_positive_inside = sum(
        statuses[index] == positive_status for index in candidate.inside_indices
    )
    n_positive_outside = sum(
        statuses[index] == positive_status for index in candidate.outside_indices
    )
    n_inside = len(candidate.inside_indices)
    n_outside = len(candidate.outside_indices)
    table = [
        [n_positive_inside, n_positive_outside],
        [n_inside - n_positive_inside, n_outside - n_positive_outside],
    ]
    odds_ratio, p_value = fishers_cutpoint_test(table)
    return {
        "start_aa": candidate.start,
        "end_aa": candidate.end,
        "width_aa": candidate.width,
        "n_events_inside": n_inside,
        "n_events_outside": n_outside,
        "n_positive_inside": n_positive_inside,
        "n_positive_outside": n_positive_outside,
        "odds_ratio": odds_ratio,
        "p_value": p_value,
        "neg_log10_p_value": _neg_log10_p(p_value),
        "event_ids_inside": tuple(sorted(candidate.event_mask)),
    }


def scan_windows(
    positions: Sequence[int],
    statuses: Sequence[str],
    event_ids: Sequence[str],
    *,
    positive_status: str = _POSITIVE_STATUS_DEFAULT,
    widths: Sequence[int] = DEFAULT_WIDTHS,
    min_events_per_window: int = DEFAULT_MIN_EVENTS_PER_WINDOW,
) -> list[dict]:
    """Compute the separation statistic at every candidate window.

    Returns one row per candidate ``(start, width)`` window that passes the
    ``min_events_per_window`` guard on both sides of the split (empty if no
    candidate does). Each row additionally carries ``event_ids_inside`` --
    the sorted tuple of event ids the window contains -- so callers can
    de-duplicate near-identical windows by event membership rather than by
    their numerically distinct ``(start, width)`` coordinates alone (see
    :func:`dedup_windows_by_event_mask`).
    """
    if not (len(positions) == len(statuses) == len(event_ids)):
        raise ValueError("positions, statuses, and event_ids must be the same length")

    candidates = _build_candidates(positions, event_ids, widths, min_events_per_window)
    return [_evaluate_candidate(candidate, statuses, positive_status) for candidate in candidates]


def _best_scan_row(scan: list[dict]) -> dict:
    """Pick the max-separation row, preferring the narrowest then earliest window on ties."""
    return max(scan, key=lambda row: (row["neg_log10_p_value"], -row["width_aa"], -row["start_aa"]))


def dedup_windows_by_event_mask(scan: list[dict]) -> list[dict]:
    """Collapse windows that contain the exact same set of events into one.

    Real intronic breakpoints commonly clamp onto a single exon-boundary
    position; a naive scan over such data reports a whole family of
    numerically distinct but functionally identical "best windows" around
    that clamped pile. This groups scan rows by ``event_ids_inside`` and
    keeps, per group, the single most-significant representative (highest
    ``neg_log10_p_value``, then narrowest width, then earliest start for a
    deterministic tie-break) -- so a human reviewer sees one candidate per
    distinct event partition, not ``N`` deceptively-precise duplicates.
    Returned sorted by descending significance.
    """
    best_by_mask: dict[frozenset[str], dict] = {}
    for row in scan:
        mask = frozenset(row["event_ids_inside"])
        current = best_by_mask.get(mask)
        if current is None or (
            row["neg_log10_p_value"],
            -row["width_aa"],
            -row["start_aa"],
        ) > (
            current["neg_log10_p_value"],
            -current["width_aa"],
            -current["start_aa"],
        ):
            best_by_mask[mask] = row
    return sorted(
        best_by_mask.values(),
        key=lambda row: (-row["neg_log10_p_value"], row["width_aa"], row["start_aa"]),
    )


def detect_window(
    positions: Sequence[int],
    statuses: Sequence[str],
    event_ids: Sequence[str],
    *,
    positive_status: str = _POSITIVE_STATUS_DEFAULT,
    seed: int = 42,
    n_permutations: int = 10_000,
    widths: Sequence[int] = DEFAULT_WIDTHS,
    min_events_per_window: int = DEFAULT_MIN_EVENTS_PER_WINDOW,
) -> dict:
    """Find the recurrence-based window that best separates ``statuses``.

    Degenerate inputs (too few events, a single distinct breakpoint
    position, a single outcome class, or no candidate window surviving the
    ``min_events_per_window`` guard) never raise: they return
    ``determinable: False`` with a human-readable ``reason`` instead, so
    callers can surface a clear "not determinable" (no-op) result.

    Returns a plain ``dict`` (JSON-serializable) with keys: ``determinable``,
    ``reason``, ``best_window`` (``{"start_aa", "end_aa", "width_aa"}`` or
    ``None``), ``observed_statistic`` (``-log10(p)`` at the best window),
    ``observed_p_value``, ``observed_odds_ratio``, ``corrected_p_value``
    (the permutation-based empirical p-value), ``scan`` (the full
    per-window scan, post min-events guard), ``top_windows`` (``scan``
    de-duplicated by event membership -- see
    :func:`dedup_windows_by_event_mask`), and ``null_max_statistics``
    (present only when determinable).

    The set of candidate windows (and which windows survive the
    ``min_events_per_window`` guard) depends only on ``positions`` and
    ``widths``, never on ``statuses`` -- so it is computed once and reused,
    unchanged, for every permutation replicate below.
    """
    if n_permutations <= 0:
        raise ValueError("n_permutations must be positive")
    if not (len(positions) == len(statuses) == len(event_ids)):
        raise ValueError("positions, statuses, and event_ids must be the same length")

    positions = [int(position) for position in positions]
    statuses = list(statuses)
    event_ids = list(event_ids)

    reason: str | None = None
    if len(positions) < 4:
        reason = "fewer than 4 mapped breakpoint events with a known outcome status"
    elif len(set(positions)) < 2:
        reason = "fewer than 2 distinct breakpoint positions to scan"
    elif len(set(statuses)) < 2:
        reason = "all events share a single outcome class; no separation is possible"

    candidates = (
        _build_candidates(positions, event_ids, widths, min_events_per_window)
        if reason is None
        else []
    )
    if reason is None and not candidates:
        reason = (
            "no candidate window of the configured widths "
            f"{tuple(widths)} met the minimum-event-count guard "
            f"(min_events_per_window={min_events_per_window}) on both sides of the split"
        )

    if reason is not None:
        return {
            "determinable": False,
            "reason": reason,
            "best_window": None,
            "observed_statistic": None,
            "observed_p_value": None,
            "observed_odds_ratio": None,
            "corrected_p_value": None,
            "scan": [],
            "top_windows": [],
        }

    scan = [_evaluate_candidate(candidate, statuses, positive_status) for candidate in candidates]
    best = _best_scan_row(scan)
    observed_statistic = best["neg_log10_p_value"]

    rng = np.random.default_rng(seed)
    null_max_statistics: list[float] = []
    for _ in range(n_permutations):
        permuted_statuses = rng.permutation(np.asarray(statuses, dtype=object)).tolist()
        permuted_max = max(
            _evaluate_candidate(candidate, permuted_statuses, positive_status)["neg_log10_p_value"]
            for candidate in candidates
        )
        null_max_statistics.append(permuted_max)

    corrected_p_value = (1 + sum(stat >= observed_statistic for stat in null_max_statistics)) / (
        n_permutations + 1
    )

    return {
        "determinable": True,
        "reason": None,
        "best_window": {
            "start_aa": best["start_aa"],
            "end_aa": best["end_aa"],
            "width_aa": best["width_aa"],
        },
        "observed_statistic": observed_statistic,
        "observed_p_value": best["p_value"],
        "observed_odds_ratio": best["odds_ratio"],
        "corrected_p_value": float(corrected_p_value),
        "scan": scan,
        "top_windows": dedup_windows_by_event_mask(scan),
        "null_max_statistics": tuple(null_max_statistics),
    }
