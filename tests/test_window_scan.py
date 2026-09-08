import math

import pytest

from cfh.stats.window_scan import (
    candidate_windows,
    dedup_windows_by_event_mask,
    detect_window,
    scan_windows,
)


def _internal_region_records():
    """Retained breakpoints cluster inside [500, 550]; lost breakpoints sit
    well outside that band on both sides -- an internal region bounded on
    both sides, unlike a single terminal cutpoint."""
    event_ids = [f"r{i}" for i in range(5)] + [f"l{i}" for i in range(5)]
    positions = [505, 510, 520, 530, 545] + [100, 150, 700, 750, 800]
    statuses = ["retained"] * 5 + ["lost"] * 5
    return event_ids, positions, statuses


def _clamped_pile_records():
    """Most events clamp onto a single position (the historical ALK-style
    intronic-clamping artifact) -- many candidate windows around that pile
    contain the exact same set of events."""
    event_ids = [f"c{i}" for i in range(8)] + [f"o{i}" for i in range(4)]
    positions = [400] * 8 + [50, 60, 900, 950]
    statuses = ["retained"] * 8 + ["lost"] * 4
    return event_ids, positions, statuses


def test_candidate_windows_anchors_both_edges_on_observed_positions():
    windows = candidate_windows([100, 300], widths=[50])
    assert (100, 150, 50) in windows
    assert (250, 300, 50) in windows
    assert (300, 350, 50) in windows
    assert (50, 100, 50) in windows


def test_candidate_windows_rejects_non_positive_width():
    with pytest.raises(ValueError):
        candidate_windows([100, 200], widths=[0])


def test_scan_windows_builds_correct_contingency_counts():
    event_ids, positions, statuses = _internal_region_records()
    scan = scan_windows(positions, statuses, event_ids, widths=[50], min_events_per_window=4)

    perfect = next(
        row for row in scan if row["n_positive_inside"] == 5 and row["n_positive_outside"] == 0
    )
    assert perfect["n_events_inside"] == 5
    assert perfect["n_events_outside"] == 5
    assert perfect["odds_ratio"] == float("inf")
    assert set(perfect["event_ids_inside"]) == {"r0", "r1", "r2", "r3", "r4"}


def test_min_events_per_window_guard_excludes_small_windows():
    event_ids, positions, statuses = _internal_region_records()
    scan_loose = scan_windows(positions, statuses, event_ids, widths=[50], min_events_per_window=1)
    scan_strict = scan_windows(positions, statuses, event_ids, widths=[50], min_events_per_window=6)
    assert len(scan_strict) < len(scan_loose)
    assert all(row["n_events_inside"] >= 6 and row["n_events_outside"] >= 6 for row in scan_strict)


def test_mismatched_lengths_raise():
    with pytest.raises(ValueError):
        scan_windows([100, 200], ["retained"], ["a", "b"])


def test_dedup_windows_by_event_mask_collapses_identical_event_sets():
    event_ids, positions, statuses = _clamped_pile_records()
    scan = scan_windows(positions, statuses, event_ids, widths=[25, 50], min_events_per_window=4)

    # Many candidate (start, width) windows contain the exact same 8-event
    # clamped pile -- a naive report would surface all of them as distinct
    # "best windows".
    masks = [frozenset(row["event_ids_inside"]) for row in scan]
    clamped_mask = frozenset(f"c{i}" for i in range(8))
    assert masks.count(clamped_mask) > 1

    deduped = dedup_windows_by_event_mask(scan)
    deduped_masks = [frozenset(row["event_ids_inside"]) for row in deduped]
    assert deduped_masks.count(clamped_mask) == 1
    assert len(deduped) == len(set(masks))


def test_dedup_keeps_the_most_significant_representative_per_mask():
    rows = [
        {
            "start_aa": 100,
            "end_aa": 150,
            "width_aa": 50,
            "neg_log10_p_value": 1.0,
            "event_ids_inside": ("a", "b"),
        },
        {
            "start_aa": 90,
            "end_aa": 140,
            "width_aa": 50,
            "neg_log10_p_value": 3.0,
            "event_ids_inside": ("a", "b"),
        },
    ]
    deduped = dedup_windows_by_event_mask(rows)
    assert len(deduped) == 1
    assert deduped[0]["neg_log10_p_value"] == 3.0


def test_detect_window_recovers_the_internal_region():
    event_ids, positions, statuses = _internal_region_records()

    result = detect_window(positions, statuses, event_ids, seed=42, n_permutations=300)

    assert result["determinable"] is True
    assert result["reason"] is None
    assert result["best_window"]["start_aa"] <= 505
    assert result["best_window"]["end_aa"] >= 545
    assert result["observed_odds_ratio"] == float("inf")
    assert math.isfinite(result["observed_statistic"]) or result["observed_statistic"] == math.inf
    assert 0.0 <= result["corrected_p_value"] <= 1.0
    assert result["scan"]
    assert result["top_windows"]
    assert len(result["top_windows"]) <= len(result["scan"])
    assert len(result["null_max_statistics"]) == 300


def test_detect_window_dedups_top_windows_for_a_clamped_pile():
    event_ids, positions, statuses = _clamped_pile_records()

    result = detect_window(positions, statuses, event_ids, seed=1, n_permutations=200)

    assert result["determinable"] is True
    # The clamped pile produces many numerically distinct but functionally
    # identical candidate windows; top_windows must not repeat the same
    # event-membership mask.
    masks = [frozenset(row["event_ids_inside"]) for row in result["top_windows"]]
    assert len(masks) == len(set(masks))
    assert len(result["top_windows"]) < len(result["scan"])


def test_permutation_correction_is_never_more_optimistic_than_the_naive_scan_minimum():
    event_ids, positions, statuses = _internal_region_records()
    n_permutations = 200

    result = detect_window(positions, statuses, event_ids, seed=7, n_permutations=n_permutations)

    assert result["corrected_p_value"] >= result["observed_p_value"]
    assert result["corrected_p_value"] >= 1 / (n_permutations + 1)


def test_permutation_is_bit_identical_for_a_fixed_seed():
    event_ids, positions, statuses = _internal_region_records()
    first = detect_window(positions, statuses, event_ids, seed=1729, n_permutations=150)
    second = detect_window(positions, statuses, event_ids, seed=1729, n_permutations=150)
    assert first == second


@pytest.mark.parametrize(
    "event_ids,positions,statuses,expected_reason_snippet",
    [
        (["a", "b", "c"], [100, 110, 120], ["retained", "lost", "retained"], "fewer than 4"),
        (
            ["a", "b", "c", "d"],
            [100, 100, 100, 100],
            ["retained", "retained", "lost", "lost"],
            "distinct breakpoint",
        ),
        (
            ["a", "b", "c", "d"],
            [100, 110, 120, 130],
            ["retained"] * 4,
            "single outcome class",
        ),
    ],
)
def test_degenerate_inputs_return_not_determinable_without_raising(
    event_ids, positions, statuses, expected_reason_snippet
):
    result = detect_window(positions, statuses, event_ids, seed=42, n_permutations=100)

    assert result["determinable"] is False
    assert expected_reason_snippet in result["reason"]
    assert result["best_window"] is None
    assert result["corrected_p_value"] is None
    assert result["scan"] == []
    assert result["top_windows"] == []


def test_no_candidate_survives_the_min_events_guard_is_a_graceful_no_op():
    event_ids = ["a", "b", "c", "d"]
    positions = [100, 200, 300, 400]
    statuses = ["retained", "lost", "retained", "lost"]

    result = detect_window(
        positions,
        statuses,
        event_ids,
        seed=42,
        n_permutations=100,
        min_events_per_window=10,
    )

    assert result["determinable"] is False
    assert "min_events_per_window" in result["reason"]


def test_non_positive_permutations_raise():
    event_ids, positions, statuses = _internal_region_records()
    with pytest.raises(ValueError):
        detect_window(positions, statuses, event_ids, n_permutations=0)
