"""Removal fragility of the frame x domain Fisher and CMH tests."""

from __future__ import annotations

import itertools
import json

import pytest

from cfh.stats.breakpoint_tests import fishers_frame_domain_test
from cfh.stats.cochran_mantel_haenszel import cochran_mantel_haenszel
from cfh.stats.fragility import cmh_fragility, fisher_fragility, removal_fragility
from cfh.stats.stratification import stratified_cmh_summary


def _exact_fisher_fragility(table: list[list[int]]) -> int | None:
    """Brute-force minimum over every split of removals between the two
    supporting cells; the oracle the greedy search is checked against."""
    (a, b), (c, d) = table
    if fishers_frame_domain_test(table)[1] >= 0.05:
        return None
    for k in range(1, a + d + 1):
        for from_a in range(min(k, a) + 1):
            from_d = k - from_a
            if from_d > d:
                continue
            if fishers_frame_domain_test([[a - from_a, b], [c, d - from_d]])[1] >= 0.05:
                return k
    return None


# Pooled retention tables from the committed 26 October MSK-IMPACT 50k run.
@pytest.mark.parametrize(
    ("gene", "table", "expected"),
    [
        ("ALK", [[186, 29], [3, 6]], 4),
        ("FGFR2", [[88, 16], [7, 8]], 4),
        ("RET", [[112, 33], [4, 8]], 4),
        ("ETV6", [[38, 8], [5, 11]], 7),
        ("BRAF", [[129, 19], [9, 5]], 1),
    ],
)
def test_real_retention_tables(gene, table, expected):
    result = fisher_fragility(table)
    assert result["significant"]
    assert result["fragility_index"] == expected == _exact_fisher_fragility(table)
    assert result["p_value_after_removal"] >= 0.05
    assert sum(row["count"] for row in result["removed"]) == expected
    assert result["n_observations"] == sum(map(sum, table))


def test_greedy_matches_exact_on_small_tables():
    for cells in itertools.product(range(0, 9, 2), repeat=4):
        table = [[cells[0], cells[1]], [cells[2], cells[3]]]
        assert fisher_fragility(table)["fragility_index"] == _exact_fisher_fragility(table), table


def test_non_significant_result_has_no_index():
    result = fisher_fragility([[88, 20], [7, 4]])
    assert not result["significant"]
    assert result["fragility_index"] is None


def test_cmh_removes_from_the_supporting_diagonal_of_any_stratum():
    tables = {"LUAD": [[20, 2], [1, 6]], "THPA": [[10, 3], [2, 4]]}
    result = cmh_fragility(tables, "positive")
    assert result["significant"] and result["fragility_index"] >= 1
    assert {row["cell"] for row in result["removed"]} <= {"hit_in_frame", "non_hit_other_frame"}
    assert {row["stratum"] for row in result["removed"]} <= set(tables)
    remaining = [[row[:] for row in table] for table in tables.values()]
    for row in result["removed"]:
        index = list(tables).index(row["stratum"])
        cell = {"hit_in_frame": (0, 0), "non_hit_other_frame": (1, 1)}[row["cell"]]
        remaining[index][cell[0]][cell[1]] -= row["count"]
    assert cochran_mantel_haenszel(remaining).p_value >= 0.05


def test_negative_association_removes_off_diagonal():
    result = removal_fragility(
        [[[2, 20], [20, 2]]],
        lambda tables: cochran_mantel_haenszel(tables).p_value,
        direction="negative",
    )
    assert result["significant"]
    assert {row["cell"] for row in result["removed"]} <= {"hit_other_frame", "non_hit_in_frame"}


def test_cmh_without_direction_has_no_fragility():
    assert cmh_fragility({"A": [[1, 1], [1, 1]]}, "none") is None
    assert cmh_fragility({"A": [[1, 1], [1, 1]]}, None) is None


def test_stratified_summary_reports_fragility_as_strict_json():
    summary = stratified_cmh_summary({"A": [[20, 1], [1, 20]], "B": [[0, 0], [0, 0]]})
    assert summary["fragility"]["significant"]
    assert summary["fragility"]["fragility_index"] >= 1
    json.dumps(summary, allow_nan=False)
    assert stratified_cmh_summary({"A": [[0, 0], [0, 0]]})["fragility"] is None
