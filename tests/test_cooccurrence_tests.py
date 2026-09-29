import pytest

from cfh.stats.cooccurrence_tests import (
    build_cooccurrence_contingency_table,
    fishers_cooccurrence_test,
)


def test_build_cooccurrence_contingency_table_partitions_the_cohort():
    cohort = [f"S{i}" for i in range(1, 11)]
    fusion_positive = {"S1", "S2", "S3"}
    comparator_altered = {"S3", "S4"}

    table = build_cooccurrence_contingency_table(fusion_positive, comparator_altered, cohort)

    # both=S3 (1), fusion_only=S1,S2 (2), comparator_only=S4 (1), neither=6
    assert table == [[1, 2], [1, 6]]


def test_build_cooccurrence_contingency_table_excludes_samples_outside_the_cohort_universe():
    cohort = ["S1", "S2"]
    fusion_positive = {"S1", "S3"}  # S3 not in cohort
    comparator_altered: set[str] = set()

    table = build_cooccurrence_contingency_table(fusion_positive, comparator_altered, cohort)

    # S3 is not eligible; the table total must equal the declared cohort size.
    assert table == [[0, 1], [0, 1]]


def test_fishers_cooccurrence_test_detects_clean_mutual_exclusivity():
    # Every fusion-positive sample lacks the comparator alteration and vice versa.
    cohort = [f"S{i}" for i in range(1, 21)]
    fusion_positive = set(cohort[:8])
    comparator_altered = set(cohort[8:16])

    table = build_cooccurrence_contingency_table(fusion_positive, comparator_altered, cohort)
    odds_ratio, p_value, direction = fishers_cooccurrence_test(table)

    assert odds_ratio == 0.0
    assert p_value < 0.05
    assert direction == "mutually_exclusive"


def test_fishers_cooccurrence_test_detects_clean_cooccurrence():
    cohort = [f"S{i}" for i in range(1, 21)]
    fusion_positive = set(cohort[:8])
    comparator_altered = set(cohort[:8])

    table = build_cooccurrence_contingency_table(fusion_positive, comparator_altered, cohort)
    odds_ratio, p_value, direction = fishers_cooccurrence_test(table)

    assert odds_ratio == float("inf")
    assert p_value < 0.05
    assert direction == "co_occurring"


def test_fishers_cooccurrence_test_is_indeterminate_for_a_perfectly_balanced_table():
    table = [[5, 5], [5, 5]]

    odds_ratio, p_value, direction = fishers_cooccurrence_test(table)

    assert odds_ratio == pytest.approx(1.0)
    assert direction == "indeterminate"


def test_fishers_cooccurrence_test_rejects_non_2x2_table():
    with pytest.raises(ValueError, match="2x2"):
        fishers_cooccurrence_test([[1, 2, 3], [4, 5, 6]])


def test_fishers_cooccurrence_test_rejects_negative_counts():
    with pytest.raises(ValueError, match="non-negative"):
        fishers_cooccurrence_test([[-1, 2], [3, 4]])
