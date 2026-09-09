"""Numerical reference and boundary cases for the uncorrected CMH test."""

import math

import pytest

from cfh.stats.cochran_mantel_haenszel import cochran_mantel_haenszel

# R's UCBAdmissions, departments A-F: rows admitted/rejected, columns male/female.
# Bickel, Hammel & O'Connell (1975), Science 187, 398-404.
# Published uncorrected R output (X²=1.5246, p=.2169, OR=.9046968):
# https://www.markirwin.net/stat149/Lecture/Lecture8.pdf
BERKELEY = [
    [[512, 89], [313, 19]],
    [[353, 17], [207, 8]],
    [[120, 202], [205, 391]],
    [[138, 131], [279, 244]],
    [[53, 94], [138, 299]],
    [[22, 24], [351, 317]],
]


def test_berkeley_published_example():
    result = cochran_mantel_haenszel(BERKELEY)
    assert result.statistic == pytest.approx(1.5246, abs=0.00005)
    assert result.p_value == pytest.approx(0.2169, abs=0.00005)
    assert result.common_odds_ratio == pytest.approx(0.9046968, abs=0.00000005)
    assert result.informative_strata == 6


def test_uninformative_strata_do_not_change_test():
    reference = cochran_mantel_haenszel(BERKELEY)
    result = cochran_mantel_haenszel(BERKELEY + [[[9, 0], [6, 0]], [[0, 0], [0, 0]]])
    assert result.statistic == reference.statistic
    assert result.p_value == reference.p_value
    assert result.common_odds_ratio == reference.common_odds_ratio
    assert result.total_strata == 8
    assert result.informative_strata == 6


@pytest.mark.parametrize("table", [[[0, 0], [0, 0]], [[1, 0], [0, 0]], [[9, 0], [6, 0]]])
def test_no_information(table):
    result = cochran_mantel_haenszel([table])
    assert result.statistic is None
    assert result.p_value is None
    assert result.common_odds_ratio is None


def test_separation_and_reversed_direction():
    positive = cochran_mantel_haenszel([[[10, 0], [0, 10]]])
    negative = cochran_mantel_haenszel([[[0, 10], [10, 0]]])
    assert positive.common_odds_ratio == math.inf
    assert negative.common_odds_ratio == 0
    assert positive.p_value == negative.p_value
    assert positive.statistic == pytest.approx(19)


def test_cancellation_is_not_concordance():
    result = cochran_mantel_haenszel([[[9, 1], [1, 9]], [[1, 9], [9, 1]]])
    assert result.statistic == 0
    assert result.p_value == 1
    assert result.common_odds_ratio == 1


@pytest.mark.parametrize(
    "tables",
    [
        [],
        [[[1, 2]]],
        [[[1, 2, 3], [4, 5, 6]]],
        [[[True, 0], [0, 1]]],
        [[[-1, 0], [0, 1]]],
        [[[1.5, 0], [0, 1]]],
        [[[float("nan"), 0], [0, 1]]],
    ],
)
def test_invalid_tables(tables):
    with pytest.raises(ValueError):
        cochran_mantel_haenszel(tables)
