"""CMH conditional association test and Mantel-Haenszel common odds ratio."""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass
from numbers import Integral

from scipy.stats import chi2


@dataclass(frozen=True)
class CMHResult:
    statistic: float | None
    p_value: float | None
    common_odds_ratio: float | None
    informative_strata: int
    total_strata: int


def cochran_mantel_haenszel(tables: Sequence[Sequence[Sequence[int]]]) -> CMHResult:
    """Test conditional independence in consistently oriented 2x2 tables.

    For [[a,b],[c,d]], E[a]=(a+b)(a+c)/n and
    Var[a]=(a+b)(c+d)(a+c)(b+d)/(n²(n-1)). The statistic is
    (sum(a-E[a]))²/sum(Var[a]), with a chi-square(1) null distribution.
    No continuity correction or pseudocounts are applied. The MH estimator
    is sum(ad/n)/sum(bc/n). Independent observations within and between
    strata are assumed; this is not a test of homogeneity of odds ratios.

    Empty or fixed-margin strata contribute zero information. If every
    stratum is uninformative, statistic and p-value are None. An undefined
    odds ratio (0/0) is None; a positive numerator over zero is infinity.
    """
    if not tables:
        raise ValueError("At least one 2x2 table is required")
    residuals, variances, numerators, denominators = [], [], [], []
    for table in tables:
        if len(table) != 2 or any(len(row) != 2 for row in table):
            raise ValueError("Each table must have shape 2x2")
        cells = [value for row in table for value in row]
        if any(isinstance(v, bool) or not isinstance(v, Integral) or v < 0 for v in cells):
            raise ValueError("Table counts must be nonnegative integers")
        a, b, c, d = map(int, cells)
        n = a + b + c + d
        if n == 0:
            continue
        numerators.append(a * d / n)
        denominators.append(b * c / n)
        margins = (a + b) * (c + d) * (a + c) * (b + d)
        if n <= 1 or margins == 0:
            continue
        residuals.append(a - (a + b) * (a + c) / n)
        variances.append(margins / (n * n * (n - 1)))
    numerator, denominator = math.fsum(numerators), math.fsum(denominators)
    odds = numerator / denominator if denominator else (math.inf if numerator else None)
    variance = math.fsum(variances)
    statistic = math.fsum(residuals) ** 2 / variance if variance else None
    return CMHResult(
        statistic=statistic,
        p_value=float(chi2.sf(statistic, 1)) if statistic is not None else None,
        common_odds_ratio=odds,
        informative_strata=len(variances),
        total_strata=len(tables),
    )
