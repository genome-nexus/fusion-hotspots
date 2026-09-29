"""Hypothesis-free fusion x alteration co-occurrence discovery.

``mutation_cooccurrence`` tests comparators named in a gene's curated config.
This module instead tests every fusion gene in a cohort scan against every
gene x alteration type (mutation, amplification, deep deletion) with enough
altered samples, and reports both co-occurring and mutually exclusive pairs.

Design choices, each addressing a known confounder of co-occurrence tests:

* **Assay eligibility.** A pair is tested only on samples whose SV panel
  covers the fusion gene *and* whose mutation/CNA panel covers the comparator
  gene. Eligibility depends only on the panel a sample was sequenced with, so
  counts are precomputed per panel group, which keeps the scan vectorized.
* **Tumor type.** The primary statistic is the Cochran-Mantel-Haenszel test
  stratified by OncoTree code, with the Mantel-Haenszel common odds ratio.
  Pooled counts and odds ratios are reported for comparison only.
* **Hypermutation.** Callers may exclude samples above a mutation-count
  threshold, since hypermutated tumors co-occur with everything.
* **Multiplicity.** All tested pairs form one Benjamini-Hochberg family,
  separate from the cohort scan's main FDR family and from curated targets.

Counting unit: samples. Repeated samples from one patient are not collapsed.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy.stats import chi2

from cfh.stats.multiple_testing import benjamini_hochberg

ALTERATION_TYPES = ("mutation", "amplification", "deep_deletion")
DEFAULT_MIN_FUSION_SAMPLES = 5
DEFAULT_MIN_COMPARATOR_SAMPLES = 20
MANTEL_FLEISS_MARGIN = 5.0
"""The CMH chi-square approximation is used only when the summed expected
``both`` count is at least this far from both of its attainable bounds
(Mantel & Fleiss, 1980). Sparse pairs -- e.g. one shared sample spread over
hundreds of tumor-type strata -- otherwise get wildly anti-conservative
p-values."""


@dataclass(frozen=True)
class SampleRecord:
    sample_id: str
    stratum: str
    sv_panel: str | None
    mutation_panel: str | None
    cna_panel: str | None


def _panel_field(alteration_type: str) -> str:
    if alteration_type == "mutation":
        return "mutation_panel"
    if alteration_type in {"amplification", "deep_deletion"}:
        return "cna_panel"
    raise ValueError(f"Unknown alteration type {alteration_type!r}")


def discover_cooccurrence(
    samples: Iterable[SampleRecord],
    panel_genes: Mapping[str, Iterable[str]],
    fusion_samples: Mapping[str, Iterable[str]],
    alterations: Mapping[tuple[str, str], Iterable[str]],
    *,
    min_fusion_samples: int = DEFAULT_MIN_FUSION_SAMPLES,
    min_comparator_samples: int = DEFAULT_MIN_COMPARATOR_SAMPLES,
    excluded_samples: Iterable[str] = (),
) -> dict[str, Any]:
    """Test every fusion gene against every comparator; return rows and counts.

    ``fusion_samples`` maps a fusion gene symbol to its fusion-positive sample
    IDs. ``alterations`` maps ``(gene symbol, alteration type)`` to altered
    sample IDs. ``panel_genes`` maps a panel ID to the gene symbols it covers.
    Rows are returned for every tested pair with a BH q-value over the CMH
    p-values of all tested pairs.
    """
    if min_fusion_samples < 1 or min_comparator_samples < 1:
        raise ValueError("minimum sample counts must be positive")
    panels = {
        panel: frozenset(gene.upper() for gene in genes) for panel, genes in panel_genes.items()
    }
    excluded = set(excluded_samples)
    records = [record for record in samples if record.sample_id not in excluded]

    strata = sorted({record.stratum for record in records})
    stratum_index = {stratum: index for index, stratum in enumerate(strata)}
    group_keys = sorted(
        {(r.sv_panel, r.mutation_panel, r.cna_panel) for r in records},
        key=lambda key: tuple("" if value is None else value for value in key),
    )
    group_index = {key: index for index, key in enumerate(group_keys)}
    sample_position = {
        record.sample_id: (
            group_index[(record.sv_panel, record.mutation_panel, record.cna_panel)],
            stratum_index[record.stratum],
        )
        for record in records
    }
    n_groups, n_strata = len(group_keys), len(strata)
    group_sizes = np.zeros((n_groups, n_strata))
    for g, t in sample_position.values():
        group_sizes[g, t] += 1

    def covers(panel: str | None, gene: str) -> bool:
        return panel is not None and gene.upper() in panels.get(panel, frozenset())

    comparators = sorted(alterations)
    n_comp = len(comparators)
    comparator_eligible = np.zeros((n_comp, n_groups), dtype=bool)
    comparator_counts = np.zeros((n_comp, n_groups, n_strata))
    altered_by_sample: dict[str, list[int]] = {}
    for j, (gene, alteration_type) in enumerate(comparators):
        field_index = ("sv_panel", "mutation_panel", "cna_panel").index(
            _panel_field(alteration_type)
        )
        for g, key in enumerate(group_keys):
            comparator_eligible[j, g] = covers(key[field_index], gene)
        for sample_id in set(alterations[(gene, alteration_type)]):
            position = sample_position.get(sample_id)
            if position is None:
                continue
            comparator_counts[j, position[0], position[1]] += 1
            altered_by_sample.setdefault(sample_id, []).append(j)

    rows: list[dict[str, Any]] = []
    skipped = {"fusion_below_minimum": 0, "pairs_below_minimum": 0}
    for fusion_gene in sorted(fusion_samples):
        fusion_eligible = np.array([covers(key[0], fusion_gene) for key in group_keys])
        fusion_counts = np.zeros((n_groups, n_strata))
        positives = []
        for sample_id in set(fusion_samples[fusion_gene]):
            position = sample_position.get(sample_id)
            if position is None or not fusion_eligible[position[0]]:
                continue
            fusion_counts[position] += 1
            positives.append((sample_id, position))
        if len(positives) < min_fusion_samples:
            skipped["fusion_below_minimum"] += 1
            continue

        universe = comparator_eligible & fusion_eligible[None, :]  # J x G
        weights = universe.astype(float)
        n = weights @ group_sizes  # J x T
        m1 = weights @ fusion_counts
        m2 = np.einsum("jg,jgt->jt", weights, comparator_counts)
        a = np.zeros((n_comp, n_strata))
        for sample_id, (g, t) in positives:
            for j in altered_by_sample.get(sample_id, ()):
                if universe[j, g]:
                    a[j, t] += 1

        testable = (m1.sum(axis=1) >= min_fusion_samples) & (
            m2.sum(axis=1) >= min_comparator_samples
        )
        skipped["pairs_below_minimum"] += int(n_comp - testable.sum())
        if not testable.any():
            continue
        statistic, p_value, mh_or, informative, adequate = _cmh(
            a[testable], m1[testable], m2[testable], n[testable]
        )
        for k, j in enumerate(int(index) for index in np.flatnonzero(testable)):
            gene, alteration_type = comparators[j]
            both = int(a[j].sum())
            fusion_total = int(m1[j].sum())
            comparator_total = int(m2[j].sum())
            universe_total = int(n[j].sum())
            pooled_table = [
                [both, fusion_total - both],
                [comparator_total - both, universe_total - fusion_total - comparator_total + both],
            ]
            rows.append(
                {
                    "fusion_gene": fusion_gene,
                    "comparator_gene": gene,
                    "alteration_type": alteration_type,
                    "same_gene": gene.upper() == fusion_gene.upper(),
                    "eligible_samples": universe_total,
                    "fusion_samples": fusion_total,
                    "comparator_samples": comparator_total,
                    "both_samples": both,
                    "pooled_table": pooled_table,
                    "pooled_odds_ratio": _odds_ratio(pooled_table),
                    "expected_both_samples": float(
                        np.divide(m1[j] * m2[j], n[j], out=np.zeros(n_strata), where=n[j] > 0).sum()
                    ),
                    "mh_common_odds_ratio": _json_number(mh_or[k]),
                    "cmh_statistic": _json_number(statistic[k]),
                    "mantel_fleiss_satisfied": bool(adequate[k]),
                    # No p-value when the chi-square approximation is invalid.
                    "cmh_p_value": _json_number(p_value[k]) if adequate[k] else None,
                    "informative_strata": int(informative[k]),
                    "direction": _direction(mh_or[k]),
                }
            )

    tested = [row for row in rows if row["cmh_p_value"] is not None]
    adjusted = benjamini_hochberg(
        [
            (
                row["fusion_gene"],
                f"{row['comparator_gene']}:{row['alteration_type']}",
                row["cmh_p_value"],
            )
            for row in tested
        ]
    )
    for row, (_, _, _, q_value) in zip(tested, adjusted, strict=True):
        row["cmh_q_value"] = q_value
    for row in rows:
        row.setdefault("cmh_q_value", None)
    return {
        "rows": rows,
        "tested_pair_count": len(tested),
        "untestable_pair_count": len(rows) - len(tested),
        "sample_count": len(records),
        "excluded_sample_count": len(excluded),
        "stratum_count": n_strata,
        "panel_group_count": n_groups,
        "comparator_count": n_comp,
        "skipped": skipped,
        "stratified_by": "oncotree_code",
        "counting_unit": "sample",
    }


def _cmh(a: np.ndarray, m1: np.ndarray, m2: np.ndarray, n: np.ndarray):
    """Vectorized CMH (no continuity correction) and MH odds ratio per row.

    Also returns whether each row meets the Mantel-Fleiss criterion.
    """
    b = m1 - a
    c = m2 - a
    d = n - m1 - m2 + a
    with np.errstate(divide="ignore", invalid="ignore"):
        informative_mask = (n > 1) & (m1 > 0) & (m1 < n) & (m2 > 0) & (m2 < n)
        expected = np.where(informative_mask, m1 * m2 / np.where(n > 0, n, 1), 0.0)
        variance = np.where(
            informative_mask,
            m1 * (n - m1) * m2 * (n - m2) / np.where(n > 1, n * n * (n - 1), 1),
            0.0,
        )
        residual = np.where(informative_mask, a - expected, 0.0).sum(axis=1)
        total_variance = variance.sum(axis=1)
        statistic = np.where(
            total_variance > 0,
            residual**2 / np.where(total_variance > 0, total_variance, 1),
            np.nan,
        )
        p_value = np.where(np.isnan(statistic), np.nan, chi2.sf(np.nan_to_num(statistic), 1))
        safe_n = np.where(n > 0, n, 1)
        numerator = np.where(n > 0, a * d / safe_n, 0.0).sum(axis=1)
        denominator = np.where(n > 0, b * c / safe_n, 0.0).sum(axis=1)
        mh_or = np.where(
            denominator > 0,
            numerator / np.where(denominator > 0, denominator, 1),
            np.where(numerator > 0, np.inf, np.nan),
        )
        expected_total = expected.sum(axis=1)
        lower = np.where(informative_mask, np.maximum(0, m1 + m2 - n), 0).sum(axis=1)
        upper = np.where(informative_mask, np.minimum(m1, m2), 0).sum(axis=1)
        adequate = (expected_total - lower >= MANTEL_FLEISS_MARGIN) & (
            upper - expected_total >= MANTEL_FLEISS_MARGIN
        )
    return statistic, p_value, mh_or, informative_mask.sum(axis=1), adequate


def _odds_ratio(table: list[list[int]]) -> float | str | None:
    (a, b), (c, d) = table
    if b * c:
        return a * d / (b * c)
    return "infinity" if a * d else None


def _json_number(value: float) -> float | str | None:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return None
    if math.isinf(value):
        return "infinity"
    return float(value)


def _direction(odds_ratio: float) -> str | None:
    if math.isnan(odds_ratio):
        return None
    if odds_ratio > 1:
        return "co_occurring"
    if odds_ratio < 1:
        return "mutually_exclusive"
    return "none"
