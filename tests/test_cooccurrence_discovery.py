"""Hypothesis-free co-occurrence discovery across fusion genes and comparators."""

from __future__ import annotations

import math

import numpy as np
import pytest

from cfh.cohort.cooccurrence_discovery import SampleRecord, discover_cooccurrence
from cfh.stats.cochran_mantel_haenszel import cochran_mantel_haenszel

PANELS = {"P1": {"FUS", "MUT", "OTHER"}, "P2": {"FUS", "OTHER"}}


def _samples(n, stratum, panel="P1", prefix=None):
    prefix = prefix or f"{stratum}-{panel}"
    return [SampleRecord(f"{prefix}-{i}", stratum, panel, panel, panel) for i in range(n)]


def test_vectorized_cmh_matches_scalar_reference():
    rng = np.random.default_rng(0)
    samples = _samples(300, "LUAD") + _samples(200, "SKCM") + _samples(100, "PRAD")
    ids = [s.sample_id for s in samples]
    fusion = {sid for sid in ids if rng.random() < 0.15}
    mutated = {sid for sid in ids if rng.random() < 0.3}
    report = discover_cooccurrence(samples, PANELS, {"FUS": fusion}, {("MUT", "mutation"): mutated})
    (row,) = report["rows"]
    tables = []
    for stratum in ("LUAD", "PRAD", "SKCM"):
        members = {s.sample_id for s in samples if s.stratum == stratum}
        f, m = fusion & members, mutated & members
        tables.append([[len(f & m), len(f - m)], [len(m - f), len(members - f - m)]])
    reference = cochran_mantel_haenszel(tables)
    assert row["cmh_statistic"] == pytest.approx(reference.statistic)
    assert row["cmh_p_value"] == pytest.approx(reference.p_value)
    assert row["mh_common_odds_ratio"] == pytest.approx(reference.common_odds_ratio)
    assert row["cmh_q_value"] == pytest.approx(row["cmh_p_value"])  # family of one


def test_tumor_type_confounding_is_removed():
    # FUS and MUT are both common in LUAD and rare in SKCM, independently
    # within each tumor type, so the pooled table looks strongly co-occurring.
    samples = _samples(100, "LUAD") + _samples(100, "SKCM")
    luad = [s.sample_id for s in samples if s.stratum == "LUAD"]
    skcm = [s.sample_id for s in samples if s.stratum == "SKCM"]
    fusion = set(luad[:50]) | set(skcm[:5])
    mutated = set(luad[25:75]) | set(skcm[5:10])  # half of LUAD fusions mutated
    report = discover_cooccurrence(
        samples, PANELS, {"FUS": fusion}, {("MUT", "mutation"): mutated}, min_comparator_samples=5
    )
    (row,) = report["rows"]
    assert row["pooled_odds_ratio"] > 1.5
    assert row["cmh_p_value"] > 0.5


def test_only_jointly_assayed_samples_enter_the_universe():
    # MUT is not on P2, so P2 samples cannot count as MUT-negative.
    samples = _samples(40, "LUAD", "P1") + _samples(60, "LUAD", "P2")
    p2 = {s.sample_id for s in samples if s.sv_panel == "P2"}
    p1 = [s.sample_id for s in samples if s.sv_panel == "P1"]
    report = discover_cooccurrence(
        samples,
        PANELS,
        {"FUS": set(p1[:10]) | p2},
        {("MUT", "mutation"): set(p1[5:30])},
        min_comparator_samples=5,
    )
    (row,) = report["rows"]
    assert row["eligible_samples"] == 40
    assert row["fusion_samples"] == 10
    assert row["both_samples"] == 5


def test_thresholds_exclusions_and_labels():
    samples = _samples(50, "LUAD")
    ids = [s.sample_id for s in samples]
    report = discover_cooccurrence(
        samples,
        PANELS,
        {"FUS": set(ids[:10]), "RARE": set(ids[:2])},
        {("FUS", "amplification"): set(ids[5:30]), ("OTHER", "mutation"): set(ids[:3])},
        min_comparator_samples=5,
        excluded_samples=ids[:1],
    )
    assert report["skipped"]["fusion_below_minimum"] == 1  # RARE
    assert report["excluded_sample_count"] == 1
    (row,) = report["rows"]  # OTHER has too few altered samples
    assert row["same_gene"] is True
    assert row["fusion_samples"] == 9
    with pytest.raises(ValueError):
        discover_cooccurrence(samples, PANELS, {}, {("X", "fusion"): set()})


def test_complete_exclusivity_reports_zero_odds_ratio_and_direction():
    samples = _samples(100, "LUAD")
    ids = [s.sample_id for s in samples]
    report = discover_cooccurrence(
        samples, PANELS, {"FUS": set(ids[:30])}, {("MUT", "mutation"): set(ids[30:60])}
    )
    (row,) = report["rows"]
    assert row["mh_common_odds_ratio"] == 0.0
    assert row["direction"] == "mutually_exclusive"
    assert row["cmh_p_value"] < 1e-4
    assert not any(isinstance(v, float) and math.isnan(v) for v in row.values())


def test_sparse_pairs_fail_mantel_fleiss_and_leave_the_fdr_family():
    # One shared sample, spread over many small tumor-type strata.
    samples = [SampleRecord(f"S{i}", f"T{i % 40}", "P1", "P1", "P1") for i in range(400)]
    ids = [s.sample_id for s in samples]
    report = discover_cooccurrence(
        samples,
        PANELS,
        {"FUS": set(ids[:5])},
        {("MUT", "mutation"): {ids[0], *ids[100:120]}},
    )
    (row,) = report["rows"]
    assert row["both_samples"] == 1
    assert row["mantel_fleiss_satisfied"] is False
    assert row["cmh_p_value"] is None and row["cmh_q_value"] is None
    assert row["cmh_statistic"] is not None
    assert report["tested_pair_count"] == 0
