"""Tumor-type-stratified CMH for the frame x domain-status table."""

from __future__ import annotations

import json

from cfh.algorithms.domain_retention import DomainRetentionAlgorithm
from cfh.gene_comparison import collect_p_values_from_algorithm_results
from cfh.genes.registry import load_gene_config
from cfh.model.fusion_event import FusionEvent
from cfh.model.fusion_feature import FusionFeature
from cfh.stats.breakpoint_tests import build_frame_domain_contingency_table
from cfh.stats.stratification import (
    UNKNOWN_STRATUM,
    stratified_cmh_summary,
    stratified_frame_domain_tables,
)

CONFIG = load_gene_config("BRAF")
KEY = CONFIG.key_domains[0].key or CONFIG.key_domains[0].name


def _cohort(counts: dict[str | None, list[list[int]]]):
    """Build events whose per-stratum frame x retention tables equal ``counts``."""
    events, features = [], []
    for stratum, table in counts.items():
        for row, status in enumerate(("retained", "lost")):
            for column, frame in enumerate(("in-frame", "out-of-frame")):
                for _ in range(table[row][column]):
                    event_id = f"E{len(events)}"
                    events.append(
                        FusionEvent(
                            Event_id=event_id,
                            Cohort="c",
                            Oncotree_code=stratum,
                            Is_protein_fusion=True,
                            Frame_status=frame,
                        )
                    )
                    features.append(
                        FusionFeature(
                            Event_id=event_id,
                            Gene="BRAF",
                            Junction_position_aa=400 + len(events),
                            Domain_retention_flags={KEY: status},
                        )
                    )
    return events, features


# Within each tumor type retention is unrelated to frame (OR = 1), but the
# pooled table shows a strong association because tumor types differ in both
# their in-frame rate and their retention rate.
SIMPSON = {"LUAD": [[40, 10], [4, 1]], "SKCM": [[1, 4], [10, 40]]}


def test_strata_sum_to_pooled_table_and_unknown_is_labelled():
    events, features = _cohort({**SIMPSON, None: [[1, 0], [0, 1]]})
    tables = stratified_frame_domain_tables(events, features, CONFIG)
    pooled = build_frame_domain_contingency_table(events, features, CONFIG)
    assert set(tables) == {"LUAD", "SKCM", UNKNOWN_STRATUM}
    assert [[sum(t[r][c] for t in tables.values()) for c in (0, 1)] for r in (0, 1)] == pooled


def test_simpsons_paradox_pooled_association_vanishes_after_stratifying():
    events, features = _cohort(SIMPSON)
    result = DomainRetentionAlgorithm().run(events, features, CONFIG, {"n_permutations": 5})
    assert result.Summary["fisher_p_value"] < 1e-6
    stratified = result.Summary["tumor_type_stratified"]
    assert stratified["mh_common_odds_ratio"] == 1.0
    assert stratified["direction"] == "none"
    assert stratified["cmh_p_value"] > 0.9
    assert stratified["informative_strata"] == 2
    strata = result.Tables["frame_domain_contingency_tables_by_tumor_type"]
    assert {row["stratum"] for row in strata} == {"LUAD", "SKCM"}


def test_stratified_p_value_is_not_in_the_fdr_family():
    events, features = _cohort(SIMPSON)
    result = DomainRetentionAlgorithm().run(events, features, CONFIG, {"n_permutations": 5})
    rows = collect_p_values_from_algorithm_results("BRAF", "s", [result], source="t")
    assert {row["test"] for row in rows} == {"fisher", "permutation"}


def test_summary_is_strict_json_for_complete_separation_and_empty_input():
    infinite = stratified_cmh_summary({"A": [[5, 0], [0, 5]]})
    assert infinite["mh_common_odds_ratio"] == "infinity"
    json.dumps(infinite, allow_nan=False)
    empty = stratified_cmh_summary({"A": [[0, 0], [0, 0]]})
    assert empty["cmh_p_value"] is None and empty["total_strata"] == 0
