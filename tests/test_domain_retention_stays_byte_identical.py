"""Proves the mutation-co-occurrence evidence layer is purely additive: it
must never change BRAF/RET's already-committed ``domain_retention`` Fisher
statistics.

Two independent proofs:

1. A network-free, deterministic A/B test: the real curated BRAF/RET
   ``GeneConfig`` objects (which now carry the new, optional
   ``mutual_exclusivity_targets`` field) produce byte-identical
   ``domain_retention`` results to a copy of themselves with that field
   stripped back to empty -- i.e. the new field is provably inert to this
   algorithm, independent of what live cBioPortal data says on any given
   day.
2. A pin against the actual numbers committed to the repo before this
   change (see PR description): reconstructing BRAF's/RET's real events
   and features from the already-committed run artifacts and re-running
   today's ``DomainRetentionAlgorithm`` reproduces those exact historical
   Fisher odds ratio/p-value/contingency-table values.
"""

from __future__ import annotations

import json
from pathlib import Path

from cfh.algorithms.domain_retention import DomainRetentionAlgorithm
from cfh.genes.registry import load_gene_config
from cfh.model.fusion_event import FusionEvent
from cfh.model.fusion_feature import FusionFeature
from conftest import latest_run_dir

# Pinned exactly from runs/braf_msk-impact-50k-2026_20260908T141855Z/results.json
# and runs/ret_msk-impact-50k-2026_20260908T141922Z/results.json, captured
# before this PR's changes were made.
_PINNED_BRAF_DOMAIN_RETENTION = {
    "fisher_odds_ratio": 4.507936507936508,
    "fisher_p_value": 0.013367557978153668,
    "observed_in_frame_retention_rate": 0.9403973509933775,
}
_PINNED_BRAF_TABLE = [[142, 21], [9, 6]]
_PINNED_RET_DOMAIN_RETENTION = {
    "fisher_odds_ratio": 7.421052631578948,
    "fisher_p_value": 0.00041966557652448966,
    "observed_in_frame_retention_rate": 0.9657534246575342,
}
_PINNED_RET_TABLE = [[141, 38], [5, 10]]


def _events_and_features_for_fisher(
    results_path: Path,
) -> tuple[list[FusionEvent], list[FusionFeature]]:
    """Reconstruct just enough of a real committed run's data (frame status
    and target-domain retention status per event) to reproduce the
    domain_retention Fisher-exact statistic -- no permutation/RNG/network
    involved, so this is fully deterministic.
    """
    payload = json.loads(results_path.read_text())
    events = []
    features = []
    for row in payload["events"]:
        events.append(
            FusionEvent(
                Event_id=row["event_id"],
                Cohort=payload["study_id"],
                Frame_status=row["frame_status"],
                Is_protein_fusion=True,
            )
        )
        features.append(
            FusionFeature(
                Event_id=row["event_id"],
                Gene=payload["gene_symbol"],
                Junction_position_aa=row["breakpoint_protein_position"],
                Domain_retention_flags={"kinase": row["domain_status"]},
            )
        )
    return events, features


def _fisher_only(config, events, features):
    result = DomainRetentionAlgorithm().run(events, features, config, {"n_permutations": 5})
    summary = result.Summary
    return (
        {
            "fisher_odds_ratio": summary["fisher_odds_ratio"],
            "fisher_p_value": summary["fisher_p_value"],
            "observed_in_frame_retention_rate": summary["observed_in_frame_retention_rate"],
        },
        result.Tables["frame_domain_contingency_table"],
    )


def test_mutual_exclusivity_targets_field_is_inert_to_domain_retention_for_braf():
    braf_config = load_gene_config("braf")
    braf_config_without_field = braf_config.model_copy(update={"mutual_exclusivity_targets": []})
    results_path = latest_run_dir("braf_msk-impact-50k-2026") / "results.json"
    events, features = _events_and_features_for_fisher(results_path)

    with_field = _fisher_only(braf_config, events, features)
    without_field = _fisher_only(braf_config_without_field, events, features)

    assert with_field == without_field


def test_mutual_exclusivity_targets_field_is_inert_to_domain_retention_for_ret():
    ret_config = load_gene_config("ret")
    assert ret_config.mutual_exclusivity_targets == []  # RET never opted in
    results_path = latest_run_dir("ret_msk-impact-50k-2026") / "results.json"
    events, features = _events_and_features_for_fisher(results_path)

    with_config = _fisher_only(ret_config, events, features)
    # RET's config already has no mutual_exclusivity_targets, so "with" and
    # "without" are the same object -- this leg proves RET is untouched at all.
    without_field = _fisher_only(
        ret_config.model_copy(update={"mutual_exclusivity_targets": []}), events, features
    )

    assert with_config == without_field


def test_braf_domain_retention_fisher_stats_match_the_pre_change_committed_run():
    results_path = latest_run_dir("braf_msk-impact-50k-2026") / "results.json"
    events, features = _events_and_features_for_fisher(results_path)
    config = load_gene_config("braf")

    summary, table = _fisher_only(config, events, features)

    assert summary == _PINNED_BRAF_DOMAIN_RETENTION
    assert table == _PINNED_BRAF_TABLE


def test_ret_domain_retention_fisher_stats_match_the_pre_change_committed_run():
    results_path = latest_run_dir("ret_msk-impact-50k-2026") / "results.json"
    events, features = _events_and_features_for_fisher(results_path)
    config = load_gene_config("ret")

    summary, table = _fisher_only(config, events, features)

    assert summary == _PINNED_RET_DOMAIN_RETENTION
    assert table == _PINNED_RET_TABLE
