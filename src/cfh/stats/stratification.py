"""Tumor-type-stratified versions of the frame x domain-status 2x2 test.

The pooled Fisher test in :mod:`cfh.stats.breakpoint_tests` mixes tumor
types. In a pan-cancer cohort, partner genes, frame calls, and domain status
all vary by tissue, so a pooled association can reflect tumor-type
composition rather than a within-tissue relationship. This module builds one
2x2 table per tumor type and summarizes them with the Cochran-Mantel-Haenszel
test and the Mantel-Haenszel common odds ratio.

The CMH statistic is a two-sided chi-square(1) test of conditional
independence; the pooled Fisher test is one-sided. ``direction`` reports which
way the pooled-over-strata association points. Strata whose table has a zero
margin contribute no information and are counted, not dropped silently.
These results are reported alongside the pooled test; they are not part of the
cross-gene FDR family.
"""

from __future__ import annotations

import math
from typing import Any

from cfh.genes.registry import GeneConfig, KeyDomain
from cfh.model.fusion_event import FusionEvent
from cfh.model.fusion_feature import FusionFeature
from cfh.stats.breakpoint_tests import RETAINED_STATUSES, build_frame_domain_contingency_table
from cfh.stats.cochran_mantel_haenszel import cochran_mantel_haenszel

UNKNOWN_STRATUM = "UNKNOWN"
STRATIFICATION_FIELD = "oncotree_code"


def tumor_type_stratum(event: FusionEvent) -> str:
    """OncoTree code, falling back to the free-text tumor type, then UNKNOWN."""
    return event.Oncotree_code or event.Tumor_type or UNKNOWN_STRATUM


def stratified_cmh_summary(tables_by_stratum: dict[str, list[list[int]]]) -> dict[str, Any]:
    """Summarize per-stratum 2x2 tables with CMH and the MH common odds ratio."""
    nonempty = {
        stratum: table
        for stratum, table in tables_by_stratum.items()
        if sum(sum(row) for row in table)
    }
    if not nonempty:
        return {
            "stratified_by": STRATIFICATION_FIELD,
            "cmh_statistic": None,
            "cmh_p_value": None,
            "mh_common_odds_ratio": None,
            "direction": None,
            "informative_strata": 0,
            "total_strata": 0,
            "unknown_stratum_observations": 0,
        }
    result = cochran_mantel_haenszel(list(nonempty.values()))
    odds_ratio = result.common_odds_ratio
    if odds_ratio is None:
        direction = None
    elif odds_ratio == 1:
        direction = "none"
    else:
        direction = "positive" if odds_ratio > 1 else "negative"
    unknown = nonempty.get(UNKNOWN_STRATUM)
    return {
        "stratified_by": STRATIFICATION_FIELD,
        "cmh_statistic": result.statistic,
        "cmh_p_value": result.p_value,
        # Strict JSON: an infinite MH odds ratio is reported as a string.
        "mh_common_odds_ratio": (
            "infinity" if odds_ratio is not None and math.isinf(odds_ratio) else odds_ratio
        ),
        "direction": direction,
        "informative_strata": result.informative_strata,
        "total_strata": result.total_strata,
        "unknown_stratum_observations": sum(sum(row) for row in unknown) if unknown else 0,
    }


def stratified_frame_domain_tables(
    events: list[FusionEvent],
    features: list[FusionFeature],
    gene_config: GeneConfig,
    *,
    domains: list[KeyDomain] | None = None,
    hit_statuses: frozenset[str] | None = None,
) -> dict[str, list[list[int]]]:
    """Per-tumor-type version of ``build_frame_domain_contingency_table``.

    Uses exactly the same inclusion and orientation rules, so the element-wise
    sum of the returned tables equals the pooled table.
    """
    event_by_id = {event.Event_id: event for event in events}
    events_by_stratum: dict[str, list[FusionEvent]] = {}
    for event in events:
        events_by_stratum.setdefault(tumor_type_stratum(event), []).append(event)
    tables: dict[str, list[list[int]]] = {}
    for stratum, stratum_events in sorted(events_by_stratum.items()):
        stratum_ids = {event.Event_id for event in stratum_events}
        stratum_features = [
            feature
            for feature in features
            if feature.Event_id in stratum_ids and feature.Event_id in event_by_id
        ]
        tables[stratum] = build_frame_domain_contingency_table(
            stratum_events,
            stratum_features,
            gene_config,
            domains=domains,
            hit_statuses=hit_statuses or RETAINED_STATUSES,
        )
    return tables
