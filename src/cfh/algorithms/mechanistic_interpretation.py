"""Counter-intuitive domain-status flagging: the capstone interpretation
algorithm for domain retention/disruption.

This is a pure aggregation algorithm, in the same spirit as
``composite_score``: it never re-runs a Fisher's-exact test or a
permutation test itself, and it never guesses at biology. It only answers
two gene-agnostic, purely statistical questions from the already-computed
``domain_retention``/``domain_disruption`` :class:`AlgorithmResult` objects
and the per-event domain-status calls the mapping layer already produced:

1. Is this gene's configured domain-retention (or domain-disruption)
   finding itself statistically supported (the same ``p < 0.05`` and
   ``odds ratio > 1`` gate :func:`cfh.reporting.domain_names.domain_interpretation_sentence`
   already uses to decide whether the report asserts anything at all)? If
   not, there is no statistically-established "expected" pattern for any
   event to be counter to, and this module says so rather than picking a
   side.

2. Among events that run counter to a *supported* pattern (e.g. domain
   *lost* despite a significant retention enrichment), is that minority
   recurrent enough -- the same partner gene appearing repeatedly among
   just the counter events -- to be statistically distinguishable from
   background noise, or is it a diffuse, non-recurrent scatter consistent
   with individual passenger events/classification noise at this cohort
   size?

Answering *why* a counter-intuitive event might still be oncogenic is
deliberately out of scope here: that is a human-curator judgment call
(``GeneConfig.mechanism_note``), never an algorithmically inferred guess.
This module's job ends at "is this pattern recurrent enough to be worth a
curator's attention," never "here is the mechanism."

No-ops (empty ``Summary``/``Tables``, no error) for a ``gene_pair`` config
-- the domain-retention/domain-disruption concept this module interprets
does not apply to promoter-swap/expression-driven pairs like TMPRSS2-ERG,
which rely on ``GeneConfig.mechanism_note`` alone (see
``real_benchmark._gene_pair_markdown_summary``).
"""

from __future__ import annotations

import math
from datetime import datetime, timezone
from typing import Any, Optional

from cfh.algorithms.base import Algorithm
from cfh.algorithms.registry import register
from cfh.genes.registry import GeneConfig, KeyDomain
from cfh.model.algorithm_result import AlgorithmResult
from cfh.model.fusion_event import FusionEvent
from cfh.model.fusion_feature import FusionFeature
from cfh.stats.breakpoint_tests import gene_breakpoint_domain_status_event_records

ALGORITHM_NAME = "mechanistic_interpretation"
ALGORITHM_VERSION = "1.0.0"

_ALPHA = 0.05
"""Same significance threshold :func:`cfh.reporting.domain_names.domain_interpretation_sentence`
uses -- kept identical so this module's "statistically supported" verdict
never contradicts that sentence's own gate."""

_RECURRENT_PARTNER_THRESHOLD = 2
"""A partner gene must appear at least this many times among the
counter-intuitive events specifically (not the whole cohort) before that
subgroup is called a candidate cluster rather than background noise. An
explicit, documented, small integer -- not tuned per gene -- exactly the
same kind of fixed, defensible threshold as ``composite_score``'s
``neg_log10_p_cap``."""

_RETENTION_CONTRADICTING_STATUSES = {"lost", "disrupted"}
_DISRUPTION_CONTRADICTING_STATUS = "retained"


def _as_algorithm_result(item: Any) -> AlgorithmResult:
    if isinstance(item, AlgorithmResult):
        return item
    return AlgorithmResult.model_validate(item)


def _results_by_name(algorithm_results: Any) -> dict[str, AlgorithmResult]:
    by_name: dict[str, AlgorithmResult] = {}
    for item in algorithm_results or []:
        result = _as_algorithm_result(item)
        by_name[result.Algorithm] = result
    return by_name


def _failed(result: Optional[AlgorithmResult]) -> bool:
    if result is None:
        return True
    return any(str(warning).startswith("Algorithm failed") for warning in (result.Warnings or []))


def _statistically_supported(
    result: Optional[AlgorithmResult],
) -> tuple[bool, float | None, float | None]:
    """``fisher_p_value < _ALPHA`` and ``fisher_odds_ratio > 1``.

    A positive-infinite odds ratio (a "clean separation" 2x2 table, e.g. a
    gene where every in-frame fusion retains the domain and no other event
    does -- see ``fishers_frame_domain_test``'s own docstring on this being
    an expected, not a numerical-bug, outcome) is treated as supported: it
    is the *most*, not least, significant possible result, so only
    ``p_value``'s finiteness is required, never ``odds_ratio``'s."""
    if result is None or _failed(result):
        return False, None, None
    summary = result.Summary or {}
    p_value = summary.get("fisher_p_value")
    odds_ratio = summary.get("fisher_odds_ratio")
    p_value = (
        p_value if isinstance(p_value, (int, float)) and not isinstance(p_value, bool) else None
    )
    odds_ratio = (
        odds_ratio
        if isinstance(odds_ratio, (int, float)) and not isinstance(odds_ratio, bool)
        else None
    )
    if (
        p_value is None
        or odds_ratio is None
        or not math.isfinite(p_value)
        or math.isnan(odds_ratio)
    ):
        return False, p_value, odds_ratio
    return (p_value < _ALPHA and odds_ratio > 1), p_value, odds_ratio


def _configured_domain_names(domains: list[KeyDomain]) -> list[str]:
    """De-duplicated, order-preserving domain display names -- a small,
    algorithms-layer-local copy of
    :func:`cfh.reporting.domain_names.configured_domain_names` (never
    imported from here: the algorithms layer stays independent of the
    reporting/presentation layer, same as every other registered
    algorithm)."""
    names: list[str] = []
    for domain in domains:
        if domain.name and domain.name != domain.accession and domain.name not in names:
            names.append(domain.name)
    return names


def _is_in_frame_protein_fusion(event: FusionEvent) -> bool:
    """Same definition ``build_frame_domain_contingency_table`` uses to
    pick the "in-frame" column. The retention/disruption claim being
    interpreted here is specifically about in-frame fusions being enriched
    for a domain status -- an out-of-frame event's domain status was never
    part of that claim's population, so it can never be "counter" to it."""
    return event.Is_protein_fusion is True and event.Frame_status == "in-frame"


def _partner_gene(event: FusionEvent, target_gene: str | None) -> str:
    target = (target_gene or "").upper()
    if str(event.Site1_gene or "").upper() == target:
        return event.Site2_gene or "unknown"
    if str(event.Site2_gene or "").upper() == target:
        return event.Site1_gene or "unknown"
    return "unknown"


def _analyze_effect(
    *,
    effect: str,
    domains: list[KeyDomain],
    algorithm_result: Optional[AlgorithmResult],
    events: list[FusionEvent],
    features: list[FusionFeature],
    gene_config: GeneConfig,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Compute one effect's (``"retention"`` or ``"disruption"``) summary
    fields and detail tables. Returns ``(summary_fields, tables)``."""
    prefix = effect
    if not domains:
        return (
            {
                f"{prefix}_domains_configured": False,
                f"{prefix}_domain_names": [],
                f"{prefix}_statistically_supported": False,
                f"{prefix}_fisher_p_value": None,
                f"{prefix}_fisher_odds_ratio": None,
                f"{prefix}_counter_intuitive_confidence": "not_applicable",
                f"{prefix}_counter_intuitive_count": 0,
                f"{prefix}_counter_intuitive_total": 0,
                f"{prefix}_counter_intuitive_percent": None,
            },
            {},
        )

    supported, p_value, odds_ratio = _statistically_supported(algorithm_result)
    domain_names = _configured_domain_names(domains)
    summary_fields: dict[str, Any] = {
        f"{prefix}_domains_configured": True,
        f"{prefix}_domain_names": domain_names,
        f"{prefix}_statistically_supported": supported,
        f"{prefix}_fisher_p_value": p_value,
        f"{prefix}_fisher_odds_ratio": odds_ratio,
    }

    if not supported:
        # Too weak to establish an "expected" pattern at all -- there is
        # nothing for any event to be counter *to*, so no per-event
        # flagging is attempted. This is the "too weak to associate
        # causation with" branch of the analysis.
        summary_fields[f"{prefix}_counter_intuitive_confidence"] = "not_applicable"
        summary_fields[f"{prefix}_counter_intuitive_count"] = 0
        summary_fields[f"{prefix}_counter_intuitive_total"] = 0
        summary_fields[f"{prefix}_counter_intuitive_percent"] = None
        return summary_fields, {}

    contradicting = (
        _RETENTION_CONTRADICTING_STATUSES
        if effect == "retention"
        else {_DISRUPTION_CONTRADICTING_STATUS}
    )
    records = gene_breakpoint_domain_status_event_records(
        events, features, gene_config, domains=domains
    )
    event_by_id = {event.Event_id: event for event in events}
    in_frame_records = []
    for record in records:
        record_event = event_by_id.get(record[0])
        if record_event is not None and _is_in_frame_protein_fusion(record_event):
            in_frame_records.append(record)
    total = len(in_frame_records)
    counter_events: list[dict[str, Any]] = []
    for event_id, breakpoint_position, status in in_frame_records:
        if status not in contradicting:
            continue
        event = event_by_id.get(event_id)
        counter_events.append(
            {
                "event_id": event_id,
                "sample_id": event.Sample_id if event else None,
                "partner_gene": _partner_gene(event, gene_config.gene_symbol)
                if event
                else "unknown",
                "breakpoint_protein_position": breakpoint_position,
                "domain_status": status,
            }
        )

    count = len(counter_events)
    partner_counts: dict[str, int] = {}
    for row in counter_events:
        partner_counts[row["partner_gene"]] = partner_counts.get(row["partner_gene"], 0) + 1
    recurrent_pairs = sorted(
        (
            (partner, partner_count)
            for partner, partner_count in partner_counts.items()
            if partner_count >= _RECURRENT_PARTNER_THRESHOLD
        ),
        key=lambda pair: (-pair[1], pair[0]),
    )
    recurrent_partners = [
        {"partner_gene": partner, "count": partner_count}
        for partner, partner_count in recurrent_pairs
    ]

    if count == 0:
        confidence = "none"
    elif recurrent_partners:
        confidence = "possible_subcluster"
    else:
        confidence = "insufficient_recurrence"

    summary_fields[f"{prefix}_counter_intuitive_confidence"] = confidence
    summary_fields[f"{prefix}_counter_intuitive_count"] = count
    summary_fields[f"{prefix}_counter_intuitive_total"] = total
    summary_fields[f"{prefix}_counter_intuitive_percent"] = (
        (100.0 * count / total) if total else None
    )
    summary_fields[f"{prefix}_recurrent_partner_threshold"] = _RECURRENT_PARTNER_THRESHOLD

    tables: dict[str, Any] = {}
    if counter_events:
        tables[f"{prefix}_counter_intuitive_events"] = counter_events
    if recurrent_partners:
        tables[f"{prefix}_counter_intuitive_recurrent_partners"] = recurrent_partners
    return summary_fields, tables


@register(ALGORITHM_NAME)
class MechanisticInterpretationAlgorithm(Algorithm):
    """Flag domain-status events that run counter to this gene's own
    statistically-supported retention/disruption pattern, and judge purely
    from partner-gene recurrence whether that minority looks like a
    candidate subcluster or background noise. See the module docstring for
    the full rationale and scope boundary."""

    VERSION = ALGORITHM_VERSION

    DEPENDS_ON = ("domain_retention", "domain_disruption")

    def run(
        self,
        events: list[FusionEvent],
        features: list[FusionFeature],
        gene_config: GeneConfig,
        params: dict,
    ) -> AlgorithmResult:
        if gene_config.gene_pair is not None:
            return AlgorithmResult(
                Algorithm=ALGORITHM_NAME,
                Algorithm_version=ALGORITHM_VERSION,
                Parameters={},
                Summary={},
                Tables={},
                Warnings=[
                    "mechanistic_interpretation does not apply to gene_pair configs "
                    "(no domain-retention concept to interpret); see "
                    "GeneConfig.mechanism_note for pair-level curated mechanism text."
                ],
                Created_at=datetime.now(timezone.utc),
            )

        params = params or {}
        results_by_name = _results_by_name(params.get("algorithm_results"))

        retention_summary, retention_tables = _analyze_effect(
            effect="retention",
            domains=gene_config.key_domains,
            algorithm_result=results_by_name.get("domain_retention"),
            events=events,
            features=features,
            gene_config=gene_config,
        )
        disruption_summary, disruption_tables = _analyze_effect(
            effect="disruption",
            domains=gene_config.disruption_required_domains,
            algorithm_result=results_by_name.get("domain_disruption"),
            events=events,
            features=features,
            gene_config=gene_config,
        )

        return AlgorithmResult(
            Algorithm=ALGORITHM_NAME,
            Algorithm_version=ALGORITHM_VERSION,
            Parameters={
                "alpha": _ALPHA,
                "recurrent_partner_threshold": _RECURRENT_PARTNER_THRESHOLD,
            },
            Summary={**retention_summary, **disruption_summary},
            Tables={**retention_tables, **disruption_tables},
            Created_at=datetime.now(timezone.utc),
        )


__all__ = ["MechanisticInterpretationAlgorithm", "ALGORITHM_NAME", "ALGORITHM_VERSION"]
