"""Composite prioritization heuristic for fusion partners.

This is the capstone aggregation algorithm: it consumes the already-computed
:class:`~cfh.model.algorithm_result.AlgorithmResult` objects produced by the
other registered algorithms (run via :func:`cfh.orchestrator.run.run_algorithms`)
and combines them into a single, interpretable, per-fusion-partner ranked
    prioritization score. It is not a calibrated probability of functional
    relevance. It never re-runs a Fisher's-exact
test, a permutation test, or an MLE/CI computation itself -- those numbers
are read straight out of the upstream results' ``Summary`` blocks. The only
new computation this module performs is (a) a documented normalization of
each upstream number onto a common ``[0, 1]`` scale, (b) joining each
partner gene's own breakpoints against an already-inferred cutpoint to
measure per-partner proximity, and (c) a weighted average of whichever
sub-scores are actually available for this gene/run.

Sub-scores and how they are combined
-------------------------------------

One row is produced per fusion-partner gene, at the same partner-gene grain
the ``frequency`` algorithm already reports at. Up to five sub-scores feed
each row's composite score, each independently scaled to ``[0, 1]``:

``recurrence`` (always present)
    ``partner_event_count / total_event_count`` from the ``frequency``
    result's ``Partner_gene_counts`` table -- the fraction of all analyzed
    events attributable to that partner. Varies per partner.

``domain_retention`` (present whenever a ``domain_retention`` result with a
p-value is supplied)
    A gene-level (not partner-specific) evidence score derived from the
    ``domain_retention`` result: ``min(-log10(p), CAP) / CAP``, using the
    permutation empirical p-value when available (falling back to the
    Fisher exact p-value). ``CAP`` defaults to 10 (i.e. p=1e-10 or smaller
    saturates the score at 1.0; p=1 scores 0.0). Applied identically to
    every partner row because the underlying test is not partner-specific
    -- it is evidence about the gene's fusions as a whole.

``domain_disruption`` (present only when the gene configures
``disruption_required_domains`` AND the algorithm produced a non-null
p-value)
    Same ``[0, 1]`` transform as ``domain_retention``, applied to the
    ``domain_disruption`` result. A gene that leaves ``disruption_required_domains``
    unset gets a graceful no-op ``domain_disruption`` result with ``None``
    statistics; this sub-score is then simply excluded from the weighted
    average for that gene -- never treated as zero evidence.

``cutpoint_proximity`` (present only when ``cutpoint_detection`` produced a
    determinable cutpoint with corrected p < 0.05, and only for partners with at least one mapped
breakpoint)
    The one genuinely partner-varying statistical sub-score. For each
    partner, the mean absolute distance (in amino acids) between that
    partner's mapped breakpoints and the already-inferred cutpoint
    (``Summary["inferred_cutpoint_aa"]``) is computed directly from the raw
    ``events``/``features`` inputs (a join, not a new statistical test), then
    scaled so the partner closest to the cutpoint scores near 1.0 and the
    partner furthest away scores near 0.0: ``1 - distance / max_distance``
    across partners with a mapped distance. A partner with no mapped
    breakpoint simply has this sub-score excluded, not zeroed.

``confidence_certainty`` (present only when ``confidence_stats`` produced an
MLE/CI block)
    Gene-level, from ``confidence_stats``'s ``Summary["mle"]["groups"]``:
    the mean confidence-interval width across groups, inverted so a
    narrower (more certain) interval scores higher: ``1 - mean(ci_high -
    ci_low)``. Applied identically to every partner row, for the same
    reason as ``domain_retention``.

Composite score
    For each partner row, ``Composite_score`` is the weighted average of
    whichever sub-scores are available for that row:
    ``sum(weight[s] * value[s] for s in available) / sum(weight[s] for s in
    available)``. A missing/not-applicable sub-score is dropped from BOTH
    the numerator and the denominator, so it contributes neither positive
    nor negative evidence -- it is excluded, not scored as zero. Default
    weights (overridable via ``params["weights"]``, summing to 1.0):
    ``recurrence=0.30``, ``domain_retention=0.25``, ``domain_disruption=0.20``,
    ``cutpoint_proximity=0.15``, ``confidence_certainty=0.10``. Overridden
    weights must be finite, non-negative, not all zero, and use only the
    five names above -- ``run`` raises ``ValueError`` otherwise, since an
    invalid weight set would silently break the ``[0, 1]`` guarantee or
    produce a fabricated score. In the residual edge case where a row's own
    applicable sub-scores happen to all carry zero weight even though the
    overall weight set is valid, ``Composite_score`` is ``None`` for that
    row (with a corresponding warning) rather than a fabricated ``0.0`` --
    it sorts after every row with a real score.

The output table is sorted by ``Composite_score`` descending (ties broken by
event count, then partner-gene name, for determinism) and each row carries
its full component breakdown plus a ``Components_applicable`` list, so the
reporting/CLI layer can render both the ranking and exactly how it was
computed -- nothing here is a black box.
"""

from __future__ import annotations

import math
from datetime import datetime, timezone
from typing import Any, Optional

from cfh.algorithms.base import Algorithm
from cfh.algorithms.frequency import _partner_gene
from cfh.algorithms.registry import register
from cfh.genes.registry import GeneConfig
from cfh.model.algorithm_result import AlgorithmResult
from cfh.model.fusion_event import FusionEvent
from cfh.model.fusion_feature import FusionFeature

ALGORITHM_NAME = "composite_score"
ALGORITHM_VERSION = "0.2.0"

DEFAULT_WEIGHTS: dict[str, float] = {
    "recurrence": 0.30,
    "domain_retention": 0.25,
    "domain_disruption": 0.20,
    "cutpoint_proximity": 0.15,
    "confidence_certainty": 0.10,
}

DEFAULT_NEG_LOG10_P_CAP = 10.0
CUTPOINT_CORRECTED_P_ALPHA = 0.05


def _validate_weight_keys(overrides: dict[str, Any]) -> None:
    """Every override key must name a real sub-score.

    A typo'd key (e.g. ``"recurence"``) would otherwise be silently
    ignored -- the intended weight never applies and the default for the
    correctly-spelled key stays in effect with no indication anything was
    wrong -- so it is rejected up front instead.
    """
    unknown = sorted(set(overrides) - set(DEFAULT_WEIGHTS))
    if unknown:
        raise ValueError(
            f"unrecognized weights override key(s) {unknown}; expected one of "
            f"{sorted(DEFAULT_WEIGHTS)}"
        )


def _validate_weights(weights: dict[str, float]) -> None:
    """Every weight must be a finite, non-negative number, and at least one
    must be positive.

    A negative or non-finite weight would break the documented ``[0, 1]``
    composite-score guarantee (the weighted average is only bounded when
    every weight and every sub-score is non-negative). A non-numeric value
    is checked explicitly first so it raises the documented ``ValueError``
    rather than a ``TypeError`` from the arithmetic below. All-zero weights
    are rejected outright: every row would silently compute a fabricated
    ``0.0`` composite score (0/0 guarded to 0.0) that is visually
    indistinguishable from a genuine "no evidence" score.
    """
    for name, weight in weights.items():
        is_numeric = isinstance(weight, (int, float)) and not isinstance(weight, bool)
        if not is_numeric or not math.isfinite(weight) or weight < 0:
            raise ValueError(
                f"weights[{name!r}] must be a finite, non-negative number; got {weight!r}"
            )
    if all(weight == 0 for weight in weights.values()):
        raise ValueError(
            "all composite_score weights are zero; at least one weight must be positive"
        )


def _validate_cap(cap: float) -> None:
    is_numeric = isinstance(cap, (int, float)) and not isinstance(cap, bool)
    if not is_numeric or not math.isfinite(cap) or cap <= 0:
        raise ValueError(f"neg_log10_p_cap must be a finite, positive number; got {cap!r}")


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
    """Treat an orchestrator-wrapped exception as an unavailable result."""
    if result is None:
        return True
    return any(str(warning).startswith("Algorithm failed") for warning in (result.Warnings or []))


def _clipped_neg_log10(p_value: float, cap: float) -> float:
    if p_value <= 0:
        return 1.0
    return max(0.0, min(1.0, -math.log10(p_value) / cap))


def _domain_test_score(
    result: Optional[AlgorithmResult], cap: float
) -> tuple[Optional[float], Optional[float]]:
    """Return ``(score, p_value)`` for a domain_retention/domain_disruption result.

    ``None`` for both when the result is absent, failed, or was a graceful
    no-op (unconfigured domains) with no p-value to report.
    """
    if result is None or _failed(result):
        return None, None
    summary = result.Summary or {}
    p_value = summary.get("permutation_empirical_p_value")
    if p_value is None:
        p_value = summary.get("fisher_p_value")
    if p_value is None:
        return None, None
    return _clipped_neg_log10(p_value, cap), p_value


def _confidence_certainty_score(result: Optional[AlgorithmResult]) -> Optional[float]:
    if result is None or _failed(result):
        return None
    mle = (result.Summary or {}).get("mle")
    if not mle:
        return None
    widths = [
        group["ci_high"] - group["ci_low"]
        for group in (mle.get("groups") or {}).values()
        if group.get("n", 0) > 0
    ]
    if not widths:
        return None
    mean_width = sum(widths) / len(widths)
    return max(0.0, min(1.0, 1.0 - mean_width))


def _cutpoint_proximity_scores(
    result: Optional[AlgorithmResult],
    events: list[FusionEvent],
    features: list[FusionFeature],
    gene_symbol: Optional[str],
) -> tuple[dict[str, float], Optional[int]]:
    """Return ``(proximity_score_by_partner, inferred_cutpoint_aa)``.

    Reads the already-inferred cutpoint position from ``result`` and joins
    it against the raw per-event breakpoint positions to measure each
    partner's mean distance -- no statistic is recomputed. An empty dict is
    returned (rather than zeros) when no determinable cutpoint or no mapped
    breakpoints are available.
    """
    if result is None or _failed(result) or not gene_symbol:
        return {}, None
    summary = result.Summary or {}
    if not summary.get("determinable"):
        return {}, None
    corrected_p = summary.get("corrected_p_value")
    if (
        isinstance(corrected_p, bool)
        or not isinstance(corrected_p, (int, float))
        or not math.isfinite(corrected_p)
        or not 0 <= corrected_p < CUTPOINT_CORRECTED_P_ALPHA
    ):
        return {}, None
    cutpoint = summary.get("inferred_cutpoint_aa")
    if cutpoint is None:
        return {}, None

    event_by_id = {event.Event_id: event for event in events}
    distances_by_partner: dict[str, list[int]] = {}
    for feature in features:
        if feature.Gene.upper() != gene_symbol.upper() or feature.Junction_position_aa is None:
            continue
        event = event_by_id.get(feature.Event_id)
        if event is None:
            continue
        partner = _partner_gene(event, gene_symbol)
        distances_by_partner.setdefault(partner, []).append(
            abs(feature.Junction_position_aa - cutpoint)
        )

    if not distances_by_partner:
        return {}, cutpoint

    mean_distance_by_partner = {
        partner: sum(distances) / len(distances)
        for partner, distances in distances_by_partner.items()
    }
    max_distance = max(mean_distance_by_partner.values())
    if max_distance == 0:
        return {partner: 1.0 for partner in mean_distance_by_partner}, cutpoint
    return {
        partner: max(0.0, 1.0 - (distance / max_distance))
        for partner, distance in mean_distance_by_partner.items()
    }, cutpoint


@register(ALGORITHM_NAME)
class CompositeScoreAlgorithm(Algorithm):
    """Aggregate other algorithms' already-computed results into one ranked
    prioritization score per fusion partner, not a significance test.

    See the module docstring for the full sub-score/combination formula.

    Expected ``params`` keys:
        algorithm_results (list, required): the ``AlgorithmResult`` objects
            (or their ``.model_dump()`` equivalents) already produced by
            running ``frequency`` and, as available, ``domain_retention``,
            ``domain_disruption``, ``cutpoint_detection``, and
            ``confidence_stats`` via the orchestrator against the same
            ``events``/``features``/``gene_config``. ``frequency`` is the
            one required entry -- everything else is optional and its
            corresponding sub-score is gracefully excluded when absent,
            failed, or a not-applicable no-op.
        weights (dict, optional): override any of ``DEFAULT_WEIGHTS``. Every
            key must be one of the five recognized sub-score names, every
            value must be a finite, non-negative number, and at least one
            weight must be positive, or ``run`` raises ``ValueError`` --
            an unrecognized key, a negative/non-finite/non-numeric weight,
            or an all-zero weight set would each break the documented
            ``[0, 1]`` composite-score guarantee or silently produce a
            fabricated score.
        neg_log10_p_cap (float, optional): saturation cap used to scale
            p-value-derived sub-scores onto ``[0, 1]``, default 10.0. Must
            be finite and positive.
    """

    VERSION = ALGORITHM_VERSION

    DEPENDS_ON = (
        "frequency",
        "domain_retention",
        "domain_disruption",
        "cutpoint_detection",
        "confidence_stats",
    )

    def run(
        self,
        events: list[FusionEvent],
        features: list[FusionFeature],
        gene_config: GeneConfig,
        params: dict,
    ) -> AlgorithmResult:
        params = params or {}
        weight_overrides = params.get("weights") or {}
        _validate_weight_keys(weight_overrides)
        weights = {**DEFAULT_WEIGHTS, **weight_overrides}
        cap = params.get("neg_log10_p_cap", DEFAULT_NEG_LOG10_P_CAP)
        _validate_weights(weights)
        _validate_cap(cap)

        results_by_name = _results_by_name(params.get("algorithm_results"))
        frequency_result = results_by_name.get("frequency")
        if frequency_result is None or _failed(frequency_result):
            raise ValueError(
                "composite_score requires a completed 'frequency' AlgorithmResult in "
                "params['algorithm_results']; recurrence is the required baseline input "
                "and every other sub-score is optional."
            )
        partner_counts = (frequency_result.Tables or {}).get("Partner_gene_counts") or []
        total_events = sum(row["Event_count"] for row in partner_counts)

        domain_retention_score, domain_retention_p = _domain_test_score(
            results_by_name.get("domain_retention"), cap
        )
        domain_disruption_result = results_by_name.get("domain_disruption")
        domain_disruption_score, domain_disruption_p = _domain_test_score(
            domain_disruption_result, cap
        )
        cutpoint_scores_by_partner, inferred_cutpoint = _cutpoint_proximity_scores(
            results_by_name.get("cutpoint_detection"),
            events,
            features,
            gene_config.gene_symbol if gene_config else None,
        )
        confidence_score = _confidence_certainty_score(results_by_name.get("confidence_stats"))

        rows: list[dict] = []
        for entry in partner_counts:
            partner = entry["Partner_gene"]
            count = entry["Event_count"]
            recurrence_score = count / total_events if total_events else 0.0
            partner_cutpoint_score = cutpoint_scores_by_partner.get(partner)

            components: dict[str, float] = {"recurrence": recurrence_score}
            if domain_retention_score is not None:
                components["domain_retention"] = domain_retention_score
            if domain_disruption_score is not None:
                components["domain_disruption"] = domain_disruption_score
            if partner_cutpoint_score is not None:
                components["cutpoint_proximity"] = partner_cutpoint_score
            if confidence_score is not None:
                components["confidence_certainty"] = confidence_score

            weight_sum = sum(weights[name] for name in components)
            composite = (
                sum(weights[name] * value for name, value in components.items()) / weight_sum
                if weight_sum > 0
                else None
            )
            rows.append(
                {
                    "Partner_gene": partner,
                    "Event_count": count,
                    "Recurrence_score": recurrence_score,
                    "Domain_retention_score": domain_retention_score,
                    "Domain_disruption_score": domain_disruption_score,
                    "Cutpoint_proximity_score": partner_cutpoint_score,
                    "Confidence_certainty_score": confidence_score,
                    "Components_applicable": sorted(components),
                    "Composite_score": composite,
                }
            )

        rows.sort(
            key=lambda row: (
                row["Composite_score"] is None,
                -(row["Composite_score"] if row["Composite_score"] is not None else 0.0),
                -row["Event_count"],
                row["Partner_gene"],
            )
        )
        for rank, row in enumerate(rows, start=1):
            row["Rank"] = rank

        warnings: list[str] = []
        if any(row["Composite_score"] is None for row in rows):
            warnings.append(
                "One or more partner rows had zero total weight among their applicable "
                "sub-scores; Composite_score is None for those rows rather than a "
                "fabricated 0.0."
            )
        if domain_retention_score is None:
            warnings.append(
                "domain_retention sub-score excluded: no applicable domain_retention "
                "result was supplied."
            )
        if domain_disruption_result is None:
            warnings.append(
                "domain_disruption sub-score excluded: no domain_disruption result was "
                "supplied (gene may not configure disruption_required_domains)."
            )
        elif domain_disruption_score is None:
            warnings.append(
                "domain_disruption sub-score excluded: the supplied result had no "
                "applicable p-value (e.g. disruption_required_domains not configured "
                "for this gene)."
            )
        if not cutpoint_scores_by_partner:
            warnings.append(
                "cutpoint_proximity sub-score excluded: corrected cutpoint p-value was "
                "missing or >= 0.05, cutpoint was not determinable, or mapped "
                "breakpoints were unavailable."
            )
        if confidence_score is None:
            warnings.append(
                "confidence_certainty sub-score excluded: no applicable confidence_stats "
                "MLE/CI result was supplied."
            )

        return AlgorithmResult(
            Algorithm=ALGORITHM_NAME,
            Algorithm_version=ALGORITHM_VERSION,
            Parameters={
                "weights": weights,
                "neg_log10_p_cap": cap,
                "cutpoint_corrected_p_alpha": CUTPOINT_CORRECTED_P_ALPHA,
            },
            Summary={
                "score_interpretation": "prioritization_heuristic_not_probability_or_significance",
                "gene_symbol": gene_config.gene_symbol if gene_config else None,
                "n_partners_ranked": len(rows),
                "total_events": total_events,
                "components_applicable": {
                    "recurrence": True,
                    "domain_retention": domain_retention_score is not None,
                    "domain_disruption": domain_disruption_score is not None,
                    "cutpoint_proximity": bool(cutpoint_scores_by_partner),
                    "confidence_certainty": confidence_score is not None,
                },
                "domain_retention_p_value_used": domain_retention_p,
                "domain_disruption_p_value_used": domain_disruption_p,
                "inferred_cutpoint_aa": inferred_cutpoint,
            },
            Tables={"composite_evidence_ranking": rows},
            Warnings=warnings,
            Created_at=datetime.now(timezone.utc),
        )
