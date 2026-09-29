"""Mutation/CNA co-occurrence vs mutual-exclusivity evidence layer.

Tests whether this gene's own fusion-positive sample set co-occurs with, or
is mutually exclusive with, a configured comparator gene+alteration-type
across the cohort (e.g. BRAF fusions vs BRAF V600E point mutations -- both
MAPK-activating, hypothesized to be near-mutually-exclusive). A two-sided
Fisher's exact test on the 2x2 contingency table (fusion status x
comparator-alteration status, over every profiled cohort sample -- see
:mod:`cfh.stats.cooccurrence_tests`) is the pre-specified statistic; the
observed odds ratio's direction (>1 co-occurring, <1 mutually exclusive) is
reported alongside the p-value rather than assumed.

Opt-in via ``GeneConfig.mutual_exclusivity_targets``: a gene that leaves
this unset (its default, empty list) produces a no-op result with no
statistics computed, the same graceful-skip pattern already used by
``domain_disruption``/``exon_retention``/``joint_partner`` for their
respective optional fields.

Even when configured, this algorithm needs data that isn't part of the
gene's own fusion ``events``/``features``: the cohort's full sample
universe, and the comparator's per-sample alteration calls. These are
supplied via ``params["cohort_sample_ids"]``/``params["comparator_alterations"]``
(populated by a live-fetch caller such as
:func:`cfh.real_benchmark.run_real_benchmark`). A gene that IS configured
but is run without that data (e.g. the default orchestrator sweep with no
live-fetch step) still no-ops gracefully, with a warning distinct from the
"not configured" one, never crashing.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from typing import Any

from cfh.algorithms.base import Algorithm
from cfh.algorithms.registry import register
from cfh.genes.registry import GeneConfig, MutualExclusivityTarget
from cfh.model.algorithm_result import AlgorithmResult
from cfh.model.fusion_event import FusionEvent
from cfh.model.fusion_feature import FusionFeature
from cfh.stats.cooccurrence_tests import (
    build_cooccurrence_contingency_table,
    fishers_cooccurrence_test,
)

ALGORITHM_NAME = "mutation_cooccurrence"
ALGORITHM_VERSION = "0.2.0"


def _target_label(target: MutualExclusivityTarget) -> str:
    label = f"{target.gene} {target.alteration_type}"
    return f"{label} ({target.protein_change})" if target.protein_change else label


def comparator_target_key(target: MutualExclusivityTarget) -> str:
    """Return the canonical sidecar key for a comparator target.

    The live-fetch layer supplies availability separately from alteration
    records: an empty record list is a valid, confirmed zero-call result,
    while an unavailable target must not be tested as though it had zero
    calls.  Use all fields that select a comparator so targets for the same
    gene remain independently addressable.
    """
    return json.dumps(
        [target.gene.upper(), target.alteration_type, target.protein_change],
        separators=(",", ":"),
    )


def _comparator_sample_ids(
    comparator_alterations: list[dict[str, Any]], target: MutualExclusivityTarget
) -> set[str]:
    gene = target.gene.upper()
    samples: set[str] = set()
    for record in comparator_alterations:
        if str(record.get("Gene", "")).upper() != gene:
            continue
        if record.get("Alteration_type") != target.alteration_type:
            continue
        if (
            target.protein_change is not None
            and record.get("Protein_change") != target.protein_change
        ):
            continue
        sample_id = record.get("Sample_id")
        if sample_id:
            samples.add(str(sample_id))
    return samples


def _no_op_result(warning: str) -> AlgorithmResult:
    return AlgorithmResult(
        Algorithm=ALGORITHM_NAME,
        Algorithm_version=ALGORITHM_VERSION,
        Parameters={},
        Summary={"targets": []},
        Tables={},
        Warnings=[warning],
        Created_at=datetime.now(timezone.utc),
    )


@register(ALGORITHM_NAME)
class MutationCooccurrenceAlgorithm(Algorithm):
    """Fisher's-exact co-occurrence/mutual-exclusivity test against each
    configured ``GeneConfig.mutual_exclusivity_targets`` entry. See the
    module docstring for the full opt-in/data-availability contract.

    Expected ``params`` keys (both required for a non-no-op result):
        cohort_sample_ids (list[str]): the caller-declared background universe.
        eligible_sample_ids_by_target (dict[str, list[str]], optional):
            samples assessed for both genes/assays, keyed by comparator_target_key.
            Missing targets fail closed. The live caller always supplies this.
            Offline callers omitting it assert that their universe is eligible.
        comparator_alterations (list[dict]): per-sample alteration records
            (``AlterationEvent.model_dump()`` shape: ``Sample_id``,
            ``Gene``, ``Alteration_type``, ``Protein_change``) for every
            configured comparator gene, already fetched by the caller.
        comparator_availability (dict[str, bool], optional): availability
            keyed by :func:`comparator_target_key`. When supplied, every
            configured target must have a ``True`` entry to be tested. This
            fails closed for partial live-fetch failures; omitting the
            sidecar preserves the original explicit-list behavior.
    """

    VERSION = ALGORITHM_VERSION

    def run(
        self,
        events: list[FusionEvent],
        features: list[FusionFeature],
        gene_config: GeneConfig,
        params: dict,
    ) -> AlgorithmResult:
        del features
        params = params or {}
        targets = gene_config.mutual_exclusivity_targets
        gene_label = gene_config.gene_symbol or gene_config.gene_pair or "This gene"
        if not targets:
            return _no_op_result(
                f"{gene_label} has no mutual_exclusivity_targets configured; "
                "co-occurrence/mutual-exclusivity analysis was skipped."
            )

        cohort_sample_ids = params.get("cohort_sample_ids")
        comparator_alterations = params.get("comparator_alterations")
        if not cohort_sample_ids:
            return _no_op_result(
                f"{gene_label} configures mutual_exclusivity_targets but no "
                "cohort_sample_ids were supplied to this run; co-occurrence/"
                "mutual-exclusivity analysis was skipped for this run."
            )

        comparator_availability = params.get("comparator_availability")
        has_availability_sidecar = comparator_availability is not None
        if comparator_alterations is None:
            return _no_op_result(
                f"{gene_label} configures mutual_exclusivity_targets but comparator "
                "alteration data were unavailable; co-occurrence/mutual-exclusivity "
                "analysis was skipped for this run."
            )

        cohort_sample_id_set = {str(sample_id) for sample_id in cohort_sample_ids}
        fusion_positive_sample_ids = {str(event.Sample_id) for event in events if event.Sample_id}
        warnings: list[str] = []
        unmatched = fusion_positive_sample_ids - cohort_sample_id_set
        if unmatched:
            warnings.append(
                f"{len(unmatched)} fusion-positive sample(s) were not found in the "
                "supplied cohort_sample_ids universe; they were excluded."
            )

        rows: list[dict[str, Any]] = []
        eligibility = params.get("eligible_sample_ids_by_target")
        for target in targets:
            target_key = comparator_target_key(target)
            if has_availability_sidecar and not (
                isinstance(comparator_availability, dict)
                and comparator_availability.get(target_key) is True
            ):
                warnings.append(
                    f"Comparator data for {_target_label(target)} were unavailable; "
                    "its co-occurrence/mutual-exclusivity test was skipped."
                )
                continue
            if eligibility is not None:
                if not isinstance(eligibility, dict) or target_key not in eligibility:
                    warnings.append(
                        f"Joint assay eligibility for {_target_label(target)} was unavailable; "
                        "its co-occurrence/mutual-exclusivity test was skipped."
                    )
                    continue
                eligible = cohort_sample_id_set & set(eligibility[target_key])
            else:
                eligible = cohort_sample_id_set
            if not eligible:
                warnings.append(
                    f"No jointly assessed samples for {_target_label(target)}; skipped."
                )
                continue
            all_comparator_samples = _comparator_sample_ids(comparator_alterations, target)
            comparator_samples = all_comparator_samples & eligible
            fusion_samples = fusion_positive_sample_ids & eligible
            table = build_cooccurrence_contingency_table(
                fusion_samples, comparator_samples, eligible
            )
            odds_ratio, p_value, direction = fishers_cooccurrence_test(table)
            rows.append(
                {
                    "comparator_gene": target.gene,
                    "alteration_type": target.alteration_type,
                    "protein_change": target.protein_change,
                    "label": _target_label(target),
                    "contingency_table": table,
                    "odds_ratio": odds_ratio,
                    "p_value": p_value,
                    "direction": direction,
                    "cohort_sample_count": len(eligible),
                    "supplied_cohort_sample_count": len(cohort_sample_id_set),
                    "excluded_cohort_sample_count": len(cohort_sample_id_set - eligible),
                    "excluded_fusion_positive_sample_count": len(
                        fusion_positive_sample_ids - eligible
                    ),
                    "excluded_comparator_sample_count": len(all_comparator_samples - eligible),
                    "counting_unit": "sample",
                    "eligibility_source": (
                        "joint_assay_metadata" if eligibility is not None else "caller_asserted"
                    ),
                    "eligible_sample_ids": sorted(eligible),
                    "fusion_positive_sample_count": len(fusion_samples),
                    "comparator_altered_sample_count": len(comparator_samples),
                    "both_altered_sample_count": table[0][0],
                }
            )
            if not comparator_samples:
                warnings.append(
                    "No cohort samples carried the configured comparator alteration "
                    f"({_target_label(target)}); its co-occurrence statistics reflect "
                    "zero comparator-positive samples."
                )

        return AlgorithmResult(
            Algorithm=ALGORITHM_NAME,
            Algorithm_version=ALGORITHM_VERSION,
            Parameters={
                "configured_targets": [target.model_dump() for target in targets],
                **(
                    {"comparator_availability": comparator_availability}
                    if has_availability_sidecar
                    else {}
                ),
            },
            Summary={"targets": rows},
            Tables={"cooccurrence_results": rows},
            Warnings=warnings,
            Created_at=datetime.now(timezone.utc),
        )
