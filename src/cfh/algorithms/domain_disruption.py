"""Frame/domain-disruption benchmark algorithm plugin -- the inverse of
:mod:`cfh.algorithms.domain_retention`.

Some oncogenic fusions require a configured domain to be LOST/DISRUPTED
rather than retained (e.g. an autoinhibitory regulatory module excised by
the breakpoint). This reuses the exact same Fisher's-exact/permutation
machinery as domain retention, just testing domain *exclusion* instead of
domain *retention* as the enriched outcome among in-frame fusions.

Opt-in per gene via ``GeneConfig.disruption_required_domains``: a gene that
leaves this field unset (its default, empty list) produces a no-op result
with no statistics computed, the same graceful-skip pattern already used by
``exon_retention``/``joint_partner``/``confidence_stats`` for genes that
don't configure their respective optional fields.
"""

from __future__ import annotations

from datetime import datetime, timezone

from cfh.algorithms.base import Algorithm
from cfh.algorithms.registry import register
from cfh.genes.registry import GeneConfig
from cfh.model.algorithm_result import AlgorithmResult
from cfh.model.fusion_event import FusionEvent
from cfh.model.fusion_feature import FusionFeature
from cfh.stats.adaptive_permutation import is_borderline, resolve_permutation_budget
from cfh.stats.breakpoint_tests import (
    DISRUPTED_STATUSES,
    build_frame_domain_contingency_table,
    domain_retention_descriptive_table,
    fishers_frame_domain_test,
    permutation_null_classifier,
    permutation_null_test,
)
from cfh.stats.stratification import stratified_cmh_summary, stratified_frame_domain_tables

VERSION = "0.4.0"
_DEFAULT_FULL_N_PERMUTATIONS = 10_000


@register("domain_disruption")
class DomainDisruptionAlgorithm(Algorithm):
    """Test whether configured-domain loss/disruption is enriched in-frame."""

    def run(
        self,
        events: list[FusionEvent],
        features: list[FusionFeature],
        gene_config: GeneConfig,
        params: dict,
    ) -> AlgorithmResult:
        params = params or {}
        domains = gene_config.disruption_required_domains
        if not domains:
            return AlgorithmResult(
                Algorithm="domain_disruption",
                Algorithm_version=VERSION,
                Parameters={},
                Summary={
                    "fisher_odds_ratio": None,
                    "fisher_p_value": None,
                    "permutation_empirical_p_value": None,
                    "observed_in_frame_disruption_rate": None,
                },
                Tables={},
                Warnings=[
                    f"{gene_config.gene_symbol or gene_config.gene_pair} has no "
                    "disruption_required_domains configured; domain-disruption "
                    "analysis was skipped."
                ],
                Created_at=datetime.now(timezone.utc),
            )

        table = build_frame_domain_contingency_table(
            events, features, gene_config, domains=domains, hit_statuses=DISRUPTED_STATUSES
        )
        odds_ratio, fisher_p_value = fishers_frame_domain_test(table)
        tables_by_tumor_type = stratified_frame_domain_tables(
            events, features, gene_config, domains=domains, hit_statuses=DISRUPTED_STATUSES
        )

        seed = params.get("seed", 42)
        genome_nexus_client = params.get("genome_nexus_client")
        budget = resolve_permutation_budget(params, default_full_n=_DEFAULT_FULL_N_PERMUTATIONS)

        n_permutations = budget["full_n"] if not budget["adaptive"] else budget["small_n"]
        permutation_p_value, observed_rate, null_rates = permutation_null_test(
            events,
            features,
            gene_config,
            domains=domains,
            hit_statuses=DISRUPTED_STATUSES,
            seed=seed,
            n_permutations=n_permutations,
            genome_nexus_client=genome_nexus_client,
        )
        escalated = False
        if budget["adaptive"] and is_borderline(
            permutation_p_value, threshold=budget["threshold"], factor=budget["factor"]
        ):
            escalated = True
            n_permutations = budget["full_n"]
            permutation_p_value, observed_rate, null_rates = permutation_null_test(
                events,
                features,
                gene_config,
                domains=domains,
                hit_statuses=DISRUPTED_STATUSES,
                seed=seed,
                n_permutations=n_permutations,
                genome_nexus_client=genome_nexus_client,
            )

        summary = {
            "fisher_odds_ratio": odds_ratio,
            "fisher_p_value": fisher_p_value,
            "permutation_empirical_p_value": permutation_p_value,
            "observed_in_frame_disruption_rate": observed_rate,
            "permutation_null_classifier": permutation_null_classifier(
                features, gene_config, domains=domains
            ),
            # Reported alongside the pooled Fisher test; not in the FDR family.
            "tumor_type_stratified": stratified_cmh_summary(tables_by_tumor_type),
        }
        if budget["adaptive"]:
            summary["adaptive_permutations"] = {
                "enabled": True,
                "n_permutations_small": budget["small_n"],
                "n_permutations_full": budget["full_n"],
                "n_permutations_used": n_permutations,
                "escalated_to_full": escalated,
            }

        return AlgorithmResult(
            Algorithm="domain_disruption",
            Algorithm_version=VERSION,
            Parameters={
                "seed": seed,
                "n_permutations": n_permutations,
                "adaptive": budget["adaptive"],
            },
            Summary=summary,
            Tables={
                "frame_domain_contingency_table": table,
                "frame_domain_contingency_tables_by_tumor_type": [
                    {"stratum": stratum, "table": stratum_table}
                    for stratum, stratum_table in tables_by_tumor_type.items()
                ],
                "domain_disruption_descriptives": domain_retention_descriptive_table(
                    features, gene_config, domains=domains
                ),
                "permutation_null_disruption_rates": list(null_rates),
            },
        )
