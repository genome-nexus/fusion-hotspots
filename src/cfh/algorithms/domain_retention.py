"""Frame/domain-retention benchmark algorithm plugin."""

from __future__ import annotations

from cfh.algorithms.base import Algorithm
from cfh.algorithms.registry import register
from cfh.genes.registry import GeneConfig
from cfh.model.algorithm_result import AlgorithmResult
from cfh.model.fusion_event import FusionEvent
from cfh.model.fusion_feature import FusionFeature
from cfh.stats.adaptive_permutation import is_borderline, resolve_permutation_budget
from cfh.stats.breakpoint_tests import (
    build_frame_domain_contingency_table,
    domain_retention_descriptive_table,
    fishers_frame_domain_test,
    permutation_null_classifier,
    permutation_null_test,
)

_DEFAULT_FULL_N_PERMUTATIONS = 10_000


@register("domain_retention")
class DomainRetentionAlgorithm(Algorithm):
    """Test whether configured-domain retention is enriched in-frame.

    See :mod:`cfh.stats.adaptive_permutation` for the optional adaptive
    permutation-budget params (``adaptive``, ``n_permutations_small``,
    ``significance_threshold``, ``borderline_factor``); every one of them
    defaults to the existing non-adaptive behavior.
    """

    def run(
        self,
        events: list[FusionEvent],
        features: list[FusionFeature],
        gene_config: GeneConfig,
        params: dict,
    ) -> AlgorithmResult:
        params = params or {}
        table = build_frame_domain_contingency_table(events, features, gene_config)
        odds_ratio, fisher_p_value = fishers_frame_domain_test(table)

        seed = params.get("seed", 42)
        genome_nexus_client = params.get("genome_nexus_client")
        budget = resolve_permutation_budget(params, default_full_n=_DEFAULT_FULL_N_PERMUTATIONS)

        n_permutations = budget["full_n"] if not budget["adaptive"] else budget["small_n"]
        permutation_p_value, observed_rate, null_rates = permutation_null_test(
            events,
            features,
            gene_config,
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
                seed=seed,
                n_permutations=n_permutations,
                genome_nexus_client=genome_nexus_client,
            )

        summary = {
            "fisher_odds_ratio": odds_ratio,
            "fisher_p_value": fisher_p_value,
            "permutation_empirical_p_value": permutation_p_value,
            "observed_in_frame_retention_rate": observed_rate,
            "permutation_null_classifier": permutation_null_classifier(features, gene_config),
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
            Algorithm="domain_retention",
            Algorithm_version="0.3.0",
            Parameters={
                "seed": seed,
                "n_permutations": n_permutations,
                "adaptive": budget["adaptive"],
            },
            Summary=summary,
            Tables={
                "frame_domain_contingency_table": table,
                "domain_retention_descriptives": domain_retention_descriptive_table(
                    features, gene_config
                ),
                "permutation_null_retention_rates": list(null_rates),
            },
        )
