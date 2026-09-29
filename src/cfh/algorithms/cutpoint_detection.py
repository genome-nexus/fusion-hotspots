"""Gene-agnostic, frequency/recurrence-only cutpoint (boundary) detection.

No external oncogenicity label and no OncoKB dependency: the outcome used
here is domain-retention status (retained vs. lost/disrupted), already
computed by the domain-retention algorithm and stored per-event via
``FusionFeature.Domain_retention_flags``. The scanning axis is each event's
already-mapped breakpoint protein position
(``FusionFeature.Junction_position_aa``); recurrence is simply how many
independent events land at or near a given position. No new network call is
introduced: an optional ``GenomeNexusClient`` may be passed in ``params``
purely to compare the inferred cutpoint against Pfam domain boundaries the
pipeline already fetches elsewhere, and the algorithm runs fully offline
without one.
"""

from __future__ import annotations

from datetime import datetime, timezone
from typing import Optional

from cfh.algorithms.base import Algorithm
from cfh.algorithms.registry import register
from cfh.genes.registry import GeneConfig
from cfh.mapping.genome_nexus_source import GenomeNexusClient, resolve_domains
from cfh.model.algorithm_result import AlgorithmResult
from cfh.model.fusion_event import FusionEvent
from cfh.model.fusion_feature import FusionFeature
from cfh.stats.adaptive_permutation import is_borderline, resolve_permutation_budget
from cfh.stats.breakpoint_tests import gene_breakpoint_domain_status_records
from cfh.stats.cutpoint_scan import detect_cutpoint

ALGORITHM_NAME = "cutpoint_detection"
ALGORITHM_VERSION = "0.1.0"
_DEFAULT_FULL_N_PERMUTATIONS = 10_000


def _nearest_boundary_comparison(cutpoint: int, boundaries: list[int]) -> Optional[dict]:
    """Nearest known domain boundary to ``cutpoint``, with a signed distance.

    Reported unconditionally whenever any ``boundaries`` are known -- including
    when ``cutpoint`` falls outside every configured domain span (e.g.
    downstream of the last domain's end), which is itself a meaningful result,
    not an absence of one. ``distance_aa`` is signed: positive means
    ``cutpoint`` is C-terminal (higher aa) of the nearest boundary, negative
    means N-terminal of it. Only ``None`` when no boundaries are known at all.
    """
    if not boundaries:
        return None
    nearest = min(boundaries, key=lambda boundary: abs(boundary - cutpoint))
    return {
        "nearest_known_domain_boundary_aa": nearest,
        "distance_aa": cutpoint - nearest,
    }


def _known_domain_boundaries(gene_config: GeneConfig, params: dict) -> list[int]:
    """Known Pfam/domain boundary positions to validate the inferred cutpoint against.

    Prefers an explicit ``params["domain_boundaries"]`` list (e.g. sourced
    from a committed fixture) so the benchmark stays reproducible offline.
    Falls back to an optional ``params["genome_nexus_client"]`` -- the same
    client type ``domain_retention`` already accepts -- resolved through the
    same :func:`~cfh.mapping.genome_nexus_source.resolve_domains` Genome
    Nexus-then-UniProt fallback chain ``domain_retention`` already relies on
    (via ``params["uniprot_source"]``, primarily for tests), so a gene whose
    canonical transcript has no Genome Nexus Pfam annotation (observed live
    for RET) still gets a real boundary comparison instead of silently
    falling back to no boundaries at all. Neither is required; with neither
    supplied (and no UniProt accession available either) the comparison is
    simply omitted.
    """
    explicit = params.get("domain_boundaries")
    if explicit:
        return sorted({int(boundary) for boundary in explicit})

    client: Optional[GenomeNexusClient] = params.get("genome_nexus_client")
    if client is None or not gene_config.gene_symbol:
        return []

    domains = resolve_domains(
        gene_config.gene_symbol,
        gene_config.protein_id,
        genome_nexus_client=client,
        uniprot_source=params.get("uniprot_source"),
    )
    boundaries: set[int] = set()
    for domain in domains:
        boundaries.add(domain.start_aa)
        boundaries.add(domain.end_aa)
    return sorted(boundaries)


@register(ALGORITHM_NAME)
class CutpointDetectionAlgorithm(Algorithm):
    """Scan a gene's observed breakpoints for the protein-position cutpoint
    that best separates domain-retained from lost/disrupted events, with a
    permutation-corrected empirical p-value for having scanned many
    candidate cutpoints.

    Expected ``params`` keys (all optional):
        seed (int): permutation RNG seed, default 42.
        n_permutations (int): number of label permutations, default 10000.
        domain_boundaries (list[int]): known domain boundary aa positions to
            compare the inferred cutpoint against.
        genome_nexus_client (GenomeNexusClient): used only as a fallback
            source for ``domain_boundaries`` when that param is omitted
            (Genome Nexus first, UniProt second -- see
            :func:`~cfh.mapping.genome_nexus_source.resolve_domains`).
        uniprot_source: optional override for the UniProt fallback used by
            ``genome_nexus_client``-based boundary resolution (mainly for
            tests); ignored otherwise.

    ``known_domain_boundary_comparison`` in the result ``Summary`` is
    ``{"nearest_known_domain_boundary_aa": ..., "distance_aa": ...}``
    whenever any boundary is known -- including when the inferred cutpoint
    falls outside every configured domain span -- and ``None`` only when no
    boundaries could be resolved at all. ``distance_aa`` is signed: positive
    means the cutpoint is C-terminal of the nearest boundary, negative means
    N-terminal of it.

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
        seed = params.get("seed", 42)
        budget = resolve_permutation_budget(params, default_full_n=10_000)

        records = gene_breakpoint_domain_status_records(events, features, gene_config)
        positions = [position for position, _ in records]
        statuses = [status for _, status in records]

        n_permutations = budget["full_n"] if not budget["adaptive"] else budget["small_n"]
        scan_result = detect_cutpoint(
            positions,
            statuses,
            seed=seed,
            n_permutations=n_permutations,
        )
        escalated = False
        if (
            budget["adaptive"]
            and scan_result["determinable"]
            and is_borderline(
                scan_result["corrected_p_value"],
                threshold=budget["threshold"],
                factor=budget["factor"],
            )
        ):
            escalated = True
            n_permutations = budget["full_n"]
            scan_result = detect_cutpoint(
                positions,
                statuses,
                seed=seed,
                n_permutations=n_permutations,
            )

        boundary_comparison = None
        warnings: list[str] = []
        if scan_result["determinable"]:
            boundaries = _known_domain_boundaries(gene_config, params)
            boundary_comparison = _nearest_boundary_comparison(
                scan_result["best_cutpoint"], boundaries
            )
        else:
            warnings.append(scan_result["reason"])

        summary = {
            "determinable": scan_result["determinable"],
            "reason": scan_result["reason"],
            "n_events_analyzed": len(records),
            "inferred_cutpoint_aa": scan_result["best_cutpoint"],
            "observed_statistic_neg_log10_p": scan_result["observed_statistic"],
            "observed_p_value": scan_result["observed_p_value"],
            "observed_odds_ratio": scan_result["observed_odds_ratio"],
            "corrected_p_value": scan_result["corrected_p_value"],
            "known_domain_boundary_comparison": boundary_comparison,
        }
        if "permutation_resolution" in budget:
            summary["permutation_resolution"] = budget["permutation_resolution"]
        if budget["adaptive"]:
            summary["adaptive_permutations"] = {
                "enabled": True,
                "n_permutations_small": budget["small_n"],
                "n_permutations_full": budget["full_n"],
                "n_permutations_used": n_permutations,
                "escalated_to_full": escalated,
            }

        return AlgorithmResult(
            Algorithm=ALGORITHM_NAME,
            Algorithm_version=ALGORITHM_VERSION,
            Parameters={
                "seed": seed,
                "n_permutations": n_permutations,
                "adaptive": budget["adaptive"],
            },
            Summary=summary,
            Tables={"cutpoint_scan": scan_result["scan"]},
            Warnings=warnings,
            Created_at=datetime.now(timezone.utc),
        )
