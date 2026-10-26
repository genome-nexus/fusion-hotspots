"""Junction-concentration test: do a gene's breakpoints pile up at one exon?

Replaces the cutpoint and window scans as the gene's breakpoint-position
test in the cohort-scan FDR family. Those scans permute retained/lost labels
whose values are fixed by breakpoint position, so a clean split at the domain
edge exists by construction and even uniformly random breakpoints score
p ~ 1e-4. This test instead compares each junction exon's count with
breakpoints redrawn along the gene; see :mod:`cfh.stats.junction_concentration`
for the null and its panel-capture caveat.
"""

from __future__ import annotations

from datetime import datetime, timezone

from cfh.algorithms.base import Algorithm
from cfh.algorithms.registry import register
from cfh.genes.registry import GeneConfig
from cfh.mapping.genome_nexus_source import cds_bounds_from_utrs, parse_canonical_transcript
from cfh.model.algorithm_result import AlgorithmResult
from cfh.model.fusion_event import FusionEvent
from cfh.model.fusion_feature import FusionFeature
from cfh.stats.adaptive_permutation import resolve_permutation_budget, run_permutation_budget
from cfh.stats.junction_concentration import (
    DEFAULT_MIN_EVENTS,
    PANEL_CAPTURE_CAVEAT,
    ROLES,
    junction_concentration_test,
    junction_probabilities,
)

ALGORITHM_NAME = "junction_concentration"
ALGORITHM_VERSION = "0.1.0"


def _not_determinable(
    reason: str, n_events: int, params: dict, n_outside_transcript: int | None = None
) -> AlgorithmResult:
    return AlgorithmResult(
        Algorithm=ALGORITHM_NAME,
        Algorithm_version=ALGORITHM_VERSION,
        Parameters=params,
        Summary={
            "determinable": False,
            "reason": reason,
            "n_events_analyzed": n_events,
            "n_outside_transcript": n_outside_transcript,
            "permutation_empirical_p_value": None,
        },
        Tables={"junction_counts": []},
        Warnings=[reason],
        Created_at=datetime.now(timezone.utc),
    )


@register(ALGORITHM_NAME)
class JunctionConcentrationAlgorithm(Algorithm):
    """Test whether any junction exon holds far more of the gene's fusions
    than breakpoints drawn uniformly along the gene would put there.

    Expected ``params`` keys (all optional except the client):
        genome_nexus_client (GenomeNexusClient): source of the canonical
            transcript's exon coordinates; without it the test is not
            determinable.
        seed (int): RNG seed, default 42.
        n_permutations (int): null replicates, default 10000; see
            :mod:`cfh.stats.adaptive_permutation` for the adaptive keys.
        min_events (int): fewest events with a known junction, default
            :data:`~cfh.stats.junction_concentration.DEFAULT_MIN_EVENTS`.
        genomic_breakpoints (dict): optional ``{event_id: {"position": int}}``;
            events whose breakpoint lies outside the transcript are left out.
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
        min_events = int(params.get("min_events", DEFAULT_MIN_EVENTS))
        budget = resolve_permutation_budget(params, default_full_n=10_000)
        if gene_config.gene_symbol is None:
            return _not_determinable(
                "junction concentration needs a single target gene", 0, {"seed": seed}
            )
        gene = gene_config.gene_symbol.upper()
        event_ids = {event.Event_id for event in events}
        candidates = [
            feature
            for feature in features
            if feature.Event_id in event_ids
            and feature.Gene.upper() == gene
            and feature.Role in ROLES
            and feature.Breakpoint_exon is not None
        ]
        parameters = {
            "seed": seed,
            "n_permutations": budget["small_n"] if budget["adaptive"] else budget["full_n"],
            "min_events": min_events,
            "adaptive": budget["adaptive"],
        }
        client = params.get("genome_nexus_client")
        if client is None:
            return _not_determinable(
                "no Genome Nexus client to fetch exon coordinates for the null",
                len(candidates),
                parameters,
            )
        transcript = parse_canonical_transcript(
            client.fetch_canonical_transcript(gene_config.gene_symbol)
        )
        if not transcript.exons:
            return _not_determinable(
                "the canonical transcript has no exon coordinates", len(candidates), parameters
            )

        # The null only places breakpoints inside the transcript, so a
        # breakpoint outside it (e.g. upstream of a short canonical isoform),
        # which the mapper clamps to the first or last exon, is left out.
        span_start = min(exon.start for exon in transcript.exons)
        span_end = max(exon.end for exon in transcript.exons)
        genomic_breakpoints = params.get("genomic_breakpoints") or {}
        outside = {
            feature.Event_id
            for feature in candidates
            if (position := (genomic_breakpoints.get(feature.Event_id) or {}).get("position"))
            is not None
            and not span_start <= int(position) <= span_end
        }
        junctions = [
            (str(feature.Role), int(feature.Breakpoint_exon))
            for feature in candidates
            if feature.Event_id not in outside and feature.Breakpoint_exon is not None
        ]
        if len(junctions) < min_events:
            return _not_determinable(
                f"{len(junctions)} event(s) with a known junction exon inside the "
                f"transcript; at least {min_events} are needed",
                len(junctions),
                parameters,
                len(outside),
            )
        cds_min, cds_max = cds_bounds_from_utrs(transcript.utrs)
        probabilities = {
            role: junction_probabilities(
                transcript.exons, role, cds_min_genomic=cds_min, cds_max_genomic=cds_max
            )
            for role in sorted({role for role, _ in junctions})
        }

        result, n_permutations, escalated = run_permutation_budget(
            lambda n: junction_concentration_test(
                junctions, probabilities, seed=seed, n_permutations=n
            ),
            lambda outcome: outcome["permutation_empirical_p_value"],
            budget,
        )
        summary = {
            "determinable": True,
            "reason": None,
            "n_events_analyzed": result["n_events"],
            "n_outside_transcript": len(outside),
            "top_junction": result["top_junction"],
            "permutation_empirical_p_value": result["permutation_empirical_p_value"],
            "null_model": "breakpoints uniform across the gene's genomic span, role kept",
            "panel_capture_caveat": PANEL_CAPTURE_CAVEAT,
        }
        if budget["adaptive"]:
            summary["adaptive_permutations"] = {
                "enabled": True,
                "n_permutations_small": budget["small_n"],
                "n_permutations_full": budget["full_n"],
                "n_permutations_used": n_permutations,
                "n_permutations_max": budget["max_n"],
                "escalated_to_full": escalated,
            }
        return AlgorithmResult(
            Algorithm=ALGORITHM_NAME,
            Algorithm_version=ALGORITHM_VERSION,
            Parameters={**parameters, "n_permutations": n_permutations},
            Summary=summary,
            Tables={"junction_counts": result["junction_counts"]},
            Warnings=[],
            Created_at=datetime.now(timezone.utc),
        )
