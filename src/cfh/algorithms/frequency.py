"""Partner-gene frequency analysis for normalized fusion events."""

from __future__ import annotations

from collections import Counter
from datetime import datetime, timezone

from cfh.algorithms.base import Algorithm
from cfh.algorithms.registry import register
from cfh.genes.registry import GeneConfig
from cfh.model.algorithm_result import AlgorithmResult
from cfh.model.fusion_event import FusionEvent
from cfh.model.fusion_feature import FusionFeature


def _partner_gene(event: FusionEvent, target_gene: str) -> str:
    """Return the gene paired with ``target_gene`` for one event.

    Events that do not identify a partner are retained as ``"unknown"`` so
    the frequency table remains an accounting of every analyzed event.
    """
    site1 = event.Site1_gene
    site2 = event.Site2_gene
    target = target_gene.upper()

    if site1 and site1.upper() == target and site2:
        return site2
    if site2 and site2.upper() == target and site1:
        return site1
    return "unknown"


def _deduplication_key(event: FusionEvent) -> str:
    """Use patient identity when present, preserving events without one."""
    return event.Patient_id or event.Sample_id or event.Event_id


@register("frequency")
class FrequencyAnalysis(Algorithm):
    """Count fusion partner genes for the configured target gene."""

    VERSION = "1.0.0"

    def run(
        self,
        events: list[FusionEvent],
        features: list[FusionFeature],
        gene_config: GeneConfig,
        params: dict,
    ) -> AlgorithmResult:
        """Return partner-gene counts, optionally retaining one event per patient."""
        del features
        dedup_by_patient = bool(params.get("dedup_by_patient", False))
        analyzed_events = events
        if dedup_by_patient:
            seen: set[str] = set()
            analyzed_events = []
            for event in events:
                key = _deduplication_key(event)
                if key not in seen:
                    seen.add(key)
                    analyzed_events.append(event)

        counts = Counter(_partner_gene(event, gene_config.gene_symbol) for event in analyzed_events)
        table = [
            {"Partner_gene": partner_gene, "Event_count": count}
            for partner_gene, count in sorted(counts.items())
        ]

        return AlgorithmResult(
            Algorithm="frequency",
            Algorithm_version=self.VERSION,
            Parameters={"dedup_by_patient": dedup_by_patient},
            Summary={
                "input_event_count": len(events),
                "analyzed_event_count": len(analyzed_events),
                "unique_partner_gene_count": len(counts),
            },
            Tables={"Partner_gene_counts": table},
            Warnings=[],
            Created_at=datetime.now(timezone.utc),
        )
