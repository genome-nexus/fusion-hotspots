"""Exon-retention enrichment analysis for normalized fusion features."""

from __future__ import annotations

import re
from datetime import datetime, timezone
from typing import Any

from cfh.algorithms.base import Algorithm
from cfh.algorithms.registry import register
from cfh.genes.registry import GeneConfig
from cfh.model.algorithm_result import AlgorithmResult
from cfh.model.fusion_event import FusionEvent
from cfh.model.fusion_feature import FusionFeature


def _target_exon(params: dict[str, Any], gene_config: GeneConfig) -> int | None:
    """Resolve an exon from explicit params or the optional config hint."""
    candidate = params.get("target_exon", params.get("exon"))
    if candidate is None:
        candidate = gene_config.expected_retained_exon_hint
    if candidate is None:
        return None
    if isinstance(candidate, int):
        return candidate
    match = re.search(r"\d+", str(candidate))
    return int(match.group()) if match else None


def _retains_exon(feature: FusionFeature, target_exon: int) -> bool:
    """Whether a feature records retention of the requested exon."""
    return any(str(exon) == str(target_exon) for exon in feature.Retained_exons or [])


@register("exon_retention")
class ExonRetentionAnalysis(Algorithm):
    """Measure the fraction of events retaining a selected exon."""

    VERSION = "1.0.0"

    def run(
        self,
        events: list[FusionEvent],
        features: list[FusionFeature],
        gene_config: GeneConfig,
        params: dict,
    ) -> AlgorithmResult:
        """Return event-level retention counts and fraction for one target exon."""
        target_exon = _target_exon(params, gene_config)
        if target_exon is None:
            return AlgorithmResult(
                Algorithm="exon_retention",
                Algorithm_version=self.VERSION,
                Parameters={},
                Summary={"retained_fraction": None},
                Tables={"Exon_retention": []},
                Warnings=["No target exon was supplied in params or GeneConfig."],
                Created_at=datetime.now(timezone.utc),
            )

        event_ids = {event.Event_id for event in events}
        retained_event_ids = {
            feature.Event_id
            for feature in features
            if feature.Event_id in event_ids
            and feature.Gene.upper() == gene_config.gene_symbol.upper()
            and _retains_exon(feature, target_exon)
        }
        retained_count = len(retained_event_ids)
        total_count = len(events)
        retained_fraction = retained_count / total_count if total_count else 0.0
        table = {
            "Exon": target_exon,
            "Retained_event_count": retained_count,
            "Total_event_count": total_count,
            "Retained_fraction": retained_fraction,
        }

        return AlgorithmResult(
            Algorithm="exon_retention",
            Algorithm_version=self.VERSION,
            Parameters={"target_exon": target_exon},
            Summary={
                "target_exon": target_exon,
                "retained_event_count": retained_count,
                "total_event_count": total_count,
                "retained_fraction": retained_fraction,
            },
            Tables={"Exon_retention": [table]},
            Warnings=[],
            Created_at=datetime.now(timezone.utc),
        )
