"""Prototype live cBioPortal-to-domain-retention benchmark pipeline."""

from __future__ import annotations

import csv
import json
import math
import re
import subprocess
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
import requests

from cfh.algorithms.confidence_stats import resolve_confidence_stats_params
from cfh.algorithms.frequency import FrequencyAnalysis
from cfh.algorithms.registry import list_algorithms
from cfh.genes.registry import GeneConfig, derive_gene_config_defaults, load_gene_config
from cfh.ingestion import cbioportal_api
from cfh.mapping.domain_source import ProteinDomain
from cfh.mapping.feature_mapper import map_event
from cfh.mapping.fusion_annotation_crosscheck import (
    FusionAnnotationCrosscheckResult,
    crosscheck_breakpoint_protein_position,
)
from cfh.mapping.genome_nexus_source import (
    CanonicalTranscript,
    GenomeNexusClient,
    GenomeNexusGeneNotFound,
    cds_bounds_from_utrs,
    exon_protein_boundaries,
    gene_track_from_canonical_transcript,
    parse_canonical_transcript,
    resolve_domains,
)
from cfh.mapping.transcript_source import resolve_breakpoint_protein_position
from cfh.model.algorithm_result import AlgorithmResult
from cfh.model.fusion_event import FusionEvent
from cfh.model.fusion_feature import FusionFeature
from cfh.normalization.alteration_normalizer import (
    normalize_discrete_copy_number,
    normalize_mutations,
)
from cfh.normalization.event_normalizer import normalize
from cfh.orchestrator.run import run_algorithms
from cfh.reporting.domain_names import (
    configured_domain_names,
    domain_interpretation_sentence,
    format_domain_names,
)
from cfh.reporting.exon_labels import exon_label_for_protein_position
from cfh.reporting.expression_association_text import expression_association_sentence
from cfh.reporting.fusion_schematic import (
    render_fusion_schematic_svg,
    render_intragenic_deletion_schematic_svg,
    render_position_axis_svg,
)
from cfh.reporting.mutual_exclusivity_text import mutual_exclusivity_report_lines
from cfh.reporting.palette import (
    AXIS_COLOR,
    BREAKPOINT_COLOR,
    DOMAIN_HIGHLIGHT_COLOR,
    LOST_COLOR,
    REFERENCE_BAR_COLOR,
    RETAINED_COLOR,
    TRUNCATED_COLOR,
    deterministic_color,
)
from cfh.reporting.pdf import render_pdf_report
from cfh.reporting.svg_utils import escape_xml_text
from cfh.stats.breakpoint_tests import build_frame_domain_contingency_table
from cfh.studies.registry import StudyConfig, load_study_config

_DELETION_EVENT_INFO_PATTERN = re.compile(
    r"deletion of (\d+) exons?\s*:\s*(in frame|out of frame)", re.IGNORECASE
)


class RealBenchmarkError(RuntimeError):
    """Base class for expected, user-actionable benchmark failures."""


class RealBenchmarkInputError(RealBenchmarkError):
    """Raised when the requested gene/study cannot define a valid run."""


class RealBenchmarkNetworkError(RealBenchmarkError):
    """Raised when a required public data service cannot be reached."""


_TARGET_LOCUS_TOLERANCE_BP = 1_000_000
"""Permit a nearby assembly/annotation offset while rejecting remote loci.

Some imported cBioPortal rows omit or mix genome-build metadata. A one-Mb
envelope retains legitimate target-gene coordinates that differ slightly
between builds, while easily rejecting partner positions many megabases away.
"""


@dataclass
class RealBenchmarkRun:
    gene_symbol: str
    study_id: str
    molecular_profile_id: str
    retrieved_at: datetime
    raw_structural_variant_count: int
    events: list[FusionEvent]
    features: list[FusionFeature]
    rows: list[dict]
    results: list[AlgorithmResult]
    summary: dict
    warnings: list[str]
    endpoints: list[str] = field(default_factory=list)
    reference: dict | None = None
    gene_track: dict | None = None
    """Gene-agnostic protein-layout summary (protein length, full domain
    map, exon-to-residue boundaries) for the schematic visualization; see
    :func:`cfh.mapping.genome_nexus_source.gene_track_from_canonical_transcript`.
    ``None`` when no canonical-transcript mapping was available for this
    gene, or there was no target-gene data to visualize."""
    intragenic_deletions: list[dict] = field(default_factory=list)
    """Same-gene (Site1==Site2==target) intragenic-deletion-style SV
    records (e.g. ``"Deletion of N exons"`` annotations), analogous to a
    panel-C schematic. Distinct from ``events``/``features``: these never
    feed any algorithm's statistics, they are visualization-only."""
    is_gene_pair: bool = False
    """``True`` for a run produced by :func:`run_gene_pair_benchmark` (a
    ``gene_pair``-configured joint-partner analysis, e.g. EML4-ALK) rather
    than the default single-gene domain-retention pipeline. ``write_outputs``
    uses this to route to the gene-pair-shaped artifact writer instead of
    the domain-retention report writer, which assumes fields (domain
    accessions, breakpoint positions, ...) a gene-pair run never has."""


class _ResolvedDomainSource:
    """Expose already-resolved Genome Nexus domains to ``map_event``."""

    def __init__(self, domains: list[ProteinDomain]):
        self.domains = domains

    def fetch(self, _accession: str) -> list[ProteinDomain]:
        return self.domains


def _retained_exon_ranks(
    canonical: CanonicalTranscript, junction_position_aa: int, role: str
) -> list[int]:
    """Return exons wholly present on the retained side of a fusion junction."""
    cds_min, cds_max = cds_bounds_from_utrs(canonical.utrs)
    boundaries = exon_protein_boundaries(
        canonical.exons,
        cds_min_genomic=cds_min,
        cds_max_genomic=cds_max,
    )
    if role == "five_prime":
        return [
            boundary.exon_rank for boundary in boundaries if boundary.end_aa <= junction_position_aa
        ]
    if role == "three_prime":
        return [
            boundary.exon_rank
            for boundary in boundaries
            if boundary.start_aa >= junction_position_aa
        ]
    return []


def _is_target_protein_fusion(event: FusionEvent, target_gene: str) -> bool:
    genes = {str(event.Site1_gene or "").upper(), str(event.Site2_gene or "").upper()}
    return (
        target_gene.upper() in genes
        and event.Is_protein_fusion is True
        and "fusion" in str(event.Event_info or "").lower()
    )


def _target_locus(canonical: CanonicalTranscript) -> tuple[int, int]:
    """Return the inclusive genomic footprint spanned by a target's exons."""
    if not canonical.exons:
        raise ValueError("Genome Nexus returned no exon coordinates for target-locus validation")
    return (
        min(exon.start for exon in canonical.exons),
        max(exon.end for exon in canonical.exons),
    )


def _target_breakpoint(row: dict, target_gene: str, target_locus: tuple[int, int]) -> int:
    """Choose the only site position compatible with the target's GN locus.

    cBioPortal structural-variant rows can pair a site label with the other
    fusion partner's coordinate. Labels alone are therefore not a safe
    breakpoint selector. The canonical transcript already fetched from
    Genome Nexus supplies the target's exon-spanned genomic footprint; use
    that to select a position independently of its reported site label.
    """
    locus_start, locus_end = target_locus
    exact_candidates: list[tuple[str, int]] = []
    nearby_candidates: list[tuple[str, int]] = []
    invalid_values: list[str] = []
    for site in ("Site1", "Site2"):
        value = row.get(f"{site}_Position")
        if value is None or pd.isna(value):
            continue
        try:
            position = int(value)
        except (TypeError, ValueError):
            invalid_values.append(f"{site}={value!r}")
            continue
        if locus_start <= position <= locus_end:
            exact_candidates.append((site, position))
        elif (
            locus_start - _TARGET_LOCUS_TOLERANCE_BP
            <= position
            <= locus_end + _TARGET_LOCUS_TOLERANCE_BP
        ):
            nearby_candidates.append((site, position))

    if len(exact_candidates) == 1:
        return exact_candidates[0][1]
    site_positions = ", ".join(
        f"{site}={row.get(f'{site}_Position')!r}" for site in ("Site1", "Site2")
    )
    if len(exact_candidates) == 2:
        raise ValueError(
            f"{target_gene} fusion has ambiguous genomic breakpoints within target locus "
            f"{locus_start}-{locus_end}: {site_positions}"
        )
    if len(nearby_candidates) == 1:
        return nearby_candidates[0][1]
    if not invalid_values and all(
        row.get(f"{site}_Position") is None or pd.isna(row.get(f"{site}_Position"))
        for site in ("Site1", "Site2")
    ):
        raise ValueError(f"{target_gene} fusion has no genomic breakpoint")
    invalid_note = f"; invalid values: {', '.join(invalid_values)}" if invalid_values else ""
    raise ValueError(
        f"{target_gene} fusion has no site breakpoint within target locus "
        f"{locus_start}-{locus_end}: {site_positions}{invalid_note}"
    )


def _target_role(event: FusionEvent, target_gene: str) -> str:
    target = target_gene.upper()
    if str(event.Five_prime_gene or "").upper() == target:
        return "five_prime"
    if str(event.Three_prime_gene or "").upper() == target:
        return "three_prime"
    raise ValueError(
        f"could not determine 5'/3' role for {target_gene} in {event.Event_id}; "
        f"Event_Info={event.Event_info!r}"
    )


def _partner(event: FusionEvent, target_gene: str) -> str:
    target = target_gene.upper()
    if str(event.Site1_gene or "").upper() == target:
        return event.Site2_gene or "unknown"
    if str(event.Site2_gene or "").upper() == target:
        return event.Site1_gene or "unknown"
    return "unknown"


def _source_frame_status(text: object) -> str | None:
    """Read an explicit frame call without applying production normalization rules."""
    value = str(text or "").strip().lower()
    if not value:
        return None
    if value in {"na", "n/a", "unknown"}:
        return None
    if "out-of-frame" in value or "out of frame" in value:
        return "out-of-frame"
    if "in-frame" in value or "in frame" in value:
        return "in-frame"
    return None


def _pretty_domain_names(config: GeneConfig) -> dict[str, str]:
    """Map a configured domain's source accession to its human-readable
    ``name`` (e.g. ``"PF07714" -> "Protein kinase domain"``), for relabeling
    the raw Pfam-id-keyed domain map in :func:`_build_gene_track`. Domains
    outside ``key_domains``/``disruption_required_domains`` keep their raw
    accession as their display name -- this only ever prettifies domains
    the gene's own config already names, never invents a label.
    """
    return {
        domain.accession: domain.name
        for domain in [*config.key_domains, *config.disruption_required_domains]
        if domain.accession
    }


def _build_gene_track(config: GeneConfig, client: GenomeNexusClient) -> dict | None:
    """Best-effort protein-layout summary for the fusion schematic.

    Reuses the same canonical-transcript endpoint (and, within one run,
    the same cached ``client``) already used by :func:`resolve_domains`
    for domain-retention classification -- this adds no additional live
    network call within a run that already looked up this gene's domains.
    Returns ``None`` (never raises) when Genome Nexus has no mapping for
    this gene, matching the existing graceful-degradation pattern used
    throughout this module.
    """
    try:
        payload = client.fetch_canonical_transcript(config.gene_symbol)
    except (GenomeNexusGeneNotFound, requests.RequestException):
        return None
    canonical = parse_canonical_transcript(payload)
    track = gene_track_from_canonical_transcript(canonical)
    pretty_names = _pretty_domain_names(config)
    for domain in track["domains"]:
        domain["name"] = pretty_names.get(domain["accession"], domain["name"])
    return track


def _intragenic_deletion_records(
    raw: pd.DataFrame,
    config: GeneConfig,
    client: GenomeNexusClient,
) -> tuple[list[dict], list[str]]:
    """Same-gene (``Site1_gene == Site2_gene == target``) intragenic
    deletion-style SV records, analogous to a paper's panel-C schematic.

    Distinct from the cross-gene fusion pipeline above: these records
    never pass ``_is_target_protein_fusion`` (no partner gene, no "fusion"
    in ``Event_Info``) and are visualization-only -- they are never fed to
    ``normalize()``/``map_event()``/any registered algorithm, so adding
    this cannot change any existing statistic. Matched generically by the
    "Deletion of N exons[ ]: in frame|out of frame" annotation convention
    already used across cBioPortal/MSK-IMPACT SV data (not specific to any
    one gene).
    """
    target = config.gene_symbol.upper()
    records: list[dict] = []
    warnings: list[str] = []
    same_gene = raw[
        (raw["Site1_Hugo_Symbol"].astype(str).str.upper() == target)
        & (raw["Site2_Hugo_Symbol"].astype(str).str.upper() == target)
    ]
    if same_gene.empty:
        return records, warnings

    # Both breakpoints belong to this same gene/transcript, so there is no
    # Five_prime_gene/Three_prime_gene partner annotation to draw a role
    # from (unlike the cross-gene fusion path below). But the deletion is
    # still directional: whichever breakpoint sits upstream of the other in
    # transcript order is the retained-up-to (five_prime-role) side, and the
    # other is the resumed-from (three_prime-role) side -- determined from
    # strand, not any hardcoded gene fact. Best-effort: if the strand can't
    # be resolved, both breakpoints fall back to the direction-unaware
    # nearest-exon behavior, same as before this fix.
    strand: int | None = None
    try:
        canonical_for_strand = parse_canonical_transcript(
            client.fetch_canonical_transcript(config.gene_symbol)
        )
        if canonical_for_strand.exons:
            strand = canonical_for_strand.exons[0].strand
    except (GenomeNexusGeneNotFound, requests.RequestException):
        strand = None

    for _, row in same_gene.iterrows():
        event_info = str(row.get("Event_Info") or "")
        match = _DELETION_EVENT_INFO_PATTERN.search(event_info)
        if match is None:
            continue
        try:
            genomic_positions = (int(row.get("Site1_Position")), int(row.get("Site2_Position")))
            positions_aa = []
            for index, genomic_position in enumerate(genomic_positions):
                role = None
                if strand is not None:
                    other_position = genomic_positions[1 - index]
                    is_upstream = (
                        genomic_position <= other_position
                        if strand == 1
                        else genomic_position >= other_position
                    )
                    role = "five_prime" if is_upstream else "three_prime"
                mapping = resolve_breakpoint_protein_position(
                    None,
                    config,
                    breakpoint_genomic=genomic_position,
                    genome_nexus_client=client,
                    role=role,
                )
                if mapping.breakpoint_protein_position is None:
                    raise ValueError("Genome Nexus returned no protein position")
                positions_aa.append(mapping.breakpoint_protein_position)
        except Exception as exc:  # noqa: BLE001 - best-effort, visualization-only
            warnings.append(
                f"Skipped intragenic-deletion record for sample "
                f"{row.get('Sample_Id')}: {type(exc).__name__}: {exc}"
            )
            continue
        retained_up_to_aa, resumed_from_aa = sorted(positions_aa)
        records.append(
            {
                "sample_id": row.get("Sample_Id"),
                "annotation": row.get("Annotation"),
                "event_info": event_info,
                "n_exons_deleted": int(match.group(1)),
                "frame_status": match.group(2).lower().replace(" ", "-"),
                "retained_up_to_aa": retained_up_to_aa,
                "resumed_from_aa": resumed_from_aa,
            }
        )
    return records, warnings


def _load_benchmark_config(gene_symbol: str) -> GeneConfig:
    try:
        config = load_gene_config(gene_symbol)
    except FileNotFoundError as exc:
        raise RealBenchmarkInputError(
            f"Unknown gene {gene_symbol!r}. Run `cfh list-genes` to see configured genes."
        ) from exc
    if config.gene_symbol is None or config.entrez_gene_id is None:
        raise RealBenchmarkInputError(
            f"Gene {gene_symbol!r} needs gene_symbol and entrez_gene_id in its config."
        )
    if not config.key_domains:
        raise RealBenchmarkInputError(
            f"Gene {gene_symbol!r} has no key domain configured for domain-retention analysis."
        )
    return config


def _maybe_load_gene_pair_config(gene_symbol: str) -> GeneConfig | None:
    """Return the loaded ``GeneConfig`` if ``gene_symbol`` resolves to a
    ``gene_pair`` config (e.g. ``EML4-ALK``), else ``None``.

    An unknown gene symbol is left to :func:`_load_benchmark_config`'s
    existing, actionable ``RealBenchmarkInputError`` -- this only ever
    short-circuits a config that loads successfully and opts into
    ``gene_pair``, so a single-gene lookup failure surfaces the same error
    it always has.
    """
    try:
        config = load_gene_config(gene_symbol)
    except FileNotFoundError:
        return None
    return config if config.gene_pair is not None else None


def _partner_component_configs(pair_config: GeneConfig) -> list[GeneConfig]:
    """Curated single-gene configs among ``pair_config.gene_pair`` that can
    be live-fetched (i.e. have their own curated YAML with a
    ``gene_symbol`` and ``entrez_gene_id``).

    A pair member with no curated single-gene config (e.g. EML4, which has
    no standalone YAML in this registry) is simply not fetched -- this
    mirrors PR #62's manual EML4-ALK workaround, which tested the
    configured pair against only the live-fetched events of the pair
    member that already had a curated config (ALK).
    """
    if pair_config.gene_pair is None:
        raise RealBenchmarkInputError("gene_pair config has no configured gene_pair")
    configs = []
    for symbol in pair_config.gene_pair:
        try:
            candidate = load_gene_config(symbol)
        except FileNotFoundError:
            continue
        if candidate.gene_symbol is not None and candidate.entrez_gene_id is not None:
            configs.append(candidate)
    if not configs:
        raise RealBenchmarkInputError(
            f"Gene pair {pair_config.gene_pair!r} has no partner with a curated "
            "single-gene config (gene_symbol + entrez_gene_id) to live-fetch structural "
            "variants from. Add one under src/cfh/genes/configs/ for at least one partner."
        )
    return configs


def run_gene_pair_benchmark(
    pair_config: GeneConfig,
    study_id: str,
    *,
    n_permutations: int = 1_000,
    algorithm_params: dict[str, dict] | None = None,
) -> RealBenchmarkRun:
    """Live joint-partner benchmark for a ``gene_pair`` config (e.g. EML4-ALK).

    Before this, a ``gene_pair`` config could not run through
    ``cfh analyze``/``cfh real-benchmark`` at all: the single-gene
    ``_load_benchmark_config`` path requires ``gene_symbol``/
    ``entrez_gene_id``/``key_domains``, none of which a ``gene_pair``
    config declares by design (see PR #62). This is the dedicated
    gene-pair path those CLI commands now route to instead of erroring.

    Reuses the existing single-gene ingestion/normalization pipeline --
    :func:`run_real_benchmark` itself -- once per curated partner gene (see
    :func:`_partner_component_configs`), pools their already-mapped
    ``events`` (deduplicated by ``Event_id``, since a fusion between two
    curated partners would otherwise be fetched twice), and tests the
    configured ordered pair for enrichment via
    :class:`~cfh.algorithms.joint_partner.JointPartnerMode`, which
    deliberately only needs gene-pair identities, not domain/breakpoint
    data. No breakpoint mapping, domain classification, or PDF report is
    produced here -- those are single-gene-domain concepts that don't apply
    to a pair-enrichment result.
    """
    if pair_config.gene_pair is None:
        raise RealBenchmarkInputError("run_gene_pair_benchmark requires a gene_pair GeneConfig")
    gene5, gene3 = pair_config.gene_pair
    component_configs = _partner_component_configs(pair_config)

    pooled_events: dict[str, FusionEvent] = {}
    pooled_rows: list[dict] = []
    warnings: list[str] = []
    endpoints: list[str] = []
    raw_count_total = 0
    profile_ids: list[str] = []
    for component in component_configs:
        component_run = run_real_benchmark(
            component.gene_symbol, study_id, n_permutations=n_permutations
        )
        for event in component_run.events:
            pooled_events[event.Event_id] = event
        pooled_rows.extend(component_run.rows)
        warnings.extend(
            f"[{component.gene_symbol}] {warning}" for warning in component_run.warnings
        )
        endpoints.extend(
            endpoint for endpoint in component_run.endpoints if endpoint not in endpoints
        )
        raw_count_total += component_run.raw_structural_variant_count
        profile_ids.append(component_run.molecular_profile_id)

    events = list(pooled_events.values())
    joint_partner_params = {
        "joint_partner": dict((algorithm_params or {}).get("joint_partner", {}))
    }
    joint_result = run_algorithms(["joint_partner"], events, [], pair_config, joint_partner_params)[
        0
    ]
    pair_results = (joint_result.Tables or {}).get("pair_results") or []
    pair_stats = pair_results[0] if pair_results else {}
    if joint_result.Warnings:
        warnings.extend(joint_result.Warnings)

    summary = {
        "gene_pair": [gene5, gene3],
        "component_genes": [component.gene_symbol for component in component_configs],
        "raw_structural_variant_count": raw_count_total,
        "total_fusions": len(events),
        "mapped_fusions": len(events),
        "in_frame_count": sum(event.Frame_status == "in-frame" for event in events),
        "eligible_event_count": pair_stats.get("eligible_event_count", 0),
        "observed_count": pair_stats.get("observed_count", 0),
        "expected_count": pair_stats.get("expected_count", 0.0),
        "fisher_p_value": pair_stats.get("p_value"),
        "fisher_odds_ratio": pair_stats.get("odds_ratio"),
        "is_enriched": joint_result.Summary.get("is_enriched"),
    }
    if pair_config.mechanism_note:
        summary["mechanism_note"] = pair_config.mechanism_note
    return RealBenchmarkRun(
        gene_symbol=f"{gene5}-{gene3}",
        study_id=study_id,
        molecular_profile_id=",".join(dict.fromkeys(profile_ids)),
        retrieved_at=datetime.now(timezone.utc),
        raw_structural_variant_count=raw_count_total,
        events=events,
        features=[],
        rows=pooled_rows,
        results=[joint_result],
        summary=summary,
        warnings=warnings,
        endpoints=endpoints,
        is_gene_pair=True,
    )


def _unavailable_domain_result(
    message: str,
    events: list[FusionEvent],
    features: list[FusionFeature],
    config: GeneConfig,
    n_permutations: int,
) -> AlgorithmResult:
    return AlgorithmResult(
        Algorithm="domain_retention",
        Algorithm_version="0.1.0",
        Parameters={"seed": 42, "n_permutations": n_permutations},
        Summary={
            "fisher_odds_ratio": None,
            "fisher_p_value": None,
            "permutation_empirical_p_value": None,
            "observed_in_frame_retention_rate": None,
        },
        Tables={
            "frame_domain_contingency_table": build_frame_domain_contingency_table(
                events, features, config
            ),
            "permutation_null_retention_rates": [],
        },
        Warnings=[message],
    )


def _no_key_domain_result(gene_config: GeneConfig) -> AlgorithmResult:
    """No-op ``domain_retention`` result for a gene with no configured key
    domain at all (e.g. an auto-generated config for a gene whose canonical
    transcript carries no annotated Pfam domain). Distinct from
    :func:`_unavailable_domain_result` (which still has a target domain,
    just no in-frame mapped observation of it): here there is no target
    domain to build a contingency table against in the first place, so
    ``build_frame_domain_contingency_table`` is never called (it would
    raise for an empty ``key_domains``). This is the same graceful-skip
    shape ``domain_disruption`` already returns when
    ``disruption_required_domains`` is unset.
    """
    return AlgorithmResult(
        Algorithm="domain_retention",
        Algorithm_version="0.1.0",
        Parameters={},
        Summary={
            "fisher_odds_ratio": None,
            "fisher_p_value": None,
            "permutation_empirical_p_value": None,
            "observed_in_frame_retention_rate": None,
        },
        Tables={
            "frame_domain_contingency_table": [[0, 0], [0, 0]],
            "permutation_null_retention_rates": [],
        },
        Warnings=[
            f"{gene_config.gene_symbol or gene_config.gene_pair} has no key_domains "
            "configured; domain-retention analysis was skipped."
        ],
    )


def analyze_structural_variant_calls_with_config(
    calls: list[dict],
    config: GeneConfig,
    study_id: str,
    *,
    molecular_profile_id: str | None = None,
    clinical_df: pd.DataFrame | None = None,
    genome_nexus_client: GenomeNexusClient | None = None,
    n_permutations: int = 1_000,
    algorithm_names: list[str] | None = None,
    algorithm_params: dict[str, dict] | None = None,
    extra_warnings: list[str] | None = None,
) -> RealBenchmarkRun:
    """Normalize and analyze already-fetched cBioPortal SV API objects
    against an already-resolved ``GeneConfig``.

    This is the config-agnostic core :func:`analyze_structural_variant_calls`
    delegates to after resolving a curated config by gene symbol. Callers
    that already have a ``GeneConfig`` in hand -- e.g. a genome-wide scan
    using an auto-generated config for a gene with no curated YAML file --
    call this directly instead, so they are not forced to write one to disk
    first. Unlike the curated lookup path, ``config.key_domains`` being
    empty is not an error here: it degrades to a graceful no-op
    ``domain_retention`` result (see :func:`_no_key_domain_result`), the
    same opt-in/no-op pattern already used for ``disruption_required_domains``,
    ``expected_retained_exon_hint``, and ``gene_pair``.

    ``algorithm_params`` lets a caller pass through additional per-algorithm
    parameters (e.g. adaptive-permutation knobs, or the
    ``mutation_cooccurrence`` live-fetch data assembled by
    :func:`run_real_benchmark`) merged under each algorithm's existing
    defaults below. ``extra_warnings`` lets a caller (e.g. a live-fetch step
    that happened before this function was called) seed the returned run's
    ``warnings`` with messages of its own -- both default to ``None``/empty
    and change nothing for a caller that doesn't pass them.
    """
    if n_permutations <= 0:
        raise RealBenchmarkInputError("n_permutations must be positive")
    if config.gene_symbol is None:
        raise RealBenchmarkInputError(
            "analyze_structural_variant_calls_with_config requires a single-gene "
            "GeneConfig (gene_symbol set, not gene_pair)"
        )
    algorithm_params = algorithm_params or {}
    profile_id = molecular_profile_id or f"{study_id}_structural_variants"
    client = genome_nexus_client or GenomeNexusClient()

    raw = cbioportal_api.structural_variants_to_dataframe(calls)
    # Restrict enrichment to annotations: patient IDs/panel metadata must not
    # change the existing statistical inputs or grouping.
    clinical_annotations = None
    if clinical_df is not None:
        clinical_annotations = clinical_df.reindex(
            columns=["Sample_id", "Tumor_type", "Oncotree_code"]
        ).astype(object)
        clinical_annotations = clinical_annotations.fillna("")
    normalized = normalize(raw, clinical_annotations, study_id)
    selected = [
        (row.to_dict(), event)
        for (_, row), event in zip(raw.iterrows(), normalized, strict=True)
        if _is_target_protein_fusion(event, config.gene_symbol)
    ]

    warnings: list[str] = list(extra_warnings or [])
    target_canonical = None
    needs_domain_lookup = bool(config.key_domains or config.disruption_required_domains)
    if selected and needs_domain_lookup:
        try:
            domains = resolve_domains(
                config.gene_symbol,
                config.protein_id,
                genome_nexus_client=client,
            )
        except requests.RequestException as exc:
            raise RealBenchmarkNetworkError(
                f"Genome Nexus/domain lookup failed for {config.gene_symbol}: {exc}. "
                "Check network access and https://www.genomenexus.org availability, then retry."
            ) from exc
    elif not selected:
        domains = []
        warnings.append(
            f"No protein-fusion records for {config.gene_symbol} were returned by "
            f"{profile_id}. Verify the gene and study ID, or try a study with SV data."
        )
    else:
        # No key_domains/disruption_required_domains configured at all, so no
        # domain would ever be classified against these results regardless
        # of what a domain lookup returned (see _combined_domains) -- skip
        # the lookup (and its network call) entirely rather than resolve
        # domains nothing will use.
        domains = []
    if selected:
        try:
            target_canonical = parse_canonical_transcript(
                client.fetch_canonical_transcript(config.gene_symbol)
            )
            config = derive_gene_config_defaults(config, target_canonical)
        except requests.RequestException as exc:
            raise RealBenchmarkNetworkError(
                f"Genome Nexus target-locus lookup failed for {config.gene_symbol}: {exc}. "
                "Check network access and https://www.genomenexus.org availability, then retry."
            ) from exc
    if target_canonical is not None:
        target_locus = _target_locus(target_canonical)
        target_cds_min_genomic, target_cds_max_genomic = cds_bounds_from_utrs(target_canonical.utrs)
    else:
        target_locus = None
        target_cds_min_genomic = None
        target_cds_max_genomic = None
    domain_source = _ResolvedDomainSource(domains)

    gene_track: dict | None = None
    intragenic_deletions: list[dict] = []
    if len(raw):
        gene_track = _build_gene_track(config, client)
        deletion_records, deletion_warnings = _intragenic_deletion_records(raw, config, client)
        intragenic_deletions = deletion_records
        warnings.extend(deletion_warnings)

    events: list[FusionEvent] = []
    features: list[FusionFeature] = []
    rows: list[dict] = []
    mapping_sensitivity: dict[str, bool | None] = {}
    has_key_domain = bool(config.key_domains)
    target_key = (
        (config.key_domains[0].key or config.key_domains[0].name) if has_key_domain else None
    )

    for row, event in selected:
        try:
            role = _target_role(event, config.gene_symbol)
            if target_locus is None:  # pragma: no cover - no selected records means no loop
                raise ValueError("target locus was not resolved")
            breakpoint = _target_breakpoint(row, config.gene_symbol, target_locus)
            mapping = resolve_breakpoint_protein_position(
                None,
                config,
                breakpoint_genomic=breakpoint,
                genome_nexus_client=client,
                role=role,
            )
            if mapping.breakpoint_protein_position is None:
                raise ValueError("Genome Nexus returned no protein position")
            feature = map_event(
                event,
                config,
                role=role,
                junction_position_aa=mapping.breakpoint_protein_position,
                domain_source=domain_source,
            ).model_copy(
                update={
                    "Breakpoint_exon": mapping.breakpoint_exon,
                    "Retained_exons": _retained_exon_ranks(
                        target_canonical,
                        mapping.breakpoint_protein_position,
                        role,
                    ),
                }
            )
        except requests.RequestException as exc:
            raise RealBenchmarkNetworkError(
                f"Genome Nexus breakpoint mapping failed for {config.gene_symbol}: {exc}. "
                "Check network access and https://www.genomenexus.org availability, then retry."
            ) from exc
        except Exception as exc:
            warnings.append(
                f"Skipped {event.Event_id} ({event.Fusion_name or 'unnamed fusion'}): "
                f"{type(exc).__name__}: {exc}"
            )
            continue

        # Independent QA cross-check (genome-nexus/fusion-annotation's ported
        # arithmetic, fed our own already-correct UTR-derived CDS bounds --
        # see cfh.mapping.fusion_annotation_crosscheck). Isolated in its own
        # try/except so a cross-check failure can never skip a row or affect
        # any existing statistic -- it only ever adds a new, additional field.
        try:
            crosscheck = crosscheck_breakpoint_protein_position(
                target_canonical.exons,
                breakpoint,
                target_canonical.exons[0].strand,
                role,
                mapping.breakpoint_protein_position,
                cds_min_genomic=target_cds_min_genomic,
                cds_max_genomic=target_cds_max_genomic,
            )
        except Exception as exc:  # noqa: BLE001 - QA cross-check must never fail the row
            crosscheck = FusionAnnotationCrosscheckResult(
                protein_position=None,
                cds_position=None,
                exon_rank=None,
                is_hybrid_codon=None,
                agrees=None,
                note=None,
                error=(
                    f"fusion-annotation cross-check raised unexpectedly: "
                    f"{type(exc).__name__}: {exc}"
                ),
            )
        if crosscheck.agrees is False:
            warnings.append(
                f"Fusion-annotation cross-check disagreement for {event.Event_id} "
                f"({config.gene_symbol}, {role}): ours={mapping.breakpoint_protein_position}, "
                f"fusion-annotation={crosscheck.protein_position}. {crosscheck.note}"
            )

        events.append(event)
        features.append(feature)
        mapping_sensitivity[event.Event_id] = mapping.is_intronic_breakpoint
        domain_detail = (
            (feature.Domain_retention_details or {}).get(target_key) if target_key else None
        )
        rows.append(
            {
                "event_id": event.Event_id,
                "sample_id": event.Sample_id,
                "patient_id": event.Patient_id,
                "tumor_type": event.Tumor_type,
                "oncotree_code": event.Oncotree_code,
                "fusion_name": event.Fusion_name,
                "partner_gene": _partner(event, config.gene_symbol),
                "frame_status": event.Frame_status,
                "target_role": role,
                "breakpoint_genomic": breakpoint,
                "breakpoint_exon": mapping.breakpoint_exon,
                "breakpoint_protein_position": mapping.breakpoint_protein_position,
                "is_intronic_breakpoint": mapping.is_intronic_breakpoint,
                "fusion_annotation_crosscheck_protein_position": crosscheck.protein_position,
                "fusion_annotation_crosscheck_agrees": crosscheck.agrees,
                "fusion_annotation_crosscheck_detail": (
                    crosscheck.note if crosscheck.agrees is not None else crosscheck.error
                ),
                "domain_status": (
                    (feature.Domain_retention_flags or {}).get(target_key, "unknown")
                    if target_key
                    else "unknown"
                ),
                "domain_retained_fraction": (
                    domain_detail.Retained_fraction if domain_detail else None
                ),
                "domain_is_truncated": domain_detail.Is_truncated if domain_detail else None,
                "retained_domains": "; ".join(feature.Retained_domains or []),
                "lost_domains": "; ".join(feature.Lost_domains or []),
                "disrupted_domains": "; ".join(feature.Disrupted_domains or []),
                "source_annotation_text": " | ".join(
                    value
                    for value in (
                        str(row.get("Annotation") or ""),
                        str(row.get("Event_Info") or ""),
                    )
                    if value
                ),
                "source_site2_effect_on_frame": row.get("Site2_Effect_On_Frame"),
            }
        )

    in_frame_mapped = has_key_domain and any(
        event.Frame_status == "in-frame"
        and (feature.Domain_retention_flags or {}).get(target_key)
        in {"retained", "lost", "disrupted"}
        for event, feature in zip(events, features, strict=True)
    )
    if not has_key_domain:
        domain_result = _no_key_domain_result(config)
        warnings.append(domain_result.Warnings[0])
    elif not in_frame_mapped:
        message = (
            "Domain-retention statistics are unavailable because no mapped in-frame "
            "protein-fusion record has a known domain state. Verify the gene/study IDs "
            "and source frame annotations."
        )
        warnings.append(message)
        domain_result = _unavailable_domain_result(
            message, events, features, config, n_permutations
        )
    else:
        try:
            domain_result = run_algorithms(
                ["domain_retention"],
                events,
                features,
                config,
                {
                    "domain_retention": {
                        "seed": 42,
                        "n_permutations": n_permutations,
                        "genome_nexus_client": client,
                        **algorithm_params.get("domain_retention", {}),
                    }
                },
            )[0]
        except requests.RequestException as exc:
            raise RealBenchmarkNetworkError(
                f"Genome Nexus permutation mapping failed for {config.gene_symbol}: {exc}. "
                "Check network access and https://www.genomenexus.org availability, then retry."
            ) from exc
        except ValueError as exc:
            message = f"Domain-retention statistics are unavailable: {exc}"
            warnings.append(message)
            domain_result = _unavailable_domain_result(
                message, events, features, config, n_permutations
            )
    selected_events = [event for _, event in selected]
    requested_algorithms = algorithm_names or ["domain_retention", "frequency"]
    other_algorithms = [name for name in requested_algorithms if name != "domain_retention"]
    other_results = run_algorithms(
        other_algorithms,
        selected_events,
        features,
        config,
        {
            "confidence_stats": {
                **resolve_confidence_stats_params(config, algorithm_params.get("confidence_stats")),
            },
            "frequency": {"dedup_by_patient": False, **algorithm_params.get("frequency", {})},
            "cutpoint_detection": {
                "n_permutations": n_permutations,
                "genome_nexus_client": client,
                **algorithm_params.get("cutpoint_detection", {}),
            },
            "window_detection": {
                "n_permutations": n_permutations,
                "genome_nexus_client": client,
                "mapping_sensitivity": mapping_sensitivity,
                **algorithm_params.get("window_detection", {}),
            },
            **{
                name: value
                for name, value in algorithm_params.items()
                if name
                not in {
                    "confidence_stats",
                    "frequency",
                    "cutpoint_detection",
                    "window_detection",
                    "domain_retention",
                }
            },
        },
        extra_results=[domain_result],
    )
    results_by_name = {result.Algorithm: result for result in [domain_result, *other_results]}
    results = [results_by_name[name] for name in requested_algorithms if name in results_by_name]
    frequency_result = results_by_name.get("frequency")
    if frequency_result is None:
        frequency_result = FrequencyAnalysis().run(
            selected_events, features, config, {"dedup_by_patient": False}
        )
    partner_counts = frequency_result.Tables["Partner_gene_counts"]
    # Frame and retention summaries describe the same post-mapping population.
    mapped_total = len(rows)
    in_frame_count = sum(row["frame_status"] == "in-frame" for row in rows)
    retained_count = sum(row["domain_status"] == "retained" for row in rows)
    in_frame_retained_count = sum(
        row["frame_status"] == "in-frame" and row["domain_status"] == "retained" for row in rows
    )
    domain_definition = (
        next(
            (
                domain
                for domain in domains
                if domain.accession == config.key_domains[0].accession
                or domain.name == config.key_domains[0].accession
            ),
            None,
        )
        if has_key_domain
        else None
    )
    # Every one of config.key_domains that Genome Nexus resolved a span for
    # -- generalizes domain_definition above (which only ever looks at
    # key_domains[0]) so a gene configured with more than one
    # retention-target domain gets all of them surfaced in summary
    # ["key_domains"] for the domain-retention-track visualization
    # (_domain_track_svg), not just the first.
    key_domain_definitions = [
        {
            "name": key_domain.name,
            "accession": key_domain.accession,
            "start_aa": matched.start_aa,
            "end_aa": matched.end_aa,
        }
        for key_domain in config.key_domains
        for matched in [
            next(
                (
                    domain
                    for domain in domains
                    if domain.accession == key_domain.accession
                    or domain.name == key_domain.accession
                ),
                None,
            )
        ]
        if matched is not None
    ]
    configured_key_domains = [domain.model_dump(mode="json") for domain in config.key_domains]
    configured_disruption_domains = [
        domain.model_dump(mode="json") for domain in config.disruption_required_domains
    ]
    total = len(selected_events)
    crosscheck_agree_count = sum(
        1 for row in rows if row["fusion_annotation_crosscheck_agrees"] is True
    )
    crosscheck_disagree_count = sum(
        1 for row in rows if row["fusion_annotation_crosscheck_agrees"] is False
    )
    crosscheck_unavailable_count = sum(
        1 for row in rows if row["fusion_annotation_crosscheck_agrees"] is None
    )
    summary = {
        "raw_structural_variant_count": len(calls),
        "total_fusions": total,
        "mapped_fusions": len(features),
        "skipped_fusions": total - len(features),
        "in_frame_count": in_frame_count,
        "in_frame_percent": 100 * in_frame_count / mapped_total if mapped_total else 0.0,
        "kinase_retained_count": retained_count,
        "kinase_retained_percent": 100 * retained_count / mapped_total if mapped_total else 0.0,
        "in_frame_kinase_retained_count": in_frame_retained_count,
        "fisher_odds_ratio": domain_result.Summary["fisher_odds_ratio"],
        "fisher_p_value": domain_result.Summary["fisher_p_value"],
        "permutation_p_value": domain_result.Summary["permutation_empirical_p_value"],
        "frame_domain_contingency_table": domain_result.Tables["frame_domain_contingency_table"],
        "partner_counts": partner_counts,
        "domain_accession": config.key_domains[0].accession if has_key_domain else None,
        "domain_start_aa": domain_definition.start_aa if domain_definition else None,
        "domain_end_aa": domain_definition.end_aa if domain_definition else None,
        "key_domains": key_domain_definitions,
        "fusion_annotation_crosscheck": {
            "agree": crosscheck_agree_count,
            "disagree": crosscheck_disagree_count,
            "unavailable": crosscheck_unavailable_count,
        },
        "configured_key_domains": configured_key_domains,
        "configured_disruption_required_domains": configured_disruption_domains,
    }
    return RealBenchmarkRun(
        gene_symbol=config.gene_symbol,
        study_id=study_id,
        molecular_profile_id=profile_id,
        retrieved_at=datetime.now(timezone.utc),
        raw_structural_variant_count=len(calls),
        events=events,
        features=features,
        rows=rows,
        results=results,
        summary=summary,
        warnings=warnings,
        endpoints=[
            "https://www.cbioportal.org/api/structural-variant/fetch",
            f"{getattr(client, 'base_url', 'https://www.genomenexus.org')}"
            f"/ensembl/canonical-transcript/hgnc/{config.gene_symbol}",
        ],
        reference=(config.benchmark_reference.model_dump() if config.benchmark_reference else None),
        gene_track=gene_track,
        intragenic_deletions=intragenic_deletions,
    )


def analyze_structural_variant_calls(
    calls: list[dict],
    gene_symbol: str,
    study_id: str,
    *,
    molecular_profile_id: str | None = None,
    clinical_df: pd.DataFrame | None = None,
    genome_nexus_client: GenomeNexusClient | None = None,
    n_permutations: int = 1_000,
    algorithm_names: list[str] | None = None,
    algorithm_params: dict[str, dict] | None = None,
    extra_warnings: list[str] | None = None,
) -> RealBenchmarkRun:
    """Normalize and analyze already-fetched cBioPortal SV API objects.

    Resolves ``gene_symbol`` against the curated
    :func:`~cfh.genes.registry.load_gene_config` registry (requiring
    ``entrez_gene_id`` and at least one configured key domain, as before)
    and delegates to :func:`analyze_structural_variant_calls_with_config`.
    """
    config = _load_benchmark_config(gene_symbol)
    return analyze_structural_variant_calls_with_config(
        calls,
        config,
        study_id,
        molecular_profile_id=molecular_profile_id,
        clinical_df=clinical_df,
        genome_nexus_client=genome_nexus_client,
        n_permutations=n_permutations,
        algorithm_names=algorithm_names,
        algorithm_params=algorithm_params,
        extra_warnings=extra_warnings,
    )


def _fetch_mutual_exclusivity_params(
    config: GeneConfig, study_id: str, study_config: StudyConfig | None
) -> tuple[dict[str, dict] | None, list[str]]:
    """Best-effort live fetch of the comparator alteration + cohort-sample-
    universe data ``mutation_cooccurrence`` needs, for a gene that opts in
    via ``GeneConfig.mutual_exclusivity_targets``.

    Never raises: a network or lookup failure here only produces a warning
    and (for a target it applies to) fewer comparator_alterations -- this
    evidence layer is additive, and must never take down an otherwise
    successful domain-retention benchmark run over a failure fetching its
    own, separate data.
    """
    warnings: list[str] = []
    sample_list_id = (
        study_config.all_sample_list_id(study_id) if study_config else f"{study_id}_all"
    )
    try:
        cohort_sample_ids = cbioportal_api.fetch_sample_list_ids(sample_list_id)
    except requests.RequestException as exc:
        warnings.append(
            f"Could not fetch cohort sample universe {sample_list_id!r} for "
            f"mutation_cooccurrence: {type(exc).__name__}: {exc}. Co-occurrence/"
            "mutual-exclusivity analysis was skipped for this run."
        )
        return None, warnings

    comparator_alterations: list[dict] = []
    for target in config.mutual_exclusivity_targets:
        if target.entrez_gene_id is None:
            warnings.append(
                f"mutual_exclusivity_targets entry for {target.gene} has no "
                "entrez_gene_id configured; it cannot be live-fetched and was skipped."
            )
            continue
        try:
            if target.alteration_type == "point_mutation":
                mutation_profile_id = (
                    study_config.mutation_profile_id(study_id)
                    if study_config
                    else f"{study_id}_mutations"
                )
                calls = cbioportal_api.fetch_mutations(
                    [target.entrez_gene_id], [mutation_profile_id]
                )
                events, row_warnings = normalize_mutations(calls, target.gene, study_id)
            elif target.alteration_type.startswith("cna_"):
                cna_profile_id = (
                    study_config.discrete_cna_profile_id(study_id)
                    if study_config
                    else f"{study_id}_cna"
                )
                calls = cbioportal_api.fetch_discrete_copy_number(
                    [target.entrez_gene_id], cna_profile_id, sample_list_id
                )
                events, row_warnings = normalize_discrete_copy_number(calls, target.gene, study_id)
            else:
                warnings.append(
                    "Unrecognized mutual_exclusivity_targets alteration_type "
                    f"{target.alteration_type!r} for {target.gene}; skipped."
                )
                continue
        except requests.RequestException as exc:
            warnings.append(
                f"Could not fetch {target.alteration_type} data for {target.gene} "
                f"(mutation_cooccurrence comparator): {type(exc).__name__}: {exc}."
            )
            continue
        warnings.extend(row_warnings)
        comparator_alterations.extend(event.model_dump() for event in events)

    return (
        {
            "mutation_cooccurrence": {
                "cohort_sample_ids": cohort_sample_ids,
                "comparator_alterations": comparator_alterations,
            }
        },
        warnings,
    )


def _fetch_expression_association_params(
    config: GeneConfig, study_config: StudyConfig | None, study_id: str
) -> tuple[dict[str, Any], str | None]:
    """Best-effort live mRNA-expression fetch feeding the
    ``expression_association`` algorithm.

    Returns an empty params dict (never a crash) when this cohort has no
    configured mRNA-expression molecular profile at all -- e.g. a targeted
    DNA panel like ``msk_impact_50k_2026``, which carries no
    ``MRNA_EXPRESSION`` profile in the first place, so
    ``load_study_config`` returns ``None`` for it and this never even
    attempts a request. A transient network failure fetching an
    otherwise-configured profile degrades the same way, surfaced as a
    warning rather than a raised exception, since expression-association
    evidence is optional corroborating evidence, not required for the run.
    """
    if study_config is None:
        return {}, None
    profile_id = study_config.mrna_expression_profile_id(study_id)
    if profile_id is None:
        return {}, None
    try:
        records = cbioportal_api.fetch_molecular_data(
            [config.entrez_gene_id],
            profile_id,
            sample_list_id=f"{study_id}_all",
        )
    except requests.RequestException as exc:
        return {}, (
            f"mRNA-expression fetch failed for {config.gene_symbol} in {profile_id}: "
            f"{type(exc).__name__}: {exc}. expression_association analysis was skipped."
        )
    expression_by_sample = cbioportal_api.molecular_data_to_expression_by_sample(records)
    if not expression_by_sample:
        return {}, (
            f"No mRNA-expression records were returned for {config.gene_symbol} from "
            f"{profile_id}; expression_association analysis was skipped."
        )
    return {
        "expression_by_sample": expression_by_sample,
        "cohort_sample_ids": list(expression_by_sample.keys()),
        "expression_field": "mRNA expression z-score",
    }, None


def run_real_benchmark(
    gene_symbol: str,
    study_id: str,
    *,
    n_permutations: int = 1_000,
    algorithm_names: list[str] | None = None,
    algorithm_params: dict[str, dict] | None = None,
) -> RealBenchmarkRun:
    """Fetch and analyze a gene's structural variants from cBioPortal.

    When the resolved gene config opts into ``mutual_exclusivity_targets``,
    this also live-fetches the comparator alteration(s) and the cohort's
    full sample universe (see :func:`_fetch_mutual_exclusivity_params`) and
    passes them through to the ``mutation_cooccurrence`` algorithm -- a
    gene that doesn't configure this makes no additional network calls and
    is completely unaffected.

    A ``gene_symbol`` that resolves to a ``gene_pair`` config (e.g.
    ``EML4-ALK``) is routed to :func:`run_gene_pair_benchmark` instead --
    see that function's docstring for why the single-gene path below
    cannot handle it.
    """
    pair_config = _maybe_load_gene_pair_config(gene_symbol)
    if pair_config is not None:
        return run_gene_pair_benchmark(
            pair_config,
            study_id,
            n_permutations=n_permutations,
            algorithm_params=algorithm_params,
        )
    config = _load_benchmark_config(gene_symbol)
    study_config = load_study_config(study_id)
    profile_id = (
        study_config.molecular_profile_id(study_id)
        if study_config
        else f"{study_id}_structural_variants"
    )
    genome_nexus_client = GenomeNexusClient(
        base_url=(
            study_config.genome_nexus_base_url if study_config else "https://www.genomenexus.org"
        )
    )
    caller_algorithm_params: dict[str, dict] = {
        name: dict(value) for name, value in (algorithm_params or {}).items()
    }
    extra_warnings: list[str] = []
    if config.mutual_exclusivity_targets:
        mutual_exclusivity_params, fetch_warnings = _fetch_mutual_exclusivity_params(
            config, study_id, study_config
        )
        extra_warnings.extend(fetch_warnings)
        if mutual_exclusivity_params:
            for name, value in mutual_exclusivity_params.items():
                caller_algorithm_params[name] = {
                    **value,
                    **caller_algorithm_params.get(name, {}),
                }
    try:
        calls = cbioportal_api.fetch_structural_variants(
            [config.entrez_gene_id],
            [profile_id],
        )
    except requests.RequestException as exc:
        raise RealBenchmarkNetworkError(
            f"cBioPortal request failed for gene {config.gene_symbol} and profile "
            f"{profile_id}: {exc}. Check the study ID, network access, and "
            "https://www.cbioportal.org availability, then retry."
        ) from exc

    resolved_algorithm_params: dict[str, dict] = dict(caller_algorithm_params)
    expression_params, expression_warning = _fetch_expression_association_params(
        config, study_config, study_id
    )
    if expression_params:
        resolved_algorithm_params["expression_association"] = {
            **expression_params,
            **resolved_algorithm_params.get("expression_association", {}),
        }

    run = analyze_structural_variant_calls(
        calls,
        gene_symbol,
        study_id,
        molecular_profile_id=profile_id,
        clinical_df=cbioportal_api.fetch_sample_tumor_types(
            study_id, [call["sampleId"] for call in calls if call.get("sampleId")]
        ),
        genome_nexus_client=genome_nexus_client,
        n_permutations=n_permutations,
        algorithm_names=algorithm_names,
        algorithm_params=resolved_algorithm_params,
        extra_warnings=extra_warnings,
    )
    if expression_warning:
        run.warnings.append(expression_warning)
    run.endpoints.append(
        f"{cbioportal_api.DEFAULT_BASE_URL}/studies/{study_id}/clinical-data/fetch"
    )
    return run


def run_analysis(
    gene_symbol: str,
    study_id: str,
    *,
    n_permutations: int = 1_000,
) -> RealBenchmarkRun:
    """Run every registered plugin against the shared live mapped input snapshot."""
    return run_real_benchmark(
        gene_symbol,
        study_id,
        n_permutations=n_permutations,
        algorithm_names=list_algorithms(),
    )


def _json_safe(value: object) -> object:
    if isinstance(value, dict):
        return {key: _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return str(value)
    return value


def _format_stat(value: float | int | None) -> str:
    return "unavailable" if value is None or not math.isfinite(value) else f"{value:.6g}"


def markdown_summary(
    run: RealBenchmarkRun,
    *,
    domain_svg_path: str = "visualizations/domain_retention_outliers.svg",
    comparison_svg_path: str = "visualizations/reference_comparison.svg",
    fusion_schematic_svg_path: str | None = None,
    intragenic_deletion_svg_path: str | None = None,
) -> str:
    """Render a concise, checked-in-friendly benchmark report.

    ``fusion_schematic_svg_path``/``intragenic_deletion_svg_path`` are the
    report-relative paths of the SVGs :func:`write_outputs` already wrote
    to ``visualizations/`` (or ``None`` when a schematic wasn't produced,
    e.g. no ``gene_track`` was resolvable, or -- for the deletion schematic
    -- this run has no intragenic-deletion-style records at all); when
    ``None`` the corresponding section is simply omitted, never rendered
    with a broken image link.
    """
    summary = run.summary
    configured_key_domains = summary.get("configured_key_domains") or []
    domain = format_domain_names(configured_domain_names(configured_key_domains))
    domain = domain or summary["domain_accession"] or "configured domain"
    results_by_name = {result.Algorithm: result for result in run.results}
    partners = ", ".join(
        f"{row['Partner_gene']} ({row['Event_count']})" for row in summary["partner_counts"]
    )
    table = summary["frame_domain_contingency_table"]
    lines = [
        f"# {run.gene_symbol} real-data fusion benchmark: {run.study_id}",
        "",
        "Retrieved from public cBioPortal and Genome Nexus on "
        f"{run.retrieved_at.date().isoformat()}.",
        "",
        "## Results",
        "",
        f"- Structural variants returned for {run.gene_symbol}: "
        f"{summary['raw_structural_variant_count']}",
        f"- Protein-fusion records found: {summary['total_fusions']}",
        f"- Protein-fusion records mapped: {summary['mapped_fusions']}",
        f"- Malformed/unmappable fusion records skipped: {summary['skipped_fusions']}",
        f"- In-frame: {summary['in_frame_count']}/{summary['mapped_fusions']} "
        f"({summary['in_frame_percent']:.1f}%)",
        f"- {domain} ({summary['domain_start_aa']}-{summary['domain_end_aa']} aa) retained: "
        f"{summary['kinase_retained_count']}/{summary['mapped_fusions']} "
        f"({summary['kinase_retained_percent']:.1f}%)",
        f"- In-frame and {domain}-retained: "
        f"{summary['in_frame_kinase_retained_count']}/{summary['in_frame_count']}",
        f"- Fisher exact test (one-sided): odds ratio "
        f"{_format_stat(summary['fisher_odds_ratio'])}, "
        f"p={_format_stat(summary['fisher_p_value'])}",
        "- Breakpoint-permutation empirical p-value: "
        f"{_format_stat(summary['permutation_p_value'])}",
        f"- Contingency table `[[retained/in-frame, retained/other], "
        f"[not-retained/in-frame, not-retained/other]]`: `{table}`",
        "",
    ]
    retention_interpretation = domain_interpretation_sentence(
        configured_key_domains,
        fisher_p_value=summary.get("fisher_p_value"),
        fisher_odds_ratio=summary.get("fisher_odds_ratio"),
        effect="retention",
    )
    disruption_summary = (
        results_by_name.get("domain_disruption").Summary
        if results_by_name.get("domain_disruption") is not None
        else {}
    )
    disruption_interpretation = domain_interpretation_sentence(
        summary.get("configured_disruption_required_domains") or [],
        fisher_p_value=disruption_summary.get("fisher_p_value"),
        fisher_odds_ratio=disruption_summary.get("fisher_odds_ratio"),
        effect="disruption",
    )
    lines.extend(
        interpretation
        for interpretation in (retention_interpretation, disruption_interpretation)
        if interpretation
    )
    if retention_interpretation or disruption_interpretation:
        lines.append("")
    cooccurrence_result = results_by_name.get("mutation_cooccurrence")
    cooccurrence_lines = mutual_exclusivity_report_lines(
        cooccurrence_result.model_dump(mode="json") if cooccurrence_result else None,
        run.gene_symbol,
    )
    if cooccurrence_lines:
        lines.extend(cooccurrence_lines)
        lines.append("")
    lines.extend(
        [
            "### Domain retention and discrepancies",
            "",
            f"![Domain retention diagram]({domain_svg_path})",
            "",
            "*Domain-retention positions for analyzed fusion events; red outlines mark "
            "reference discrepancies.*",
            "",
        ]
    )
    if fusion_schematic_svg_path:
        lines.extend(
            [
                "### Fusion-transcript schematic",
                "",
                f"![{run.gene_symbol} fusion-transcript schematic]({fusion_schematic_svg_path})",
                "",
                "*One row per recurrent partner/breakpoint group, sharing one amino-acid "
                f"x-axis for {run.gene_symbol}'s full protein length; the partner-contributed "
                "portion is colored per partner, the retained target-gene portion is "
                "colored by domain-retention status, and a red line marks the breakpoint.*",
                "",
            ]
        )
    if intragenic_deletion_svg_path:
        lines.extend(
            [
                "### Intragenic-deletion schematic",
                "",
                f"![{run.gene_symbol} intragenic-deletion schematic]"
                f"({intragenic_deletion_svg_path})",
                "",
                "*Same-gene (Site1==Site2=="
                f"{run.gene_symbol}) intragenic-deletion-style SV records: a retained "
                "N-terminal block, a plain connector line for the deleted span, and a "
                "resumed C-terminal block.*",
                "",
            ]
        )
    lines.extend(
        [
            "## Method",
            "",
            f"The cBioPortal `{run.molecular_profile_id}` structural-variant profile was "
            "queried by the configured Entrez gene ID. Fusion-annotated records were "
            "adapted to the production SV schema and normalized; when "
            "`site2EffectOnFrame=NA`, frame status was resolved from `Event_Info`, not "
            "copied into `FusionEvent.Frame_status`.",
            "",
            f"{run.gene_symbol} genomic breakpoints were mapped against the Genome Nexus "
            f"canonical transcript, and retention was classified against its returned {domain} "
            "coordinates. Counts are event-level with no patient deduplication. The "
            "Fisher comparison's `other` column combines out-of-frame and unknown-frame "
            "events, as pre-specified by the domain-retention algorithm.",
            "",
            "For each fusion, breakpoint selection preferred the Genome Nexus canonical "
            "transcript's exon-spanned target locus over cBioPortal site labels; malformed "
            "rows with no unambiguous target-locus coordinate were skipped and listed in Warnings.",
            "",
            "## Full-suite highlights",
            "",
            "- Registered algorithms executed: "
            + ", ".join(result.Algorithm for result in run.results),
        ]
    )
    cutpoint = results_by_name.get("cutpoint_detection")
    if cutpoint and cutpoint.Summary.get("determinable"):
        cutpoint_summary = cutpoint.Summary
        cutpoint_exon_label = exon_label_for_protein_position(
            (run.gene_track or {}).get("exon_boundaries_aa"),
            cutpoint_summary["inferred_cutpoint_aa"],
        )
        lines.append(
            "- Cutpoint detection: inferred breakpoint "
            f"{cutpoint_summary['inferred_cutpoint_aa']} aa ({cutpoint_exon_label}); "
            f"corrected permutation p={_format_stat(cutpoint_summary['corrected_p_value'])}."
        )
    elif cutpoint:
        lines.append(
            "- Cutpoint detection: not determinable "
            f"({cutpoint.Summary.get('reason') or 'no reason reported'})."
        )
    composite = results_by_name.get("composite_score")
    ranking = (composite.Tables or {}).get("composite_evidence_ranking", []) if composite else []
    if ranking:
        top = ranking[0]
        lines.append(
            "- Top composite score: "
            f"{top['Partner_gene']} ({top['Event_count']} events), "
            f"{top['Composite_score']:.6g}."
        )
    expression_result = results_by_name.get("expression_association")
    if expression_result is not None:
        expression_sentence = expression_association_sentence(
            run.gene_symbol, expression_result.Summary
        )
        if expression_sentence:
            lines.append(f"- Expression association: {expression_sentence}")
        elif expression_result.Warnings:
            lines.append(f"- Expression association: {expression_result.Warnings[0]}")
    lines.extend(
        [
            "",
            "## Partners",
            "",
            partners or "None",
            "",
        ]
    )
    if run.warnings:
        lines.extend(["## Warnings", ""])
        lines.extend(f"- {warning}" for warning in run.warnings)
        lines.append("")
    if run.reference:
        citation = run.reference["citation"]
        lines.extend(
            [
                "## Reference comparison",
                "",
                f"| Metric | {citation} | This run |",
                "|---|---:|---:|",
                f"| In-frame | {run.reference['in_frame_percent']:.1f}% | "
                f"{summary['in_frame_percent']:.1f}% |",
                f"| Domain retained | {run.reference['domain_retained_percent']:.1f}% | "
                f"{summary['kinase_retained_percent']:.1f}% |",
                "",
                f"![Reference comparison]({comparison_svg_path})",
                "",
                "*Published reference percentages compared with this run.*",
                "",
            ]
        )
    else:
        lines.extend(
            [
                "## Reference comparison",
                "",
                f"![Reference comparison]({comparison_svg_path})",
                "",
                "*Configured reference percentages compared with this run.*",
                "",
            ]
        )
    lines.extend(
        [
            "## Interpretation",
            "",
        ]
    )
    if not summary["total_fusions"]:
        lines.append("No fusion records were available, so comparison is not possible.")
    elif run.gene_symbol.upper() == "BRAF" and run.study_id == "msk_impact_50k_2026":
        lines.extend(
            [
                "This does **not** reproduce the Zehir et al. (PMC5461196) report of "
                "33/33 BRAF fusions being in-frame with the kinase domain retained: "
                f"this live successor cohort has {summary['in_frame_count']}/"
                f"{summary['mapped_fusions']} mapped fusions in-frame and "
                f"{summary['in_frame_kinase_retained_count']}/"
                f"{summary['in_frame_count']} in-frame fusions retaining {domain}.",
                "",
                "`msk_impact_50k_2026` is a newer successor cohort, not the paper's "
                "original `msk_impact_2017` cohort. This is therefore replication in a "
                "related cohort, not a reanalysis of the paper's original 33 cases.",
            ]
        )
    elif run.gene_symbol.upper() == "ALK" and run.study_id == "msk_impact_50k_2026":
        eml4_count = next(
            (
                row["Event_count"]
                for row in summary["partner_counts"]
                if row["Partner_gene"] == "EML4"
            ),
            0,
        )
        lines.append(
            f"EML4 is the recurrent partner ({eml4_count}/{summary['total_fusions']} events), "
            f"and PF07714 retention is {summary['kinase_retained_percent']:.1f}%. These "
            "directions are consistent with the well-known EML4-ALK fusion pattern; this "
            "is a cohort-specific live measurement, not a comparison forced to a literature value."
        )
    elif run.gene_symbol.upper() == "NTRK1" and run.study_id == "msk_impact_50k_2026":
        partner_counts = {
            row["Partner_gene"]: row["Event_count"] for row in summary["partner_counts"]
        }
        lines.append(
            "LMNA and TPM3 each occur in "
            f"{partner_counts.get('LMNA', 0)}/{summary['total_fusions']} and "
            f"{partner_counts.get('TPM3', 0)}/{summary['total_fusions']} events, respectively; "
            f"{summary['in_frame_kinase_retained_count']}/{summary['in_frame_count']} in-frame "
            "events retain PF07714. This is directionally consistent with LMNA-NTRK1/TPM3-NTRK1 "
            "biology, while the all-event in-frame percentage is reported as observed rather than "
            "treated as a literature replication."
        )
    else:
        lines.append("These values describe the live study named above.")
    lines.append("")
    return "\n".join(lines)


def _git_sha() -> str | None:
    try:
        return subprocess.run(
            ["git", "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def _run_id(run: RealBenchmarkRun) -> str:
    gene = re.sub(r"[^a-z0-9]+", "-", run.gene_symbol.lower()).strip("-")
    study = re.sub(r"[^a-z0-9]+", "-", run.study_id.lower()).strip("-")
    timestamp = run.retrieved_at.astimezone(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    return f"{gene}_{study}_{timestamp}"


def _write_tsv(path: Path, rows: list[dict], fieldnames: list[str] | None = None) -> None:
    names = fieldnames or (list(rows[0]) if rows else [])
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=names, delimiter="\t", lineterminator="\n")
        if names:
            writer.writeheader()
            writer.writerows(rows)


def _discrepancies(run: RealBenchmarkRun) -> list[dict]:
    discrepancies = []
    for row in run.rows:
        common = {
            "event_id": row.get("event_id"),
            "partner_gene": row.get("partner_gene"),
            "frame_status": row.get("frame_status"),
            "breakpoint_protein_position": row.get("breakpoint_protein_position"),
            "breakpoint_exon": row.get("breakpoint_exon"),
            "is_intronic_breakpoint": row.get("is_intronic_breakpoint"),
            "retained_domains": row.get("retained_domains", ""),
            "lost_domains": row.get("lost_domains", ""),
            "disrupted_domains": row.get("disrupted_domains", ""),
            "source_annotation_text": row.get("source_annotation_text", ""),
            "source_site2_effect_on_frame": row.get("source_site2_effect_on_frame"),
        }
        if run.reference and (
            row.get("frame_status") != "in-frame" or row.get("domain_status") != "retained"
        ):
            discrepancies.append({"discrepancy_type": "reference_discrepancy", **common})
        source_calls = {
            status
            for status in (
                _source_frame_status(row.get("source_site2_effect_on_frame")),
                _source_frame_status(row.get("source_annotation_text")),
            )
            if status in {"in-frame", "out-of-frame"}
        }
        if any(source_call != row.get("frame_status") for source_call in source_calls):
            discrepancies.append({"discrepancy_type": "source_vs_derived_qa_mismatch", **common})
    return discrepancies


def _domain_track_key_domains(run: RealBenchmarkRun) -> list[dict]:
    """Every key (retention-target) domain with a resolved amino-acid span
    to highlight on the domain-retention track, gene-agnostic and not
    limited to one domain.

    Prefers ``run.summary["key_domains"]`` (all of a gene's configured
    ``key_domains``, see :func:`analyze_structural_variant_calls`). Falls
    back to deriving a single-entry list from the older
    ``domain_accession``/``domain_start_aa``/``domain_end_aa`` summary
    fields (looking the accession's human-readable name up in
    ``run.gene_track["domains"]`` when available) so this also works
    against a ``results.json`` written before ``key_domains`` existed --
    the same span, just not yet carrying its own name/accession.
    """
    key_domains = run.summary.get("key_domains")
    if key_domains:
        return key_domains
    accession = run.summary.get("domain_accession")
    start = run.summary.get("domain_start_aa")
    end = run.summary.get("domain_end_aa")
    if start is None or end is None:
        return []
    name = accession
    for domain in (run.gene_track or {}).get("domains") or []:
        if domain.get("accession") == accession:
            name = domain.get("name") or accession
            break
    return [
        {
            "name": name or "target domain",
            "accession": accession,
            "start_aa": start,
            "end_aa": end,
        }
    ]


def _domain_highlight_color(index: int, domain_name: str) -> str:
    """Color for one configured key-domain highlight span.

    The first key domain always gets ``DOMAIN_HIGHLIGHT_COLOR`` (documented
    in the track's own legend as "configured key-domain span"); every
    additional key domain -- when a gene configures more than one, e.g. to
    also highlight an autoinhibitory or RAS-binding domain alongside the
    kinase domain -- gets its own stable color derived the same way
    :func:`cfh.reporting.fusion_schematic.partner_color` derives a
    partner's, so multiple domain highlights stay visually distinguishable
    without hardcoding a fixed-size color list. The legend entry covers
    this whole scheme (one shared meaning -- "a configured key-domain
    span" -- regardless of exactly which stable shade a given domain name
    hashes to), the same way the fusion schematic's legend documents
    partner colors as decorative without a swatch per partner.
    """
    if index == 0:
        return DOMAIN_HIGHLIGHT_COLOR
    return deterministic_color(domain_name, lightness=0.6, saturation=0.5)


def _domain_track_svg(run: RealBenchmarkRun, outlier_ids: set[str]) -> str:
    """Render breakpoints by quantitative domain-retention state.

    Shares its position-axis and exon-boundary convention
    (:func:`cfh.reporting.fusion_schematic.render_position_axis_svg`) with
    the fusion-transcript schematic so the two renderers can't drift apart
    on transcript-end or exon labels; both read the same gene-agnostic
    ``gene_track["exon_boundaries_aa"]`` field.
    """
    axis_left = 60.0
    axis_width = 800.0
    key_domains = _domain_track_key_domains(run)
    positions = [
        row["breakpoint_protein_position"] for row in run.rows if row["breakpoint_protein_position"]
    ]
    protein_length = (run.gene_track or {}).get("protein_length")
    domain_ends = [d["end_aa"] for d in key_domains if d.get("end_aa") is not None]
    maximum = max([*domain_ends, *positions, protein_length or 0, 1])
    scale = axis_width / maximum

    domain_band_height = 16.0
    domain_band_gap = 3.0
    top_margin = 34.0
    domains_block_height = len(key_domains) * (domain_band_height + domain_band_gap)
    backbone_y = top_margin + domains_block_height + 10.0

    domain_elements = []
    for index, domain in enumerate(key_domains):
        d_start = domain.get("start_aa")
        d_end = domain.get("end_aa")
        if d_start is None or d_end is None:
            continue
        band_top = top_margin + index * (domain_band_height + domain_band_gap)
        rect_x0 = axis_left + d_start * scale
        rect_width = max(2.0, (d_end - d_start) * scale)
        rect_x1 = rect_x0 + rect_width
        color = _domain_highlight_color(index, domain.get("name") or "domain")
        domain_elements.append(
            f'<rect x="{rect_x0:.1f}" y="{band_top:.1f}" '
            f'width="{rect_width:.1f}" height="{domain_band_height:.1f}" '
            f'fill="{color}" opacity="0.55"/>'
        )
        label = f"{domain.get('name') or 'domain'} ({d_start}-{d_end})"
        label_y = band_top + domain_band_height - 4
        near_right_edge = rect_x0 > axis_left + axis_width / 2
        if near_right_edge:
            domain_elements.append(
                f'<text x="{rect_x1:.1f}" y="{label_y:.1f}" font-family="sans-serif" '
                f'font-size="9.5" text-anchor="end">{label}</text>'
            )
        else:
            domain_elements.append(
                f'<text x="{rect_x0:.1f}" y="{label_y:.1f}" font-family="sans-serif" '
                f'font-size="9.5">{label}</text>'
            )

    dots_top = backbone_y + 8.0
    dots = []
    for index, row in enumerate(run.rows):
        position = row["breakpoint_protein_position"]
        if position is None:
            continue
        fraction = row.get("domain_retained_fraction")
        status = row.get("domain_status")
        is_outlier = row["event_id"] in outlier_ids
        # `domain_is_truncated` is derived as exactly `0.0 < fraction < 1.0`
        # (cfh.mapping.feature_mapper.calculate_domain_retention), so a
        # fraction strictly between 0 and 1 always implies truncated and
        # vice versa -- there is no real "partially retained but not
        # truncated" case to separately color.
        if row.get("domain_is_truncated") or status == "disrupted":
            color = TRUNCATED_COLOR
            status_word = "truncated"
        elif fraction == 0.0 or status == "lost":
            color = LOST_COLOR
            status_word = "lost"
        elif fraction == 1.0 or status == "retained":
            color = RETAINED_COLOR
            status_word = "retained"
        else:
            # domain_status == "unknown": no resolvable domain/breakpoint
            # coordinates for this event's target domain (see
            # classify_domain_retention), a genuine "insufficient data"
            # state distinct from a confirmed "lost" call. Never observed
            # in the committed BRAF/RET real-cohort runs (every row there
            # resolves to retained/lost/disrupted), so rather than invent
            # an extra legended color for an unreached case, skip drawing
            # a dot instead of implying a retention outcome that was never
            # actually determined.
            continue
        stroke = BREAKPOINT_COLOR if is_outlier else "none"
        stroke_width = "1.5" if is_outlier else "0"
        y = dots_top + (index % 5) * 7
        title_parts = [f"event {row['event_id']}"]
        if row.get("sample_id"):
            title_parts.append(f"sample {row['sample_id']}")
        if row.get("partner_gene"):
            title_parts.append(f"partner {row['partner_gene']}")
        title_parts.append(f"breakpoint aa {position}")
        title_parts.append(f"domain {status_word}")
        if is_outlier:
            title_parts.append("reference discrepancy")
        title = escape_xml_text("; ".join(title_parts))
        dots.append(
            f'<circle cx="{axis_left + position * scale:.1f}" cy="{y:.1f}" r="3" fill="{color}" '
            f'stroke="{stroke}" stroke-width="{stroke_width}"><title>{title}</title></circle>'
        )
    dots_bottom = dots_top + 4 * 7 + 3

    position_axis_y = dots_bottom + 12.0
    position_axis_elements = render_position_axis_svg(
        position_axis_y,
        int(maximum),
        scale,
        axis_left=axis_left,
        exon_boundaries=(run.gene_track or {}).get("exon_boundaries_aa") or [],
    )
    for tick in range(100, int(maximum), 100):
        x = axis_left + tick * scale
        position_axis_elements.append(
            f'<line x1="{x:.1f}" y1="{position_axis_y:.1f}" x2="{x:.1f}" '
            f'y2="{position_axis_y + 5:.1f}" stroke="{AXIS_COLOR}" stroke-width="1"/>'
        )
        position_axis_elements.append(
            f'<text x="{x:.1f}" y="{position_axis_y + 15:.1f}" font-family="sans-serif" '
            f'font-size="8" text-anchor="middle">{tick}</text>'
        )

    exon_boundaries = (run.gene_track or {}).get("exon_boundaries_aa") or []
    exon_tick_y = position_axis_y + 18.0
    exon_tick_label_height = 20.0 if exon_boundaries else 0.0

    legend_y = exon_tick_y + exon_tick_label_height + 15.0
    height = legend_y + 33.0

    return "\n".join(
        [
            f'<svg xmlns="http://www.w3.org/2000/svg" width="920" height="{height:.0f}" '
            f'viewBox="0 0 920 {height:.0f}">',
            f'<rect width="920" height="{height:.0f}" fill="white"/>',
            f'<text x="60" y="28" font-family="sans-serif" font-size="16">'
            f"{run.gene_symbol} domain-retention track</text>",
            *domain_elements,
            f'<line x1="{axis_left:.1f}" y1="{backbone_y:.1f}" x2="{axis_left + axis_width:.1f}" '
            f'y2="{backbone_y:.1f}" stroke="{AXIS_COLOR}" stroke-width="4"/>',
            *dots,
            *position_axis_elements,
            f'<circle cx="60" cy="{legend_y:.1f}" r="4" fill="{RETAINED_COLOR}"/>'
            f'<text x="70" y="{legend_y + 5:.1f}" font-family="sans-serif" font-size="12">'
            "fully retained</text>",
            f'<circle cx="180" cy="{legend_y:.1f}" r="4" fill="{TRUNCATED_COLOR}"/>'
            f'<text x="190" y="{legend_y + 5:.1f}" font-family="sans-serif" font-size="12">'
            "truncated</text>",
            f'<circle cx="275" cy="{legend_y:.1f}" r="4" fill="{LOST_COLOR}"/>'
            f'<text x="285" y="{legend_y + 5:.1f}" font-family="sans-serif" font-size="12">'
            "fully lost</text>",
            f'<circle cx="365" cy="{legend_y:.1f}" r="4" fill="white" stroke="{BREAKPOINT_COLOR}" '
            'stroke-width="1.5"/>'
            f'<text x="375" y="{legend_y + 5:.1f}" font-family="sans-serif" font-size="12">'
            "reference discrepancy</text>",
            f'<rect x="60" y="{legend_y + 13:.1f}" width="12" height="10" '
            f'fill="{DOMAIN_HIGHLIGHT_COLOR}" opacity="0.55"/>'
            f'<text x="80" y="{legend_y + 22:.1f}" font-family="sans-serif" font-size="12">'
            "configured key-domain span (name/range labeled above; a second "
            "or later configured domain gets its own distinct shade)</text>",
            "</svg>",
        ]
    )


def _comparison_svg(run: RealBenchmarkRun) -> str:
    reference = run.reference or {}
    metrics = [
        ("In-frame", reference.get("in_frame_percent", 0), run.summary["in_frame_percent"]),
        (
            "Domain retained",
            reference.get("domain_retained_percent", 0),
            run.summary["kinase_retained_percent"],
        ),
    ]
    elements = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="620" height="210" viewBox="0 0 620 210">',
        '<rect width="620" height="210" fill="white"/>',
        f'<text x="20" y="25" font-family="sans-serif" font-size="16">'
        f"Reference vs {run.study_id}</text>",
    ]
    for index, (label, ref_value, run_value) in enumerate(metrics):
        y = 55 + index * 70
        elements.extend(
            [
                f'<text x="20" y="{y}" font-family="sans-serif" font-size="12">{label}</text>',
                f'<rect x="140" y="{y - 14}" width="{ref_value * 3.8:.1f}" '
                f'height="16" fill="{REFERENCE_BAR_COLOR}"/>',
                f'<text x="530" y="{y}" font-family="sans-serif" font-size="12">'
                f"reference {ref_value:.1f}%</text>",
                f'<rect x="140" y="{y + 10}" width="{run_value * 3.8:.1f}" height="16" '
                f'fill="{RETAINED_COLOR}"/>',
                f'<text x="530" y="{y + 24}" font-family="sans-serif" font-size="12">'
                f"run {run_value:.1f}%</text>",
            ]
        )
    elements.extend(
        [
            f'<rect x="140" y="190" width="12" height="10" fill="{REFERENCE_BAR_COLOR}"/>',
            '<text x="157" y="199" font-family="sans-serif" font-size="11">reference</text>',
            f'<rect x="245" y="190" width="12" height="10" fill="{RETAINED_COLOR}"/>',
            '<text x="262" y="199" font-family="sans-serif" font-size="11">current run</text>',
            "</svg>",
        ]
    )
    return "\n".join(elements)


def _gene_pair_markdown_summary(run: RealBenchmarkRun) -> str:
    """Concise, checked-in-friendly report for a ``gene_pair`` (joint-partner)
    run -- the pair-enrichment counterpart to :func:`markdown_summary`, which
    assumes single-gene domain-retention fields this kind of run never has.
    """
    summary = run.summary
    gene5, gene3 = summary["gene_pair"]
    lines = [
        f"# {gene5}-{gene3} joint-partner fusion benchmark: {run.study_id}",
        "",
        f"Retrieved from public cBioPortal on {run.retrieved_at.date().isoformat()}.",
        "",
    ]
    mechanism_note = summary.get("mechanism_note")
    if mechanism_note:
        lines.extend(["## Mechanism", "", mechanism_note, ""])
    lines += [
        "## Method",
        "",
        "Structural-variant records were live-fetched for "
        f"{', '.join(summary['component_genes'])} (the pair member(s) with a curated "
        "single-gene config) through the same cBioPortal/Genome Nexus ingestion, "
        "normalization, and breakpoint-mapping pipeline used for single-gene analysis "
        "(see `run_real_benchmark`), then pooled (deduplicated by event id) and tested "
        f"via `JointPartnerMode` for whether the configured ordered pair {gene5}->{gene3} "
        "is enriched relative to a marginal-independence null.",
        "",
        "## Results",
        "",
        f"- Structural variants returned: {summary['raw_structural_variant_count']}",
        "- Eligible fusion events (determinable 5'/3' orientation): "
        f"{summary['eligible_event_count']}",
        f"- Observed {gene5}->{gene3} count: {summary['observed_count']}",
        f"- Expected under independence: {summary['expected_count']:.2f}",
        "- Fisher exact test (one-sided, greater): odds ratio "
        f"{_format_stat(summary['fisher_odds_ratio'])}, "
        f"p={_format_stat(summary['fisher_p_value'])}",
        f"- Enriched (p < 0.05): {summary['is_enriched']}",
        "",
    ]
    if run.warnings:
        lines.extend(["## Warnings", ""])
        lines.extend(f"- {warning}" for warning in run.warnings)
        lines.append("")
    return "\n".join(lines)


def _write_gene_pair_outputs(
    run: RealBenchmarkRun,
    output_dir: str | Path,
    *,
    cli_args: list[str] | None = None,
    run_id: str | None = None,
) -> dict[str, Path]:
    """Write a run directory for a ``gene_pair`` (joint-partner) benchmark.

    Distinct from :func:`write_outputs`: a gene-pair run has no domain data,
    breakpoint-retention classification, or PDF report -- it is the
    pair-enrichment result computed directly by ``JointPartnerMode``, so
    this writes only the artifacts meaningful for it (results.tsv/json,
    report.md, manifest.json), following the same
    ``runs/<type>_<ISO8601-timestamp>/`` naming convention as every other
    run type.
    """
    destination = Path(output_dir) / (run_id or _run_id(run))
    destination.mkdir(parents=True, exist_ok=True)
    tsv_path = destination / "results.tsv"
    json_path = destination / "results.json"
    markdown_path = destination / "report.md"
    manifest_path = destination / "manifest.json"

    _write_tsv(tsv_path, run.rows)

    algorithm_results = [result.model_dump(mode="json") for result in run.results]
    payload = _json_safe(
        {
            "gene_symbol": run.gene_symbol,
            "study_id": run.study_id,
            "molecular_profile_id": run.molecular_profile_id,
            "retrieved_at": run.retrieved_at.isoformat(),
            "summary": run.summary,
            "warnings": run.warnings,
            "events": run.rows,
            "algorithm_results": algorithm_results,
        }
    )
    json_path.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")
    markdown_path.write_text(_gene_pair_markdown_summary(run))
    manifest_path.write_text(
        json.dumps(
            {
                "gene": run.gene_symbol,
                "study_id": run.study_id,
                "endpoints_used": run.endpoints,
                "git_sha": _git_sha(),
                "cli_args": cli_args or [],
                "timestamp": run.retrieved_at.isoformat(),
            },
            indent=2,
        )
        + "\n"
    )
    return {
        "run_directory": destination,
        "manifest": manifest_path,
        "tsv": tsv_path,
        "json": json_path,
        "markdown": markdown_path,
    }


def write_outputs(
    run: RealBenchmarkRun,
    output_dir: str | Path,
    *,
    output_stem: str | None = None,
    cli_args: list[str] | None = None,
    run_id: str | None = None,
    pdf: bool = True,
) -> dict[str, Path]:
    """Write a complete, provenance-bearing run directory."""
    del output_stem  # retained as a compatibility-only keyword for older callers
    if run.is_gene_pair:
        return _write_gene_pair_outputs(run, output_dir, cli_args=cli_args, run_id=run_id)
    destination = Path(output_dir) / (run_id or _run_id(run))
    visualization_dir = destination / "visualizations"
    visualization_dir.mkdir(parents=True, exist_ok=True)
    tsv_path = destination / "results.tsv"
    json_path = destination / "results.json"
    markdown_path = destination / "report.md"

    _write_tsv(tsv_path, run.rows)

    algorithm_results = [result.model_dump(mode="json") for result in run.results]
    for result in algorithm_results:
        null_rates = (result.get("Tables") or {}).get("permutation_null_retention_rates")
        if null_rates is not None:
            result["Tables"]["permutation_null_retention_rates"] = {
                "omitted_from_artifact": True,
                "count": len(null_rates),
            }
    payload = _json_safe(
        {
            "gene_symbol": run.gene_symbol,
            "study_id": run.study_id,
            "molecular_profile_id": run.molecular_profile_id,
            "retrieved_at": run.retrieved_at.isoformat(),
            "summary": run.summary,
            "warnings": run.warnings,
            "events": run.rows,
            "algorithm_results": algorithm_results,
            "reference": run.reference,
            "gene_track": run.gene_track,
            "intragenic_deletions": run.intragenic_deletions,
        }
    )
    json_path.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")

    fusion_schematic_svg = None
    schematic_markup = render_fusion_schematic_svg(payload)
    if schematic_markup is not None:
        fusion_schematic_svg = visualization_dir / "fusion_schematic.svg"
        fusion_schematic_svg.write_text(schematic_markup + "\n")

    intragenic_deletion_svg = None
    deletion_markup = render_intragenic_deletion_schematic_svg(payload)
    if deletion_markup is not None:
        intragenic_deletion_svg = visualization_dir / "intragenic_deletion_schematic.svg"
        intragenic_deletion_svg.write_text(deletion_markup + "\n")

    discrepancies = _discrepancies(run)
    outliers_path = destination / "outliers.tsv"
    _write_tsv(
        outliers_path,
        discrepancies,
        [
            "discrepancy_type",
            "event_id",
            "partner_gene",
            "frame_status",
            "breakpoint_protein_position",
            "breakpoint_exon",
            "is_intronic_breakpoint",
            "retained_domains",
            "lost_domains",
            "disrupted_domains",
            "source_annotation_text",
            "source_site2_effect_on_frame",
        ],
    )
    reference_ids = {
        row["event_id"]
        for row in discrepancies
        if row["discrepancy_type"] == "reference_discrepancy"
    }
    domain_svg = visualization_dir / "domain_retention_outliers.svg"
    comparison_svg = visualization_dir / "reference_comparison.svg"
    markdown_path.write_text(
        markdown_summary(
            run,
            domain_svg_path=domain_svg.relative_to(destination).as_posix(),
            comparison_svg_path=comparison_svg.relative_to(destination).as_posix(),
            fusion_schematic_svg_path=(
                fusion_schematic_svg.relative_to(destination).as_posix()
                if fusion_schematic_svg
                else None
            ),
            intragenic_deletion_svg_path=(
                intragenic_deletion_svg.relative_to(destination).as_posix()
                if intragenic_deletion_svg
                else None
            ),
        )
    )
    domain_svg.write_text(_domain_track_svg(run, reference_ids) + "\n")
    comparison_svg.write_text(_comparison_svg(run) + "\n")
    manifest_path = destination / "manifest.json"
    manifest_path.write_text(
        json.dumps(
            {
                "gene": run.gene_symbol,
                "study_id": run.study_id,
                "endpoints_used": run.endpoints,
                "git_sha": _git_sha(),
                "cli_args": cli_args or [],
                "timestamp": run.retrieved_at.isoformat(),
            },
            indent=2,
        )
        + "\n"
    )
    paths = {
        "run_directory": destination,
        "manifest": manifest_path,
        "tsv": tsv_path,
        "json": json_path,
        "markdown": markdown_path,
        "outliers": outliers_path,
        "domain_svg": domain_svg,
        "comparison_svg": comparison_svg,
    }
    if fusion_schematic_svg:
        paths["fusion_schematic_svg"] = fusion_schematic_svg
    if intragenic_deletion_svg:
        paths["intragenic_deletion_svg"] = intragenic_deletion_svg
    if pdf:
        pdf_path = destination / "report.pdf"
        try:
            render_pdf_report(
                payload,
                pdf_path,
                results_tsv_path=tsv_path,
                visualizations_dir=visualization_dir,
            )
        except Exception as exc:
            # report.pdf is a convenience rendering of data already fully
            # captured in results.json/results.tsv/report.md -- a failure
            # here (e.g. a reportlab table-layout edge case for a
            # many-domain gene's wide per-event table) must not take down
            # an otherwise-successful benchmark run, especially inside a
            # genome-wide cohort scan where one gene's PDF failing would
            # otherwise abort every other gene's already-completed results.
            run.warnings.append(
                f"report.pdf could not be rendered ({type(exc).__name__}: {exc}); "
                "all other outputs (results.json/results.tsv/report.md) were written "
                "successfully and are unaffected."
            )
            if pdf_path.exists():
                pdf_path.unlink()
        else:
            paths["pdf"] = pdf_path
    return paths
