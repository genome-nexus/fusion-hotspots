"""Per-gene biological configuration registry.

Gene-specific facts (canonical transcript, protein accession, domain
boundaries, ...) live in YAML files under ``genes/configs/`` and are loaded
into a validated :class:`GeneConfig`. Generic pipeline code must never
hardcode a gene's biology directly; it should always receive a
``GeneConfig`` instance instead.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
from typing import TYPE_CHECKING, Optional

import yaml
from pydantic import BaseModel, ConfigDict, model_validator

if TYPE_CHECKING:
    from cfh.mapping.genome_nexus_source import CanonicalTranscript, PfamDomain

CONFIGS_DIR = Path(__file__).parent / "configs"


class KeyDomain(BaseModel):
    model_config = ConfigDict(extra="forbid")

    name: str
    source: str
    key: Optional[str] = None
    accession: Optional[str] = None
    """Source-native identifier, such as a Pfam accession from Genome Nexus."""
    start_aa: Optional[int] = None
    end_aa: Optional[int] = None
    """Opt-in curated coordinate override (one-based, inclusive). When both
    are set, ``cfh.mapping.feature_mapper.map_event`` uses these coordinates
    directly instead of looking the domain up by name/accession against a
    live domain source (UniProt) -- for a region with no live-matchable
    feature at all (e.g. a regulatory/autoinhibitory element that isn't
    annotated as a discrete "domain"/"region" by UniProt, only described in
    the primary literature by residue range). Leave unset for the normal,
    live-matched path; only set both together, never just one."""

    @model_validator(mode="after")
    def _validate_coordinate_pair(self) -> "KeyDomain":
        if (self.start_aa is None) != (self.end_aa is None):
            raise ValueError(
                "KeyDomain.start_aa and end_aa must be set together, or both left unset"
            )
        if self.start_aa is not None and self.end_aa is not None and self.end_aa < self.start_aa:
            raise ValueError("KeyDomain.end_aa must be >= start_aa")
        return self


class MutualExclusivityTarget(BaseModel):
    """One opt-in comparator alteration to test a gene's fusion status
    against for cohort-wide co-occurrence/mutual-exclusivity (Fisher's
    exact), via the ``mutation_cooccurrence`` algorithm."""

    model_config = ConfigDict(extra="forbid")

    gene: str
    alteration_type: str
    """E.g. ``point_mutation``, ``cna_amp``, ``cna_del`` -- see
    ``cfh.normalization.alteration_normalizer.CNA_ALTERATION_TYPE_BY_CODE``
    for the full copy-number vocabulary."""
    entrez_gene_id: Optional[int] = None
    """Required only for a caller (e.g. ``run_real_benchmark``) that wants
    to live-fetch this target's alteration data automatically; a target
    left without it is still valid config but cannot be auto-fetched."""
    protein_change: Optional[str] = None
    """Optional point-mutation filter, e.g. ``"V600E"``. Only applicable
    when ``alteration_type == "point_mutation"``; ignored otherwise."""


class BenchmarkReference(BaseModel):
    """Optional literature baseline used by reports and discrepancy artifacts."""

    model_config = ConfigDict(extra="forbid")

    citation: str
    fusion_count: int
    in_frame_percent: float
    domain_retained_percent: float


class GeneConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")

    gene_symbol: Optional[str] = None
    gene_pair: Optional[tuple[str, str]] = None
    canonical_transcript_id: Optional[str] = None
    protein_id: Optional[str] = None
    key_domains: list[KeyDomain] = []
    disruption_required_domains: list[KeyDomain] = []
    """Domains a fusion must LOSE/DISRUPT (not retain) to be functionally
    oncogenic, e.g. an autoinhibitory N-terminal module. Opt-in per gene,
    parallel to ``key_domains``; the domain-disruption algorithm no-ops
    when this is left empty."""
    expected_retained_exon_hint: Optional[str] = None
    analysis_modes: list[str] = []
    entrez_gene_id: Optional[int] = None
    benchmark_reference: Optional[BenchmarkReference] = None
    """NCBI Entrez gene id, e.g. for cBioPortal structural-variant API queries."""
    mutual_exclusivity_targets: list[MutualExclusivityTarget] = []
    """Opt-in comparator alteration(s) -- e.g. this gene's own hotspot point
    mutation, or another pathway member's alteration -- to test this gene's
    fusion-positive sample set against for cohort-wide co-occurrence/mutual-
    exclusivity (Fisher's exact), via the ``mutation_cooccurrence``
    algorithm. Empty by default; the algorithm gracefully no-ops for genes
    that don't configure this, the same opt-in pattern as
    ``disruption_required_domains``/``gene_pair``."""
    mechanism_note: Optional[str] = None
    """Opt-in free-text description of *why* this fusion is oncogenic, surfaced
    verbatim in a ``gene_pair`` run's report (see ``_gene_pair_markdown_summary``)
    so the report never has to guess at or default to a mechanism (e.g.
    domain-retention) that this gene/pair's algorithm results don't actually
    support. Left unset by default -- the same opt-in, no-op-when-absent
    pattern as every other field above."""
    mechanism_note_verified: bool = True
    """Whether ``mechanism_note`` has been cross-checked by a human against
    primary literature (the standard every currently-curated gene's note
    meets), as opposed to drafted from an AI literature-synthesis tool (e.g.
    OpenEvidence) and not yet independently verified. Defaults to ``True`` so
    every existing curated config keeps its current "Curated mechanism"
    report/viewer treatment without having to set this explicitly; a config
    that pastes in an unverified AI-drafted note must set this to ``False``
    so the report/viewer render it as "AI-suggested mechanism (unverified)"
    instead, never conflating the two trust levels. Meaningless when
    ``mechanism_note`` is unset."""

    @model_validator(mode="after")
    def _validate_config_target(self) -> "GeneConfig":
        """Require either a single gene or an ordered fusion-gene pair."""
        if self.gene_symbol is None and self.gene_pair is None:
            raise ValueError("GeneConfig requires gene_symbol or gene_pair")
        if self.gene_pair is None and (
            self.canonical_transcript_id is None or self.protein_id is None
        ):
            raise ValueError(
                "single-gene GeneConfig requires canonical_transcript_id and protein_id"
            )
        return self


def _config_path(gene_symbol: str) -> Path:
    return CONFIGS_DIR / f"{gene_symbol.lower()}.yaml"


def load_gene_config(gene_symbol: str) -> GeneConfig:
    """Load and validate a gene's YAML config by symbol (e.g. ``"braf"``)."""
    path = _config_path(gene_symbol)
    if not path.exists():
        raise FileNotFoundError(f"No gene config found for {gene_symbol!r} at {path}")
    with path.open() as fh:
        raw = yaml.safe_load(fh)
    return GeneConfig.model_validate(raw)


def _pfam_domain_for_key_domain(
    key_domain: KeyDomain, pfam_domains: list[PfamDomain]
) -> PfamDomain | None:
    """Match a configured domain to its coordinate-bearing Pfam record."""
    if key_domain.accession:
        return next(
            (domain for domain in pfam_domains if domain.pfam_id == key_domain.accession), None
        )
    return next((domain for domain in pfam_domains if domain.pfam_id == key_domain.name), None)


def derive_gene_config_defaults(
    config: GeneConfig,
    canonical: CanonicalTranscript,
    *,
    domain_name_resolver: Callable[[str], str | None] | None = None,
) -> GeneConfig:
    """Fill coordinate-derived analysis defaults without replacing curated values.

    The most N-terminal configured key domain is the retention boundary when a
    gene has several key domains.  This is the loosest (earliest) constraint:
    retaining its first exon preserves the N-terminus of every more C-terminal
    key domain as well.

    An explicit ``disruption_required_domains`` YAML value, including an empty
    list, is authoritative.  ``expected_retained_exon_hint: null`` is treated as
    unset, while any non-null human-supplied hint wins.
    """
    from cfh.mapping.genome_nexus_source import cds_bounds_from_utrs, exon_protein_boundaries

    matched_key_domains = [
        match
        for key_domain in config.key_domains
        if (match := _pfam_domain_for_key_domain(key_domain, canonical.pfam_domains)) is not None
    ]
    if not matched_key_domains:
        return config
    retention_domain = min(
        matched_key_domains,
        key=lambda domain: (domain.start_aa, domain.end_aa, domain.pfam_id),
    )

    updates: dict[str, object] = {}
    if config.expected_retained_exon_hint is None:
        cds_min, cds_max = cds_bounds_from_utrs(canonical.utrs)
        boundaries = exon_protein_boundaries(
            canonical.exons,
            cds_min_genomic=cds_min,
            cds_max_genomic=cds_max,
        )
        containing = [
            boundary
            for boundary in boundaries
            if boundary.start_aa <= retention_domain.start_aa <= boundary.end_aa
        ]
        preceding = [
            boundary for boundary in boundaries if boundary.end_aa <= retention_domain.start_aa
        ]
        chosen_boundary = (
            min(containing, key=lambda boundary: boundary.exon_rank)
            if containing
            else max(
                preceding,
                key=lambda boundary: (boundary.end_aa, boundary.exon_rank),
                default=None,
            )
        )
        if chosen_boundary is not None:
            updates["expected_retained_exon_hint"] = str(chosen_boundary.exon_rank)

    if "disruption_required_domains" not in config.model_fields_set:
        candidates = sorted(
            (
                domain
                for domain in canonical.pfam_domains
                if domain.end_aa <= retention_domain.start_aa
            ),
            key=lambda domain: (domain.start_aa, domain.end_aa, domain.pfam_id),
        )
        # GeneConfig domains are addressed by Pfam accession, so repeated
        # occurrences of one family (for example tandem domains) are one
        # configurable disruption target rather than duplicate entries.
        candidates_by_accession = {domain.pfam_id: domain for domain in candidates}
        updates["disruption_required_domains"] = [
            KeyDomain(
                name=(domain_name_resolver(domain.pfam_id) if domain_name_resolver else None)
                or domain.pfam_id,
                source="genome_nexus",
                key=f"auto_disruption_{domain.pfam_id.lower()}",
                accession=domain.pfam_id,
            )
            for domain in candidates_by_accession.values()
        ]

    return config.model_copy(update=updates) if updates else config


def available_genes() -> list[str]:
    """List gene symbols with a registered config."""
    return sorted(p.stem.upper() for p in CONFIGS_DIR.glob("*.yaml"))
