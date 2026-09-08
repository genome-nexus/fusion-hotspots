"""Resolve a fusion breakpoint to a transcript/exon position.

Primary source: exon numbers already present in the cBioPortal SV
annotation text. Fallbacks, in order: Genome Nexus's canonical-transcript
exon structure (real breakpoint-to-protein-position mapping from a genomic
coordinate, keyed only by gene symbol), then Ensembl's REST overlap API
for the protein's features (keyed by an Ensembl protein id), used when the
annotation doesn't mention an exon at all.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Any

import requests

from cfh.genes.registry import GeneConfig
from cfh.mapping.genome_nexus_source import (
    GenomeNexusClient,
    GenomeNexusGeneNotFound,
    cds_bounds_from_utrs,
    map_genomic_breakpoint_to_protein_position,
    parse_canonical_transcript,
)

_EXON_PATTERN = re.compile(r"exon\s*(\d+)", re.IGNORECASE)


class EnsemblFallbackUnavailable(RuntimeError):
    """Raised when the Ensembl-fallback mapping cannot even be attempted.

    This is the case when the SV annotation has no exon number AND the
    caller has no ``ensembl_protein_id`` to query Ensembl with. Callers
    should treat this as "fallback unavailable" and must not confuse it
    with "queried Ensembl and it had nothing to say" (which is a real,
    non-exceptional result -- see :func:`resolve_transcript_mapping`).
    """


@dataclass
class TranscriptMapping:
    transcript_id: str
    breakpoint_exon: int | None
    source: str  # "annotation" | "genome_nexus_fallback" | "ensembl_fallback"
    ensembl_features: list[dict] | None = None
    breakpoint_protein_features: list[dict] | None = None
    """Ensembl protein features (e.g. domains) whose range contains the
    supplied breakpoint amino-acid coordinate, when one was given. This is
    the real, computed result of the Ensembl fallback: the
    ``overlap/translation`` endpoint exposes protein-feature boundaries,
    not an exon-to-protein-coordinate map, so ``breakpoint_exon`` genuinely
    stays ``None`` on this path -- that is a data-source limitation, not a
    silently-dropped answer.
    """
    breakpoint_protein_position: int | None = None
    """Estimated protein residue for the breakpoint, computed from Genome
    Nexus exon coordinates via :func:`~cfh.mapping.genome_nexus_source.
    map_genomic_breakpoint_to_protein_position` (source ==
    ``"genome_nexus_fallback"``). Unlike the Ensembl overlap fallback, this
    path produces a real position estimate, not just an "is there a
    feature here" check.
    """
    is_intronic_breakpoint: bool | None = None


class EnsemblClient:
    BASE_URL = "https://rest.ensembl.org"

    def __init__(self, session: "requests.Session | None" = None):
        self.session = session or requests.Session()

    def fetch_protein_features(self, ensembl_protein_id: str) -> list[dict]:
        url = f"{self.BASE_URL}/overlap/translation/{ensembl_protein_id}"
        response = self.session.get(
            url,
            params={"feature": "protein_feature"},
            headers={"Content-Type": "application/json"},
            timeout=30,
        )
        response.raise_for_status()
        return response.json()


def extract_exon_from_annotation(annotation: str | None) -> int | None:
    """Pull an explicit exon number out of free-text SV annotation, if present."""
    if not annotation:
        return None
    match = _EXON_PATTERN.search(annotation)
    return int(match.group(1)) if match else None


def _features_overlapping(features: list[dict], breakpoint_aa: int) -> list[dict]:
    overlapping = []
    for feature in features:
        start, end = feature.get("start"), feature.get("end")
        if start is None or end is None:
            continue
        if start <= breakpoint_aa <= end:
            overlapping.append(feature)
    return overlapping


def resolve_transcript_mapping(
    annotation: str | None,
    gene_config: GeneConfig,
    *,
    breakpoint_aa: int | None = None,
    ensembl_protein_id: str | None = None,
    ensembl_client: "EnsemblClient | Any | None" = None,
) -> TranscriptMapping:
    """Resolve breakpoint exon info, preferring the SV annotation text.

    When the annotation lacks an explicit exon number, this falls back to
    fetching Ensembl protein features for ``ensembl_protein_id`` and, when a
    ``breakpoint_aa`` protein-coordinate is supplied, actually uses those
    fetched features to report which ones overlap the breakpoint
    (``breakpoint_protein_features``) -- a real result computed from the
    fallback data source, not a silent ``None``.

    Raises:
        EnsemblFallbackUnavailable: the annotation has no exon number and
            no ``ensembl_protein_id`` was supplied, so the fallback cannot
            be attempted at all.
    """
    exon = extract_exon_from_annotation(annotation)
    if exon is not None:
        return TranscriptMapping(
            transcript_id=gene_config.canonical_transcript_id,
            breakpoint_exon=exon,
            source="annotation",
        )

    if not ensembl_protein_id:
        raise EnsemblFallbackUnavailable(
            "annotation has no exon number and no ensembl_protein_id was "
            "supplied, so the Ensembl fallback cannot be attempted"
        )

    client = ensembl_client or EnsemblClient()
    features = client.fetch_protein_features(ensembl_protein_id)

    breakpoint_protein_features = (
        _features_overlapping(features, breakpoint_aa) if breakpoint_aa is not None else None
    )

    return TranscriptMapping(
        transcript_id=gene_config.canonical_transcript_id,
        breakpoint_exon=None,
        source="ensembl_fallback",
        ensembl_features=features,
        breakpoint_protein_features=breakpoint_protein_features,
    )


def resolve_breakpoint_protein_position(
    annotation: str | None,
    gene_config: GeneConfig,
    *,
    breakpoint_genomic: int | None = None,
    breakpoint_aa: int | None = None,
    genome_nexus_client: "GenomeNexusClient | None" = None,
    ensembl_protein_id: str | None = None,
    ensembl_client: "EnsemblClient | Any | None" = None,
    role: str | None = None,
) -> TranscriptMapping:
    """Resolve a breakpoint to a real protein-position estimate where possible.

    Fallback order when the annotation has no exon number:

    1. Genome Nexus's canonical-transcript endpoint (keyed only by
       ``gene_config.gene_symbol`` -- no pre-existing Ensembl protein id
       needed). When a ``breakpoint_genomic`` coordinate is supplied, this
       uses the transcript's real exon coordinates -- clipped to their
       coding-only portion via the transcript's own UTR annotations
       (:func:`~cfh.mapping.genome_nexus_source.cds_bounds_from_utrs`) --
       to compute an actual protein-position estimate
       (:func:`~cfh.mapping.genome_nexus_source.map_genomic_breakpoint_to_protein_position`),
       not just a "some feature overlaps here" flag.
    2. If Genome Nexus has no canonical-transcript mapping for this gene at
       all (:class:`~cfh.mapping.genome_nexus_source.GenomeNexusGeneNotFound`),
       or found the gene but had no exon data to compute a position from,
       falls back to the Ensembl ``overlap/translation`` protein-feature
       check from :func:`resolve_transcript_mapping` -- but only when a
       protein-coordinate ``breakpoint_aa`` is supplied, since that
       endpoint has no way to consume a genomic coordinate (it exposes
       protein-feature boundaries, not an exon/CDS map). Silently
       returning an apparently-successful result with no computed position
       when only ``breakpoint_genomic`` is available would misrepresent
       "we don't have the data for this" as "we checked and found
       nothing". The Ensembl protein id for this fallback is
       ``ensembl_protein_id`` if the caller supplied one, otherwise the
       ``protein_id`` from the Genome Nexus canonical-transcript payload
       (when one was found) -- no need for a caller to redundantly supply
       an id that Genome Nexus already gave us.
    3. If neither can produce or attempt a mapping, raises
       :class:`EnsemblFallbackUnavailable` -- the genuinely unmappable
       case.

    ``role`` is the target gene's role in the fusion ("five_prime" or
    "three_prime"), when known -- forwarded to
    :func:`~cfh.mapping.genome_nexus_source.map_genomic_breakpoint_to_protein_position`
    to disambiguate which side of an intronic breakpoint to snap to.
    Leaving it ``None`` (the default) keeps the direction-unaware
    nearest-exon fallback for callers with no fusion-partner context.
    """
    exon = extract_exon_from_annotation(annotation)
    if exon is not None:
        return TranscriptMapping(
            transcript_id=gene_config.canonical_transcript_id,
            breakpoint_exon=exon,
            source="annotation",
        )

    client = genome_nexus_client or GenomeNexusClient()
    try:
        payload = client.fetch_canonical_transcript(gene_config.gene_symbol)
    except GenomeNexusGeneNotFound:
        payload = None

    canonical = None
    if payload is not None:
        canonical = parse_canonical_transcript(payload)
        if breakpoint_genomic is not None and canonical.exons:
            strand = canonical.exons[0].strand
            cds_min_genomic, cds_max_genomic = cds_bounds_from_utrs(canonical.utrs)
            position = map_genomic_breakpoint_to_protein_position(
                canonical.exons,
                breakpoint_genomic,
                strand,
                cds_min_genomic=cds_min_genomic,
                cds_max_genomic=cds_max_genomic,
                role=role,
            )
            return TranscriptMapping(
                transcript_id=canonical.refseq_mrna_id or canonical.transcript_id,
                breakpoint_exon=position.exon_rank,
                source="genome_nexus_fallback",
                breakpoint_protein_position=position.protein_position,
                is_intronic_breakpoint=position.is_intronic,
            )
        # Genome Nexus knows the gene, but we either have no genomic
        # breakpoint to map or it returned no exon data -- fall through to
        # the Ensembl protein-feature check (if a protein coordinate is
        # available) rather than returning a hollow "success" here.

    effective_ensembl_protein_id = ensembl_protein_id or (
        canonical.protein_id if canonical is not None else None
    )

    if effective_ensembl_protein_id and breakpoint_aa is not None:
        return resolve_transcript_mapping(
            annotation,
            gene_config,
            breakpoint_aa=breakpoint_aa,
            ensembl_protein_id=effective_ensembl_protein_id,
            ensembl_client=ensembl_client,
        )

    raise EnsemblFallbackUnavailable(
        "could not compute a real breakpoint position: Genome Nexus either has no "
        "canonical-transcript mapping for this gene or returned no usable exon data, "
        "and no protein-coordinate breakpoint_aa with a usable Ensembl protein id "
        "(explicit ensembl_protein_id, or protein_id from a found Genome Nexus "
        "canonical transcript) was available for the Ensembl protein-feature "
        "fallback -- a genomic breakpoint coordinate alone cannot be mapped through "
        "that endpoint"
    )
