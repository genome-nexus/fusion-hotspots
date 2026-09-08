"""Genome Nexus integration: canonical-transcript/Pfam-domain lookup by gene
symbol, and genomic-breakpoint-to-protein-position mapping using a
transcript's exon structure.

Generic by construction: every function here takes a gene symbol,
accession, or exon list supplied by the caller (typically from a
``GeneConfig`` or a fetched Genome Nexus payload) -- nothing gene-specific
is hardcoded.

Genome Nexus's own variant-annotation endpoints are for point mutations
expressed as HGVS/genomic-change strings, not structural-variant fusion
breakpoints -- there is no endpoint that maps a fusion breakpoint to a
protein position directly. This module does that arithmetic itself (see
:func:`map_genomic_breakpoint_to_protein_position`) from the exon
coordinates Genome Nexus does provide.
"""

from __future__ import annotations

import math
import time
from dataclasses import dataclass

import requests

from cfh.mapping.domain_source import ProteinDomain

DEFAULT_BASE_URL = "https://www.genomenexus.org"
DEFAULT_ISOFORM_OVERRIDE_SOURCE = "mskcc"

_RETRYABLE_STATUS_CODES = {429, 500, 502, 503, 504}


class GenomeNexusGeneNotFound(RuntimeError):
    """Raised when Genome Nexus has no canonical-transcript mapping at all
    for a requested gene symbol (a 404 from the canonical-transcript
    endpoint). Callers should treat this as "try the next domain/transcript
    source" (e.g. UniProt-direct), not as a transient failure.
    """


@dataclass
class PfamDomain:
    pfam_id: str
    start_aa: int
    end_aa: int


@dataclass
class ExonRecord:
    exon_id: str
    start: int
    end: int
    rank: int
    strand: int


@dataclass
class UtrRecord:
    utr_type: str  # "five_prime_UTR" | "three_prime_UTR"
    start: int
    end: int
    strand: int


@dataclass
class CanonicalTranscript:
    transcript_id: str
    refseq_mrna_id: str | None
    protein_id: str | None
    protein_length: int | None
    uniprot_id: str | None
    pfam_domains: list[PfamDomain]
    exons: list[ExonRecord]
    utrs: list[UtrRecord]


@dataclass
class BreakpointProteinPosition:
    protein_position: int
    exon_rank: int
    is_intronic: bool
    """True if the breakpoint fell in an intron and had to be clamped to
    the nearest exon boundary rather than landing inside a coding exon.
    """


class GenomeNexusClient:
    """Thin REST client with basic retry/backoff and a per-process cache.

    The public Genome Nexus instance documents no rate limit and requires
    no auth, but repeatedly querying the same gene/transcript in a batch
    pipeline is still wasteful and impolite, so results are cached by
    request key for the life of the client instance.
    """

    def __init__(
        self,
        base_url: str = DEFAULT_BASE_URL,
        session: "requests.Session | None" = None,
        max_retries: int = 3,
        backoff_seconds: float = 0.5,
    ):
        self.base_url = base_url.rstrip("/")
        self.session = session or requests.Session()
        self.max_retries = max_retries
        self.backoff_seconds = backoff_seconds
        self._canonical_transcript_cache: dict[str, dict] = {}
        self._transcript_cache: dict[str, dict] = {}

    def _get(self, path: str, *, params: dict | None = None) -> requests.Response:
        url = f"{self.base_url}{path}"
        attempt = 0
        while True:
            response = self.session.get(url, params=params, timeout=30)
            if response.status_code not in _RETRYABLE_STATUS_CODES or attempt >= self.max_retries:
                return response
            time.sleep(self.backoff_seconds * (2**attempt))
            attempt += 1

    def fetch_canonical_transcript(
        self, hugo_symbol: str, *, isoform_override_source: str = DEFAULT_ISOFORM_OVERRIDE_SOURCE
    ) -> dict:
        """GET ``/ensembl/canonical-transcript/hgnc/{hugoSymbol}``, cached by symbol."""
        cache_key = f"{hugo_symbol}:{isoform_override_source}"
        if cache_key in self._canonical_transcript_cache:
            return self._canonical_transcript_cache[cache_key]

        response = self._get(
            f"/ensembl/canonical-transcript/hgnc/{hugo_symbol}",
            params={"isoformOverrideSource": isoform_override_source},
        )
        if response.status_code == 404:
            raise GenomeNexusGeneNotFound(
                f"Genome Nexus has no canonical-transcript mapping for gene symbol {hugo_symbol!r}"
            )
        response.raise_for_status()
        payload = response.json()
        self._canonical_transcript_cache[cache_key] = payload
        return payload

    def fetch_transcript(self, transcript_id: str) -> dict:
        """GET ``/ensembl/transcript/{transcriptId}``, cached by transcript id."""
        if transcript_id in self._transcript_cache:
            return self._transcript_cache[transcript_id]

        response = self._get(f"/ensembl/transcript/{transcript_id}")
        response.raise_for_status()
        payload = response.json()
        self._transcript_cache[transcript_id] = payload
        return payload


def parse_canonical_transcript(payload: dict) -> CanonicalTranscript:
    """Parse a canonical-transcript (or transcript-by-id) Genome Nexus payload."""
    # Some real Genome Nexus canonical-transcript payloads (observed live on
    # a genome-wide gene set, e.g. PALB2/BBC3/INPP4B) include an empty {}
    # entry in pfamDomains alongside otherwise well-formed ones; skip any
    # entry missing a required field rather than crash the whole parse.
    pfam_domains = [
        PfamDomain(
            pfam_id=domain["pfamDomainId"],
            start_aa=domain["pfamDomainStart"],
            end_aa=domain["pfamDomainEnd"],
        )
        for domain in payload.get("pfamDomains", [])
        if domain.get("pfamDomainId") is not None
        and domain.get("pfamDomainStart") is not None
        and domain.get("pfamDomainEnd") is not None
    ]
    exons = [
        ExonRecord(
            exon_id=exon["exonId"],
            start=exon["exonStart"],
            end=exon["exonEnd"],
            rank=exon["rank"],
            strand=exon["strand"],
        )
        for exon in payload.get("exons", [])
    ]
    utrs = [
        UtrRecord(
            utr_type=utr["type"],
            start=utr["start"],
            end=utr["end"],
            strand=utr["strand"],
        )
        for utr in payload.get("utrs", [])
    ]
    return CanonicalTranscript(
        transcript_id=payload["transcriptId"],
        refseq_mrna_id=payload.get("refseqMrnaId"),
        protein_id=payload.get("proteinId"),
        protein_length=payload.get("proteinLength"),
        uniprot_id=payload.get("uniprotId"),
        pfam_domains=pfam_domains,
        exons=exons,
        utrs=utrs,
    )


def pfam_domains_to_protein_domains(canonical: CanonicalTranscript) -> list[ProteinDomain]:
    """Adapt Genome Nexus Pfam domains to the shared :class:`ProteinDomain` shape."""
    return [
        ProteinDomain(
            name=domain.pfam_id,
            start_aa=domain.start_aa,
            end_aa=domain.end_aa,
            source="genome_nexus",
            accession=domain.pfam_id,
        )
        for domain in canonical.pfam_domains
    ]


def resolve_domains(
    gene_symbol: str,
    uniprot_accession: str | None,
    *,
    genome_nexus_client: "GenomeNexusClient | None" = None,
    uniprot_source=None,
) -> list[ProteinDomain]:
    """Domain-annotation source chain: Genome Nexus first, UniProt-direct fallback.

    Genome Nexus's canonical-transcript endpoint is queried by gene symbol
    and needs no accession lookup step, so it's tried first. If it has no
    mapping for this gene at all (:class:`GenomeNexusGeneNotFound`), this
    falls back to a direct UniProt lookup by accession rather than raising,
    so a gene missing from Genome Nexus still gets domain annotations.
    """
    client = genome_nexus_client or GenomeNexusClient()
    try:
        payload = client.fetch_canonical_transcript(gene_symbol)
    except GenomeNexusGeneNotFound:
        payload = None

    if payload is not None:
        canonical = parse_canonical_transcript(payload)
        domains = pfam_domains_to_protein_domains(canonical)
        if domains:
            return domains

    if uniprot_accession is None:
        return []

    from cfh.mapping.domain_source import UniProtDomainSource

    source = uniprot_source or UniProtDomainSource()
    return source.fetch(uniprot_accession)


def cds_bounds_from_utrs(utrs: list[UtrRecord]) -> tuple[int | None, int | None]:
    """Derive ``(cds_min_genomic, cds_max_genomic)`` genomic bounds from a
    transcript's annotated UTRs, so an exon that partly overlaps a UTR can
    be clipped to its coding-only portion.

    A UTR of a given type can be split across multiple exons and so appear
    as several, non-adjacent segments in the payload (e.g. a plus-strand
    gene whose first two exons are both entirely or partly 5' UTR). The
    correct CDS boundary is adjacent to the OUTERMOST such segment -- the
    one farthest from the CDS -- not merely the first one encountered in
    the payload's (unordered) list, so every segment of each type must be
    aggregated, not just one.

    Strand-generic: which UTR type sits at the genomic-low vs. genomic-high
    end depends on strand, not on which gene this is.

    * Plus strand: transcription runs low-to-high genomic coordinates, so
      all 5' UTR segments sit below the CDS and all 3' UTR segments sit
      above it. The CDS start is one past the highest (max) end among all
      5' UTR segments; the CDS end is one before the lowest (min) start
      among all 3' UTR segments.
    * Minus strand: transcription runs high-to-low, so the roles invert --
      5' UTR segments sit above the CDS, 3' UTR segments below it.

    Either bound is ``None`` (no clipping on that end) when no UTR segment
    of the corresponding type is present in the payload.
    """
    five_prime_utrs = [u for u in utrs if u.utr_type == "five_prime_UTR"]
    three_prime_utrs = [u for u in utrs if u.utr_type == "three_prime_UTR"]
    any_utr = five_prime_utrs or three_prime_utrs
    strand = any_utr[0].strand if any_utr else None

    cds_min_genomic: int | None = None
    cds_max_genomic: int | None = None

    if strand == 1:
        if five_prime_utrs:
            cds_min_genomic = max(u.end for u in five_prime_utrs) + 1
        if three_prime_utrs:
            cds_max_genomic = min(u.start for u in three_prime_utrs) - 1
    elif strand == -1:
        if three_prime_utrs:
            cds_min_genomic = max(u.end for u in three_prime_utrs) + 1
        if five_prime_utrs:
            cds_max_genomic = min(u.start for u in five_prime_utrs) - 1

    return cds_min_genomic, cds_max_genomic


def _coding_interval(
    exon: ExonRecord, cds_min_genomic: int | None, cds_max_genomic: int | None
) -> tuple[int, int] | None:
    start, end = exon.start, exon.end
    if cds_min_genomic is not None:
        start = max(start, cds_min_genomic)
    if cds_max_genomic is not None:
        end = min(end, cds_max_genomic)
    if start > end:
        return None  # exon is entirely outside the CDS (e.g. a pure-UTR exon)
    return start, end


@dataclass
class ExonProteinBoundary:
    exon_rank: int
    start_aa: int
    end_aa: int


def exon_protein_boundaries(
    exons: list[ExonRecord],
    *,
    cds_min_genomic: int | None = None,
    cds_max_genomic: int | None = None,
) -> list[ExonProteinBoundary]:
    """Convert a transcript's exon structure into per-exon protein-residue
    ranges, ordered 5'->3' by rank.

    Reuses the same coding-interval clipping (:func:`_coding_interval`) and
    nucleotide-to-codon rounding (``ceil(nt / 3)``) arithmetic as
    :func:`map_genomic_breakpoint_to_protein_position`, so an exon boundary
    reported here is consistent with any individual breakpoint mapped
    against the same transcript. An exon entirely outside the CDS (e.g. a
    pure-UTR exon) contributes no boundary, same as that function's
    per-exon coding-interval clipping.
    """
    ordered = sorted(exons, key=lambda e: e.rank)
    boundaries: list[ExonProteinBoundary] = []
    cds_nt_cursor = 0
    for exon in ordered:
        interval = _coding_interval(exon, cds_min_genomic, cds_max_genomic)
        if interval is None:
            continue
        start, end = interval
        length_nt = end - start + 1
        start_aa = math.ceil((cds_nt_cursor + 1) / 3)
        cds_nt_cursor += length_nt
        end_aa = math.ceil(cds_nt_cursor / 3)
        boundaries.append(
            ExonProteinBoundary(exon_rank=exon.rank, start_aa=start_aa, end_aa=end_aa)
        )
    return boundaries


def gene_track_from_canonical_transcript(canonical: CanonicalTranscript) -> dict:
    """Build a JSON-serializable, gene-agnostic protein-layout summary
    (protein length, full domain map, exon-to-residue boundaries) from an
    already-fetched :class:`CanonicalTranscript`.

    Every field here comes from the same canonical-transcript payload
    already fetched for domain resolution elsewhere in the pipeline (see
    :func:`resolve_domains`); this is purely a different projection of
    that same response, not a new data source.
    """
    cds_min_genomic, cds_max_genomic = cds_bounds_from_utrs(canonical.utrs)
    domains = pfam_domains_to_protein_domains(canonical)
    boundaries = exon_protein_boundaries(
        canonical.exons, cds_min_genomic=cds_min_genomic, cds_max_genomic=cds_max_genomic
    )
    return {
        "protein_length": canonical.protein_length,
        "domains": [
            {
                "name": domain.name,
                "accession": domain.accession,
                "start_aa": domain.start_aa,
                "end_aa": domain.end_aa,
            }
            for domain in domains
        ],
        "exon_boundaries_aa": [
            {"exon_rank": b.exon_rank, "start_aa": b.start_aa, "end_aa": b.end_aa}
            for b in boundaries
        ],
    }


def map_genomic_breakpoint_to_protein_position(
    exons: list[ExonRecord],
    breakpoint_genomic: int,
    strand: int,
    *,
    cds_min_genomic: int | None = None,
    cds_max_genomic: int | None = None,
    role: str | None = None,
) -> BreakpointProteinPosition:
    """Estimate the protein residue a genomic breakpoint falls in.

    This is standard breakpoint-to-CDS-position arithmetic, not something
    Genome Nexus (or any variant-annotation API) computes for a fusion
    breakpoint directly:

    1. Order exons 5'->3' along the mature transcript by ``rank`` ascending
       (transcript order, independent of genomic strand).
    2. Optionally clip each exon to ``[cds_min_genomic, cds_max_genomic]``
       so a first/last exon that partly overlaps the 5'/3' UTR only
       contributes its coding portion (an exon entirely outside this
       range contributes nothing).
    3. Find the exon whose clipped genomic range contains the breakpoint.
       If none does (an intronic breakpoint), disambiguate directionally:

       - When ``role`` is ``"five_prime"`` (this gene retains everything
         *before* the breakpoint in the chimeric transcript -- everything
         after is spliced away), snap to the nearest exon strictly
         upstream of the breakpoint in transcript order -- the last
         fully-retained exon -- never a downstream exon just because it
         happens to be geometrically closer.
       - When ``role`` is ``"three_prime"`` (this gene retains everything
         *after* the breakpoint), snap to the nearest exon strictly
         downstream of the breakpoint in transcript order.
       - When ``role`` is ``None``/unrecognized (no fusion-partner context
         available, e.g. a same-gene intragenic event, or a caller that
         predates this parameter), fall back to whichever exon boundary is
         nearest by raw genomic distance, regardless of side -- the
         original, direction-unaware approximation.

       "Upstream"/"downstream" here means transcript order, which flips
       relative to genomic coordinate depending on ``strand``: on the plus
       strand transcript order tracks increasing genomic position, on the
       minus strand it runs in decreasing genomic position. Either way,
       the breakpoint is clamped to the chosen exon's nearest boundary --
       an approximation, flagged via ``is_intronic``.
    4. Convert the breakpoint's position within that exon to a
       CDS-relative nucleotide offset: on the plus strand this counts up
       from the exon's clipped start; on the minus strand -- transcribed
       high-to-low in genomic coordinate space -- it counts down from the
       exon's clipped end.
    5. The 1-indexed CDS nucleotide position is the sum of the coding
       lengths of every exon ranked before this one, plus that offset,
       plus one.
    6. The protein residue number is that nucleotide position divided by
       3, rounded up (any nucleotide within a codon's 3 bases maps to
       that codon's amino acid).
    """
    ordered = sorted(exons, key=lambda e: e.rank)
    coding_intervals = {
        exon.rank: _coding_interval(exon, cds_min_genomic, cds_max_genomic) for exon in ordered
    }

    containing = None
    for exon in ordered:
        interval = coding_intervals[exon.rank]
        if interval is None:
            continue
        start, end = interval
        if start <= breakpoint_genomic <= end:
            containing = exon
            break

    is_intronic = containing is None
    if containing is None:

        def _distance(exon: ExonRecord) -> float:
            interval = coding_intervals[exon.rank]
            if interval is None:
                return math.inf
            start, end = interval
            if breakpoint_genomic < start:
                return start - breakpoint_genomic
            if breakpoint_genomic > end:
                return breakpoint_genomic - end
            return 0

        def _is_upstream_of_breakpoint(exon: ExonRecord) -> bool:
            """True if ``exon`` lies strictly before the breakpoint in
            transcript order (i.e. is retained by a 5'-partner gene)."""
            interval = coding_intervals[exon.rank]
            if interval is None:
                return False
            start, end = interval
            return end <= breakpoint_genomic if strand == 1 else start >= breakpoint_genomic

        directional_candidates: list[ExonRecord] = []
        if role == "five_prime" or role == "three_prime":
            want_upstream = role == "five_prime"
            directional_candidates = [
                exon
                for exon in ordered
                if coding_intervals[exon.rank] is not None
                and _is_upstream_of_breakpoint(exon) == want_upstream
            ]

        containing = min(directional_candidates or ordered, key=_distance)
        start, end = coding_intervals[containing.rank]
        clamped_breakpoint = min(max(breakpoint_genomic, start), end)
    else:
        start, end = coding_intervals[containing.rank]
        clamped_breakpoint = breakpoint_genomic

    offset = (clamped_breakpoint - start) if strand == 1 else (end - clamped_breakpoint)

    preceding_length = sum(
        (coding_intervals[exon.rank][1] - coding_intervals[exon.rank][0] + 1)
        for exon in ordered
        if exon.rank < containing.rank and coding_intervals[exon.rank] is not None
    )

    cds_nt_position = preceding_length + offset + 1
    protein_position = math.ceil(cds_nt_position / 3)

    return BreakpointProteinPosition(
        protein_position=protein_position,
        exon_rank=containing.rank,
        is_intronic=is_intronic,
    )
