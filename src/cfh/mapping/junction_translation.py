"""Predict the chimeric CDS and protein a fusion's two breakpoints produce.

This module is pure and source-independent. It does no I/O and knows nothing
about MSK-IMPACT, cBioPortal or Genome Nexus: a caller builds two
:class:`TranscriptModel` objects from any source and passes each partner's
genomic breakpoint with an explicit genome build. Nothing here defaults a
build, a transcript or a cohort.

Coordinate convention
    Genomic positions are one-based and inclusive. A reported breakpoint is the
    last base kept on the 5' side and the first base kept on the 3' side, so
    both reported bases stay in the product
    (:data:`CONVENTION_BOTH_BASES_RETAINED`, the default).
    :data:`CONVENTION_THREE_PRIME_BASE_DROPPED` removes the 3' reported base
    instead: the 3' side starts one base further along the 3' transcript
    (strand-aware), and that position is located and joined exactly like a
    reported one, so the convention applies to every join, exon edges
    included. SV callers do not agree on this, and microhomology at the
    junction makes the true position ambiguous, so the convention is always
    recorded on the result: moving one base changes the frame.

How a breakpoint becomes a product
    - Both breaks inside coding exons, exon edges included: the DNA joins
      directly and the two exon pieces form one chimeric exon
      (``mid_exon_dna``). An exact last-base/first-base join is still a DNA
      join of those two bases, not a splice.
    - Both breaks in introns: the 5' donor splices to the 3' acceptor
      (``splice_mediated``).
    - One break intronic and the other inside a coding exon: the broken exon
      has lost its splice site, so splicing can only use the nearest intact
      exon past it. The exon is dropped and recorded in ``notes``; when no such
      exon exists there is no credible product.
    - An intronic break inside a splice site it relies on, i.e. within
      :data:`DONOR_SITE_NT` nt after the 5' partner's donor exon or within
      :data:`ACCEPTOR_SITE_NT` nt before the 3' partner's acceptor exon:
      the joined sequence there comes from the other partner, which this
      module does not have, so whether the site survives is unknown. The
      event gets no product and ``splice_site_unresolved``, with the reason.

    Intronic breakpoints are never silently moved to an exon edge. The exon
    used and the intronic distance are in the breakpoint's
    :class:`BreakpointLocation` and in ``notes``.

A break outside the coding region gives no chimeric reading frame, and is
reported as having no credible product with the reason.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace

CONVENTION_BOTH_BASES_RETAINED = "both_breakpoint_bases_retained"
CONVENTION_THREE_PRIME_BASE_DROPPED = "three_prime_breakpoint_base_dropped"

FRAME_IN = "in_frame"
FRAME_OUT = "out_of_frame"
FRAME_UNKNOWN = "unknown"

MECHANISM_SPLICE = "splice_mediated"
MECHANISM_MID_EXON = "mid_exon_dna"
MECHANISM_NONE = "no_credible_product"
MECHANISM_SPLICE_SITE_UNRESOLVED = "splice_site_unresolved"
"""An intronic break falls inside the splice site the product would need, so
no product is predicted until the joined genomic sequence is checked."""

DONOR_SITE_NT = 6
"""Intronic bases after a donor exon (+1..+6, including the GT) that a 5'
break must keep for the donor to be assumed intact."""
ACCEPTOR_SITE_NT = 20
"""Intronic bases before an acceptor exon that a 3' break must keep: the AG
(-1, -2), -3, and the polypyrimidine tract and branch-point side."""

REGION_CODING_EXON = "coding_exon"
REGION_FIVE_PRIME_UTR = "five_prime_utr"
REGION_THREE_PRIME_UTR = "three_prime_utr"
REGION_INTRON = "intron"
REGION_OUTSIDE = "outside_transcript"

DOMAIN_INTACT = "intact"
DOMAIN_TRUNCATED = "truncated"
DOMAIN_LOST = "lost"

JUNCTION_CONTEXT_RESIDUES = 15

_GENETIC_CODE = {
    a + b + c: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"[i]
    for i, (a, b, c) in enumerate((a, b, c) for a in "TCAG" for b in "TCAG" for c in "TCAG")
}

_BUILD_ALIASES = {"hg19": "GRCh37", "grch37": "GRCh37", "hg38": "GRCh38", "grch38": "GRCh38"}


def normalize_genome_build(build: str | None) -> str | None:
    """Canonical build spelling, or ``None`` when no build was supplied."""
    if build is None or not str(build).strip():
        return None
    text = str(build).strip()
    return _BUILD_ALIASES.get(text.lower(), text)


def translate(cds: str) -> str:
    """Translate whole codons of ``cds``; unknown bases give ``X``."""
    return "".join(
        _GENETIC_CODE.get(cds[i : i + 3].upper(), "X") for i in range(0, len(cds) - 2, 3)
    )


@dataclass(frozen=True)
class TranscriptExon:
    rank: int
    start: int
    end: int


@dataclass(frozen=True)
class TranscriptDomain:
    family: str
    start_aa: int
    end_aa: int


@dataclass(frozen=True)
class TranscriptModel:
    """One transcript's structure, from any source.

    ``cds_start`` and ``cds_end`` are the genomic bounds of the coding region
    (low and high coordinate, regardless of strand), including the stop codon.
    ``cds_sequence`` is the spliced coding sequence in transcript orientation.
    """

    gene: str
    transcript_id: str
    genome_build: str
    strand: int
    exons: tuple[TranscriptExon, ...]
    cds_start: int
    cds_end: int
    cds_sequence: str
    domains: tuple[TranscriptDomain, ...] = ()

    @property
    def ordered_exons(self) -> tuple[TranscriptExon, ...]:
        return tuple(sorted(self.exons, key=lambda exon: exon.rank))

    @property
    def span(self) -> tuple[int, int]:
        return (min(e.start for e in self.exons), max(e.end for e in self.exons))

    @property
    def protein(self) -> str:
        return translate(self.cds_sequence).rstrip("*")


@dataclass(frozen=True)
class _CodingSegment:
    rank: int
    start: int
    end: int
    cds_first: int
    cds_last: int


def _coding_segments(model: TranscriptModel) -> list[_CodingSegment]:
    """Each exon's coding part with its one-based CDS interval, 5' to 3'."""
    segments: list[_CodingSegment] = []
    cursor = 0
    for exon in model.ordered_exons:
        start, end = max(exon.start, model.cds_start), min(exon.end, model.cds_end)
        if start > end:
            continue
        length = end - start + 1
        segments.append(_CodingSegment(exon.rank, start, end, cursor + 1, cursor + length))
        cursor += length
    return segments


def transcript_model_problem(model: TranscriptModel) -> str | None:
    """Why ``model`` cannot be used for translation, or ``None`` if it can."""
    if not model.exons:
        return f"{model.transcript_id} has no exon structure"
    if model.strand not in (1, -1):
        return f"{model.transcript_id} has invalid strand {model.strand!r}"
    if not model.cds_sequence:
        return f"{model.transcript_id} has no CDS sequence"
    if len(model.cds_sequence) % 3:
        return f"{model.transcript_id} CDS length {len(model.cds_sequence)} is not a whole codon"
    segments = _coding_segments(model)
    length = segments[-1].cds_last if segments else 0
    if length != len(model.cds_sequence):
        return (
            f"{model.transcript_id} exon structure gives a {length} nt CDS but its sequence "
            f"is {len(model.cds_sequence)} nt"
        )
    return None


@dataclass(frozen=True)
class BreakpointLocation:
    """Where a breakpoint falls on one partner's transcript."""

    gene: str
    transcript_id: str
    genomic_position: int
    region: str
    exon_rank: int | None = None
    """The exon containing the breakpoint, for an exonic break."""
    cds_position: int | None = None
    """One-based CDS position, for a break in a coding exon."""
    at_exon_edge: str | None = None
    """``"first_base"`` or ``"last_base"`` of the exon, in transcript order."""
    upstream_exon_rank: int | None = None
    downstream_exon_rank: int | None = None
    distance_from_upstream_exon_nt: int | None = None
    """For an intronic break: nt past the end of the upstream exon."""
    distance_to_downstream_exon_nt: int | None = None
    """For an intronic break: nt before the start of the downstream exon."""
    describes: str = ""


def locate_breakpoint(model: TranscriptModel, position: int) -> BreakpointLocation:
    """Place a genomic position on ``model`` in transcript orientation."""
    base = BreakpointLocation(
        gene=model.gene,
        transcript_id=model.transcript_id,
        genomic_position=position,
        region=REGION_OUTSIDE,
    )
    exons = model.ordered_exons
    span_start, span_end = model.span
    if position < span_start or position > span_end:
        return replace(
            base,
            describes=f"outside {model.transcript_id} ({span_start}-{span_end})",
        )

    segments = {segment.rank: segment for segment in _coding_segments(model)}
    for exon in exons:
        if not exon.start <= position <= exon.end:
            continue
        low_edge, high_edge = position == exon.start, position == exon.end
        first_base = low_edge if model.strand == 1 else high_edge
        last_base = high_edge if model.strand == 1 else low_edge
        edge = "first_base" if first_base else "last_base" if last_base else None
        segment = segments.get(exon.rank)
        if segment is not None and segment.start <= position <= segment.end:
            offset = position - segment.start if model.strand == 1 else segment.end - position
            cds_position = segment.cds_first + offset
            return replace(
                base,
                region=REGION_CODING_EXON,
                exon_rank=exon.rank,
                cds_position=cds_position,
                at_exon_edge=edge,
                describes=f"exon {exon.rank}, c.{cds_position}",
            )
        before_cds = (position < model.cds_start) == (model.strand == 1)
        region = REGION_FIVE_PRIME_UTR if before_cds else REGION_THREE_PRIME_UTR
        label = "5' UTR" if before_cds else "3' UTR"
        return replace(
            base,
            region=region,
            exon_rank=exon.rank,
            at_exon_edge=edge,
            describes=f"exon {exon.rank} ({label})",
        )

    upstream = next(
        exon
        for exon in reversed(exons)
        if (exon.end < position if model.strand == 1 else exon.start > position)
    )
    downstream = next(
        exon
        for exon in exons
        if (exon.start > position if model.strand == 1 else exon.end < position)
    )
    from_upstream = position - upstream.end if model.strand == 1 else upstream.start - position
    to_downstream = downstream.start - position if model.strand == 1 else position - downstream.end
    return replace(
        base,
        region=REGION_INTRON,
        upstream_exon_rank=upstream.rank,
        downstream_exon_rank=downstream.rank,
        distance_from_upstream_exon_nt=from_upstream,
        distance_to_downstream_exon_nt=to_downstream,
        describes=(
            f"intron between exons {upstream.rank} and {downstream.rank}, {from_upstream} nt "
            f"past exon {upstream.rank}, {to_downstream} nt before exon {downstream.rank}"
        ),
    )


@dataclass(frozen=True)
class DomainSpan:
    """What one partner keeps of one domain, in residues of its own protein."""

    partner: str
    gene: str
    family: str
    start_aa: int
    end_aa: int
    retained_start_aa: int | None
    retained_end_aa: int | None
    status: str


@dataclass(frozen=True)
class JunctionPrediction:
    five_prime_gene: str
    three_prime_gene: str
    genome_build: str | None
    coordinate_convention: str
    five_prime_location: BreakpointLocation | None = None
    three_prime_location: BreakpointLocation | None = None
    mechanism: str = MECHANISM_NONE
    frame_status: str = FRAME_UNKNOWN
    no_product_reason: str | None = None
    five_prime_cds_retained: int | None = None
    """Number of 5' CDS bases kept (c.1 through this position)."""
    three_prime_cds_start: int | None = None
    """One-based CDS position of the first 3' base kept."""
    five_prime_last_residue: int | None = None
    """Last 5' residue encoded by a complete 5' codon."""
    three_prime_first_residue: int | None = None
    """First 3' residue encoded by a complete 3' codon."""
    hybrid_codon: str | None = None
    """The junction-spanning codon, when the 5' side ends mid-codon and in frame."""
    chimeric_cds: str | None = None
    chimeric_protein: str | None = None
    """Translation from the 5' start through the first stop codon, without it."""
    junction_cds: str | None = None
    junction_peptide: str | None = None
    """Residues around the junction: 5' side, ``|``, 3' side."""
    premature_stop_at_junction: bool = False
    domain_spans: tuple[DomainSpan, ...] = ()
    notes: tuple[str, ...] = field(default=())

    @property
    def has_product(self) -> bool:
        return self.chimeric_cds is not None


def _domain_spans(
    model: TranscriptModel, partner: str, last_residue: int | None, first_residue: int | None
) -> tuple[DomainSpan, ...]:
    spans: list[DomainSpan] = []
    for domain in model.domains:
        if partner == "five_prime":
            kept_start, kept_end = domain.start_aa, min(domain.end_aa, last_residue or 0)
        else:
            kept_start, kept_end = max(domain.start_aa, first_residue or 0), domain.end_aa
        if kept_end < kept_start:
            status, kept_start_value, kept_end_value = DOMAIN_LOST, None, None
        else:
            whole = kept_start == domain.start_aa and kept_end == domain.end_aa
            status = DOMAIN_INTACT if whole else DOMAIN_TRUNCATED
            kept_start_value, kept_end_value = kept_start, kept_end
        spans.append(
            DomainSpan(
                partner,
                model.gene,
                domain.family,
                domain.start_aa,
                domain.end_aa,
                kept_start_value,
                kept_end_value,
                status,
            )
        )
    return tuple(spans)


def _five_prime_candidates(
    model: TranscriptModel, location: BreakpointLocation
) -> tuple[int | None, bool, str | None, str | None]:
    """``(cds_bases_kept, join_is_direct, note, problem)`` for the 5' partner.

    ``join_is_direct`` means the break leaves a piece of a coding exon that
    can only continue by DNA, not by its own splice donor.
    """
    segments = _coding_segments(model)
    if location.region == REGION_CODING_EXON:
        direct = location.at_exon_edge != "last_base"
        return location.cds_position, direct, None, None
    if location.region == REGION_INTRON:
        upstream_rank = location.upstream_exon_rank
        assert upstream_rank is not None  # set for every intronic location
        by_rank = {segment.rank: segment for segment in segments}
        upstream = by_rank.get(upstream_rank)
        if upstream is None:
            if segments and upstream_rank > segments[-1].rank:
                return None, False, None, f"5' break in {model.gene} is past its stop codon"
            return None, False, None, f"5' break in {model.gene} is before its coding region"
        return (
            upstream.cds_last,
            False,
            (
                f"{model.gene} 5' break is intronic ({location.describes}); product uses "
                f"exon {upstream.rank}'s splice donor, not the break position"
            ),
            None,
        )
    if location.region == REGION_THREE_PRIME_UTR:
        return None, False, None, f"5' break in {model.gene} is past its stop codon"
    if location.region == REGION_FIVE_PRIME_UTR:
        return None, False, None, f"5' break in {model.gene} is before its coding region"
    return None, False, None, f"5' break in {model.gene} is outside the transcript"


def _three_prime_candidates(
    model: TranscriptModel, location: BreakpointLocation
) -> tuple[int | None, bool, str | None, str | None]:
    """``(first_cds_base_kept, join_is_direct, note, problem)`` for the 3' partner."""
    segments = _coding_segments(model)
    if location.region == REGION_CODING_EXON:
        direct = location.at_exon_edge != "first_base"
        return location.cds_position, direct, None, None
    if location.region == REGION_INTRON:
        downstream_rank = location.downstream_exon_rank
        assert downstream_rank is not None  # set for every intronic location
        by_rank = {segment.rank: segment for segment in segments}
        downstream = by_rank.get(downstream_rank)
        if downstream is None:
            if segments and downstream_rank < segments[0].rank:
                return None, False, None, f"3' break in {model.gene} is before its coding region"
            return None, False, None, f"3' break in {model.gene} is past its stop codon"
        return (
            downstream.cds_first,
            False,
            (
                f"{model.gene} 3' break is intronic ({location.describes}); product uses "
                f"exon {downstream.rank}'s splice acceptor, not the break position"
            ),
            None,
        )
    if location.region == REGION_FIVE_PRIME_UTR:
        return None, False, None, f"3' break in {model.gene} is before its coding region"
    if location.region == REGION_THREE_PRIME_UTR:
        return None, False, None, f"3' break in {model.gene} is past its stop codon"
    return None, False, None, f"3' break in {model.gene} is outside the transcript"


def _splice_five_prime_past_broken_exon(
    model: TranscriptModel, location: BreakpointLocation
) -> tuple[int | None, str]:
    """Nearest intact 5' donor upstream of a broken coding exon."""
    broken_rank = location.exon_rank
    assert broken_rank is not None  # set for every exonic location
    upstream = [s for s in _coding_segments(model) if s.rank < broken_rank]
    if not upstream:
        return None, (
            f"{model.gene} 5' break is inside its first coding exon "
            f"({location.describes}), so no upstream splice donor can be used"
        )
    segment = upstream[-1]
    return segment.cds_last, (
        f"{model.gene} 5' break is inside exon {location.exon_rank} "
        f"({location.describes}); its donor is lost, so the product uses exon "
        f"{segment.rank}'s donor and drops the exon {location.exon_rank} piece"
    )


def _splice_three_prime_past_broken_exon(
    model: TranscriptModel, location: BreakpointLocation
) -> tuple[int | None, str]:
    """Nearest intact 3' acceptor downstream of a broken coding exon."""
    broken_rank = location.exon_rank
    assert broken_rank is not None  # set for every exonic location
    downstream = [s for s in _coding_segments(model) if s.rank > broken_rank]
    if not downstream:
        return None, (
            f"{model.gene} 3' break is inside its last coding exon "
            f"({location.describes}), so no downstream splice acceptor exists"
        )
    segment = downstream[0]
    return segment.cds_first, (
        f"{model.gene} 3' break is inside exon {location.exon_rank} "
        f"({location.describes}); its acceptor is lost, so the product uses exon "
        f"{segment.rank}'s acceptor and drops the exon {location.exon_rank} piece"
    )


def _splice_site_problems(
    five_gene: str,
    five_location: BreakpointLocation,
    three_gene: str,
    three_location: BreakpointLocation,
) -> list[str]:
    """Why an intronic break may have destroyed the splice site it relies on.

    A 5' intronic break keeps intron bases +1..+d after its donor exon; the
    rest of the donor comes from the 3' partner's sequence. A 3' intronic
    break keeps bases -d..-1 before its acceptor exon. When ``d`` is short of
    :data:`DONOR_SITE_NT` or :data:`ACCEPTOR_SITE_NT` the site may or may not
    be rebuilt by the joined sequence, which is not available here.
    """
    problems: list[str] = []
    if five_location.region == REGION_INTRON:
        kept = five_location.distance_from_upstream_exon_nt
        assert kept is not None  # set for every intronic location
        if kept < DONOR_SITE_NT:
            problems.append(
                f"{five_gene} 5' break keeps only {kept} nt of the exon "
                f"{five_location.upstream_exon_rank} splice donor (+1..+{DONOR_SITE_NT}); "
                "whether the joined sequence restores it is unresolved"
            )
    if three_location.region == REGION_INTRON:
        kept = three_location.distance_to_downstream_exon_nt
        assert kept is not None  # set for every intronic location
        if kept < ACCEPTOR_SITE_NT:
            region = "the essential AG" if kept < 3 else "the polypyrimidine tract / branch region"
            problems.append(
                f"{three_gene} 3' break keeps only {kept} nt before the exon "
                f"{three_location.downstream_exon_rank} splice acceptor "
                f"(-{ACCEPTOR_SITE_NT}..-1), inside {region}; whether the joined "
                "sequence restores it is unresolved"
            )
    return problems


def _no_product(
    base: JunctionPrediction, reason: str, notes: list[str] | tuple[str, ...] = ()
) -> JunctionPrediction:
    return replace(base, mechanism=MECHANISM_NONE, no_product_reason=reason, notes=tuple(notes))


def predict_junction(
    five_prime_gene: str,
    five_prime_breakpoint: int,
    five_prime_transcript: TranscriptModel | None,
    three_prime_gene: str,
    three_prime_breakpoint: int,
    three_prime_transcript: TranscriptModel | None,
    genome_build: str | None,
    *,
    coordinate_convention: str = CONVENTION_BOTH_BASES_RETAINED,
) -> JunctionPrediction:
    """Predict the chimeric product of one fusion.

    ``genome_build`` is required and is compared with both transcript models.
    A missing or mismatched build, or an unusable transcript model, gives a
    prediction with ``no_product_reason`` set rather than an exception, so a
    cohort screen can report each event's reason.
    """
    if coordinate_convention not in (
        CONVENTION_BOTH_BASES_RETAINED,
        CONVENTION_THREE_PRIME_BASE_DROPPED,
    ):
        raise ValueError(f"unknown coordinate convention {coordinate_convention!r}")
    build = normalize_genome_build(genome_build)
    base = JunctionPrediction(
        five_prime_gene=five_prime_gene,
        three_prime_gene=three_prime_gene,
        genome_build=build,
        coordinate_convention=coordinate_convention,
    )
    if build is None:
        return _no_product(base, "the event has no genome build")
    for label, model in (("5'", five_prime_transcript), ("3'", three_prime_transcript)):
        if model is None:
            gene = five_prime_gene if label == "5'" else three_prime_gene
            return _no_product(base, f"no transcript model for {label} partner {gene}")
        if normalize_genome_build(model.genome_build) != build:
            return _no_product(
                base,
                f"{label} transcript {model.transcript_id} is {model.genome_build}, "
                f"the event is {build}",
            )
        problem = transcript_model_problem(model)
        if problem:
            return _no_product(base, problem)
    assert five_prime_transcript is not None and three_prime_transcript is not None

    convention_notes: list[str] = []
    three_position = three_prime_breakpoint
    if coordinate_convention == CONVENTION_THREE_PRIME_BASE_DROPPED:
        three_position += three_prime_transcript.strand
        convention_notes.append(
            f"{three_prime_gene} reported 3' base {three_prime_breakpoint} is dropped; "
            f"the 3' side starts at {three_position}"
        )
    five_location = locate_breakpoint(five_prime_transcript, five_prime_breakpoint)
    three_location = locate_breakpoint(three_prime_transcript, three_position)
    base = replace(base, five_prime_location=five_location, three_prime_location=three_location)

    five_kept, five_direct, five_note, five_problem = _five_prime_candidates(
        five_prime_transcript, five_location
    )
    three_first, three_direct, three_note, three_problem = _three_prime_candidates(
        three_prime_transcript, three_location
    )
    notes = [*convention_notes, *(note for note in (five_note, three_note) if note)]
    if five_problem or three_problem:
        return _no_product(base, "; ".join(p for p in (five_problem, three_problem) if p), notes)
    assert five_kept is not None and three_first is not None

    disrupted = _splice_site_problems(
        five_prime_gene, five_location, three_prime_gene, three_location
    )
    if disrupted:
        return replace(
            _no_product(base, "; ".join(disrupted), notes),
            mechanism=MECHANISM_SPLICE_SITE_UNRESOLVED,
        )

    five_exonic = five_location.region == REGION_CODING_EXON
    three_exonic = three_location.region == REGION_CODING_EXON
    if five_exonic and three_exonic:
        # Both breaks are in coding exons: the DNA joins directly, exon edges
        # included, so the coordinate convention moves the 3' start here too.
        mechanism = MECHANISM_MID_EXON
    elif five_direct:
        # 3' break is intronic, so the 5' exon piece has no donor of its own.
        rescued, note = _splice_five_prime_past_broken_exon(five_prime_transcript, five_location)
        notes.append(note)
        if rescued is None:
            return _no_product(base, note, notes)
        five_kept, mechanism = rescued, MECHANISM_SPLICE
    elif three_direct:
        # 5' break is intronic, so the 3' exon piece has no acceptor of its own.
        rescued, note = _splice_three_prime_past_broken_exon(three_prime_transcript, three_location)
        notes.append(note)
        if rescued is None:
            return _no_product(base, note, notes)
        three_first, mechanism = rescued, MECHANISM_SPLICE
    else:
        mechanism = MECHANISM_SPLICE

    frame = FRAME_IN if five_kept % 3 == (three_first - 1) % 3 else FRAME_OUT
    five_cds = five_prime_transcript.cds_sequence[:five_kept]
    three_cds = three_prime_transcript.cds_sequence[three_first - 1 :]
    chimeric_cds = five_cds + three_cds
    full_translation = translate(chimeric_cds)
    stop = full_translation.find("*")
    chimeric_protein = full_translation[:stop] if stop >= 0 else full_translation

    last_residue = five_kept // 3
    first_residue = (three_first - 1) // 3 + (1 if (three_first - 1) % 3 == 0 else 2)
    hybrid = five_kept % 3 != 0 and frame == FRAME_IN
    codon_start = five_kept - five_kept % 3
    hybrid_codon: str | None = chimeric_cds[codon_start : codon_start + 3] if hybrid else None

    junction_aa_index = five_kept // 3
    left = chimeric_protein[
        max(0, junction_aa_index - JUNCTION_CONTEXT_RESIDUES) : junction_aa_index
    ]
    if frame == FRAME_IN:
        right_start = junction_aa_index + (1 if hybrid else 0)
        hybrid_residue = _GENETIC_CODE.get(hybrid_codon.upper(), "X") if hybrid_codon else ""
        premature_stop = hybrid_residue == "*"
        right = chimeric_protein[right_start : right_start + JUNCTION_CONTEXT_RESIDUES]
        peptide = f"{left}{hybrid_residue}|{right}" if hybrid_residue else f"{left}|{right}"
    else:
        premature_stop = False
        right = chimeric_protein[junction_aa_index : junction_aa_index + JUNCTION_CONTEXT_RESIDUES]
        peptide = f"{left}|{right}"
    junction_cds = chimeric_cds[max(0, five_kept - 45) : five_kept + 45]

    return replace(
        base,
        mechanism=mechanism,
        frame_status=frame,
        five_prime_cds_retained=five_kept,
        three_prime_cds_start=three_first,
        five_prime_last_residue=last_residue,
        three_prime_first_residue=first_residue,
        hybrid_codon=hybrid_codon,
        chimeric_cds=chimeric_cds,
        chimeric_protein=chimeric_protein,
        junction_cds=junction_cds,
        junction_peptide=peptide,
        premature_stop_at_junction=premature_stop,
        domain_spans=(
            *_domain_spans(five_prime_transcript, "five_prime", last_residue, None),
            *_domain_spans(three_prime_transcript, "three_prime", None, first_residue),
        ),
        notes=tuple(notes),
    )
