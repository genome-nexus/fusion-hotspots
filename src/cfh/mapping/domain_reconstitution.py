"""Does the fusion partner supply the rest of a truncated domain?

When a breakpoint falls inside a gene's key domain, the chimeric protein
keeps only part of it. If the partner carries the same Pfam family and its
breakpoint keeps the complementary part, the two pieces can join into one
hybrid domain. ETS-ETS fusions are the motivating case: in FLI1::ERG the
5' FLI1 keeps the N-terminal 82% of its ETS domain and the 3' ERG keeps the
C-terminal 28% of its own, together spanning the whole domain.

Both pieces are placed in domain-relative coordinates (0 = domain start,
1 = domain end) so domains of different absolute positions line up. That
assumes the two copies align end to end, which holds for paralogs whose
Pfam hits span the full model; it is an approximation otherwise.

What this does not establish:

- Frame. A complementary pair only yields a hybrid domain if the junction
  is in frame. The caller reports the frame call alongside.
- The join itself. Breakpoints here come from DNA. A mid-exon breakpoint
  cannot be spliced as called, so the real junction in the transcript may
  differ; a junction inside sequence the two paralogs share is also
  ambiguous in position, and complements by construction wherever it is
  placed. Mid-exon flags are returned so readers can weigh this.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, fields, replace

from cfh.mapping.genome_nexus_source import PfamDomain

DEFAULT_TOLERANCE_FRACTION = 0.10
"""Gap or overlap, as a fraction of the domain, still read as complementary.
0.10 is about 8 residues of an 83-residue ETS domain: room for breakpoint
imprecision and Pfam boundary differences between paralogs."""

COMPLEMENTARY = "complementary"
OVERLAPPING = "overlapping"
GAP = "gap"
PARTNER_DOMAIN_INTACT = "partner_domain_intact"
PARTNER_DOMAIN_LOST = "partner_domain_lost"
PARTNER_LACKS_DOMAIN = "partner_lacks_domain"
PARTNER_UNMAPPED = "partner_unmapped"
TARGET_DOMAIN_NOT_IN_PFAM = "target_domain_not_in_pfam"


@dataclass(frozen=True)
class DomainReconstitution:
    status: str
    pfam_accession: str | None = None
    five_prime_span: tuple[float, float] | None = None
    """Domain-relative interval kept by the 5' gene."""
    three_prime_span: tuple[float, float] | None = None
    """Domain-relative interval kept by the 3' gene."""
    junction_gap_fraction: float | None = None
    """Start of the 3' piece minus end of the 5' piece; negative = overlap."""
    junction_gap_residues: float | None = None
    """``junction_gap_fraction`` scaled by the mean of the two domain lengths."""
    partner_breakpoint_protein_position: int | None = None
    partner_domain_start_aa: int | None = None
    partner_domain_end_aa: int | None = None
    partner_breakpoint_mid_exon: bool | None = None
    detail: str | None = None


def reconstitution_columns(result: DomainReconstitution | None) -> dict:
    """Flat ``reconstitution_*`` event-row columns, present (blank) on every
    row so TSV output keeps one column set; spans become ``start-end``."""
    values = asdict(result) if result is not None else {}
    columns = {}
    for item in fields(DomainReconstitution):
        if item.name == "status":
            continue
        value = values.get(item.name)
        if isinstance(value, tuple):
            value = f"{value[0]:.4f}-{value[1]:.4f}"
        columns[f"reconstitution_{item.name}"] = value
    return columns


def kept_span(domain: PfamDomain, breakpoint_aa: int, role: str) -> tuple[float, float] | None:
    """Domain-relative interval a gene keeps, or ``None`` if it keeps none.

    Matches :func:`cfh.mapping.feature_mapper.calculate_domain_retention`:
    one-based inclusive coordinates, and the breakpoint residue stays with
    the surviving fragment.
    """
    length = domain.end_aa - domain.start_aa + 1
    if role == "five_prime":
        end = min(domain.end_aa, breakpoint_aa)
        if end < domain.start_aa:
            return None
        return (0.0, (end - domain.start_aa + 1) / length)
    if role == "three_prime":
        start = max(domain.start_aa, breakpoint_aa)
        if start > domain.end_aa:
            return None
        return ((start - domain.start_aa) / length, 1.0)
    raise ValueError(f"role must be 'five_prime' or 'three_prime', not {role!r}")


def pfam_hit_for_interval(
    pfam_domains: list[PfamDomain], start_aa: int | None, end_aa: int | None
) -> PfamDomain | None:
    """The Pfam hit overlapping a key-domain interval the most.

    Curated key domains may be UniProt features with no Pfam accession; the
    overlapping Pfam hit names the family to look for in the partner.
    """
    if start_aa is None or end_aa is None:
        return None
    best: tuple[int, PfamDomain] | None = None
    for domain in pfam_domains:
        overlap = min(end_aa, domain.end_aa) - max(start_aa, domain.start_aa) + 1
        if overlap > 0 and (best is None or overlap > best[0]):
            best = (overlap, domain)
    return best[1] if best else None


def assess_reconstitution(
    *,
    target_domain: PfamDomain,
    target_breakpoint_aa: int,
    target_role: str,
    partner_pfam_domains: list[PfamDomain],
    partner_breakpoint_aa: int,
    partner_breakpoint_mid_exon: bool | None = None,
    tolerance_fraction: float = DEFAULT_TOLERANCE_FRACTION,
) -> DomainReconstitution:
    """Compare what the target and partner each keep of one Pfam family."""
    partner_role = "three_prime" if target_role == "five_prime" else "five_prime"
    target_span = kept_span(target_domain, target_breakpoint_aa, target_role)
    if target_span is None:
        raise ValueError("target keeps none of its domain; it is lost, not truncated")
    base = DomainReconstitution(
        status=PARTNER_LACKS_DOMAIN,
        pfam_accession=target_domain.pfam_id,
        partner_breakpoint_protein_position=partner_breakpoint_aa,
        partner_breakpoint_mid_exon=partner_breakpoint_mid_exon,
    )
    candidates = [d for d in partner_pfam_domains if d.pfam_id == target_domain.pfam_id]
    if not candidates:
        return base

    def ordered(partner_span: tuple[float, float]) -> tuple[tuple[float, float], ...]:
        """(5' piece, 3' piece)."""
        if target_role == "five_prime":
            return (target_span, partner_span)
        return (partner_span, target_span)

    def distance(domain: PfamDomain) -> float:
        span = kept_span(domain, partner_breakpoint_aa, partner_role)
        if span is None:
            return float("inf")
        five, three = ordered(span)
        return abs(three[0] - five[1])

    # A partner with repeated copies of the family: judge the copy whose
    # kept piece sits closest to complementing the target's.
    domain = min(candidates, key=distance)
    base = replace(
        base, partner_domain_start_aa=domain.start_aa, partner_domain_end_aa=domain.end_aa
    )
    partner_span = kept_span(domain, partner_breakpoint_aa, partner_role)
    if partner_span is None:
        return replace(base, status=PARTNER_DOMAIN_LOST)
    five, three = ordered(partner_span)
    gap = three[0] - five[1]
    mean_length = (
        (target_domain.end_aa - target_domain.start_aa + 1) + (domain.end_aa - domain.start_aa + 1)
    ) / 2
    if partner_span == (0.0, 1.0):
        status = PARTNER_DOMAIN_INTACT
    elif gap > tolerance_fraction:
        status = GAP
    elif gap < -tolerance_fraction:
        status = OVERLAPPING
    else:
        status = COMPLEMENTARY
    return replace(
        base,
        status=status,
        five_prime_span=five,
        three_prime_span=three,
        junction_gap_fraction=round(gap, 4),
        junction_gap_residues=round(gap * mean_length, 1),
    )
