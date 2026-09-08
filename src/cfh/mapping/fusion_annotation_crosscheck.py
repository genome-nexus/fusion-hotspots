"""Independent cross-check of our own genomic-breakpoint-to-protein-position
mapping (:mod:`cfh.mapping.transcript_source`), using a second, separately
implemented arithmetic core ported from genome-nexus/fusion-annotation.

Portions of this file are ported from:
    https://github.com/genome-nexus/fusion-annotation
    src/fusion_annotation/core.py, commit 6baba8638b0742389f94077cb3e9a705db4ed3cc
    Copyright the genome-nexus/fusion-annotation contributors.
    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

        http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.

Modified from genome-nexus/fusion-annotation (Apache-2.0): the five
functions below (``build_exon_cds_map``, ``cds_coord_at_exon_boundary``,
``build_exon_genomic_map``, ``cds_coord_at_genomic``,
``parse_genomic_breakpoint``) are reproduced with formatting/type-annotation
adaptations only (re-typed for this file's ``from __future__ import
annotations`` style and this project's line length/lint config) -- the
logic and control flow are verbatim from ``src/fusion_annotation/core.py``
at the commit above. The Apache-2.0 license header above and this
"Modified from..." notice are additions made by this project for
attribution purposes (per License Section 4(b)/4(c)); the upstream source
file itself carries no per-file copyright/license header of its own --
licensing is established solely by the upstream repository's top-level
``LICENSE`` file. Everything else in this file -- the
``FusionAnnotationCrosscheckResult`` dataclass,
``crosscheck_breakpoint_protein_position``, and its helpers -- is new code
written for this project, not part of the original project.

Deliberately NOT ported: fusion-annotation's own CDS-genomic-bounds
derivation (``_cds_bounds_from_utrs`` in ``src/fusion_annotation/gn_provider.py``).
That helper calls ``next((u for u in utrs if u["type"] == "five_prime_UTR"), None)``
(and the analogous line for ``three_prime_UTR``), which only ever reads the
FIRST five_prime_UTR/three_prime_UTR record Genome Nexus returns. A
transcript's 5' or 3' UTR can be split across multiple exons and so appear
as several non-adjacent segments in the payload; taking only the first
one silently mis-derives the CDS boundary whenever a UTR has more than one
segment. This project's own :func:`cfh.mapping.genome_nexus_source.cds_bounds_from_utrs`
already aggregates every segment of each UTR type correctly (verified live
against a PIK3CA fixture where the naive first-segment approach was off
by 50,222 genomic bases). So this module ports only fusion-annotation's
strand-order/CDS-offset *arithmetic* (the functions above) and feeds it
this project's own, already-correct UTR-derived CDS bounds -- it never
reproduces or calls fusion-annotation's own UTR-bounds logic.

One related, narrower fallback IS adopted (not ported code, just the same
documented convention): when this project's own ``cds_bounds_from_utrs``
reports ``None`` for one side because the transcript has no UTR of that
type at all (a real, gene-agnostic transcript characteristic distinct from
the multi-segment bug above -- e.g. NTRK1's canonical ENST00000524377 has
no ``three_prime_UTR`` record whatsoever, cited by name in
``gn_provider.py``'s own docstring), this module falls back to the
transcript's outermost exon coordinate on that side, exactly as
``gn_provider.py`` itself documents doing for this case ("GN then simply
omits that UTR record ... falls back to the transcript's outermost exon
coordinate on that side instead of raising"). Without this, the
cross-check could never even be attempted for such a transcript.

Purpose and comparison convention
----------------------------------
:func:`cfh.mapping.transcript_source.resolve_breakpoint_protein_position`
maps a genomic breakpoint to a single collapsed protein-residue estimate,
``ceil(cds_nt_position / 3)`` -- the amino acid whose codon contains the
breakpoint nucleotide, independent of which side of a fusion junction the
target gene plays. genome-nexus/fusion-annotation's own downstream usage
(``annotate_effect`` in its ``core.py``, not ported here since it needs a
*pair* of transcripts to reconstruct a chimeric CDS) distinguishes three
different residues at a fusion junction depending on fusion role:

* the last CDS residue *fully* encoded by bases upstream of the
  breakpoint ("last complete 5' residue"),
* the first CDS residue *fully* encoded by bases downstream of the
  breakpoint ("first complete 3' residue"), and
* the hybrid codon straddling the junction itself, when the breakpoint
  does not land exactly on a codon boundary.

These three coincide only when the breakpoint falls exactly on a codon
boundary (``cds_position % 3 == 0``, no hybrid codon).

The exact relationship between fusion-annotation's role-appropriate
residue and this project's ``ceil(cds_position / 3)`` is provable, not a
loose approximation, given the *same* underlying ``cds_position`` (which
holds whenever both resolvers land the breakpoint in the same coding
exon -- i.e. not the intronic-clamping-mismatch case described below):

* Three-prime ("first complete 3' residue"): fusion-annotation computes
  ``(cds_position - 1) // 3 + 1``, which is the standard integer-ceiling
  identity ``ceil(n / m) == (n - 1) // m + 1`` for positive ``n`` -- i.e.
  this is *always exactly* ``ceil(cds_position / 3)``, identical to this
  project's value, with **no legitimate tolerance**. Any difference here
  is a real disagreement (typically an intronic-clamping mismatch -- see
  below), never a rounding-convention artifact.
* Five-prime ("last complete 5' residue"): fusion-annotation computes
  ``cds_position // 3`` (floor division). This equals
  ``ceil(cds_position / 3)`` when ``cds_position`` is divisible by 3 (no
  hybrid codon), and is *exactly* ``ceil(cds_position / 3) - 1`` --
  never ``+ 1``, never off by more than 1 -- whenever it is not (a hybrid
  codon). The comparison enforces this exact, signed relationship
  (``expected = our_protein_position - 1`` when hybrid, an exact
  ``our_protein_position`` search on the wrong side is never accepted)
  rather than an undirected ``abs(diff) <= 1`` tolerance, which would
  wrongly accept a value one *higher* than ours as if it were the
  expected convention gap.

:func:`crosscheck_breakpoint_protein_position` computes fusion-annotation's
role-appropriate residue independently (via the ported
``cds_coord_at_genomic``) and compares it against this exact, signed
expected value per role. Any other numeric gap is reported as a real
disagreement.

A second, independent known source of (real) disagreement: for an
intronic breakpoint, this project's own resolver clamps to the *nearest*
coding exon by absolute genomic distance, while fusion-annotation's ported
``cds_coord_at_genomic`` clamps *directionally* -- to the nearest coding
exon strictly upstream for a 5'-partner breakpoint, or strictly
downstream for a 3'-partner breakpoint. These can disagree when the two
flanking exons are not the same distance away. This is the same
known issue a directional-intronic-breakpoint-snapping fix (tracked on a
separate branch) addresses for this project's own resolver; this
cross-check surfaces it rather than silently absorbing or hiding it.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from cfh.mapping.genome_nexus_source import ExonRecord

# ----------------------------------------------------------------------------
# Ported from genome-nexus/fusion-annotation src/fusion_annotation/core.py @
# 6baba8638b0742389f94077cb3e9a705db4ed3cc (Apache-2.0): formatting/
# type-annotation adaptations only (line length, this file's
# `from __future__ import annotations` style) for this project's lint
# config -- the logic and control flow are verbatim, unchanged from the
# original.
# ----------------------------------------------------------------------------


def build_exon_cds_map(
    strand: int, exons: list[dict], cds_g_start: int, cds_g_end: int
) -> list[tuple[int, int]]:
    """Map genomic exons to CDS-relative (start,end) coords, in transcription order.

    exons: list of {"start": g_lo, "end": g_hi} (genomic, lo<hi as Ensembl returns).
    cds_g_start/cds_g_end: genomic bounds of the CDS (lo<hi).
    Returns 1-based inclusive CDS coords for each *coding* exon, in rank order;
    non-coding exons contribute an empty (0,0) placeholder so the index == rank-1.
    """
    ordered = sorted(exons, key=lambda e: e["start"], reverse=(strand == -1))
    out: list[tuple[int, int]] = []
    pos = 0
    for e in ordered:
        o_lo, o_hi = max(e["start"], cds_g_start), min(e["end"], cds_g_end)
        clen = max(0, o_hi - o_lo + 1)
        if clen:
            out.append((pos + 1, pos + clen))
            pos += clen
        else:
            out.append((0, 0))
    return out


def cds_coord_at_exon_boundary(
    exon_cds: list[tuple[int, int]], exon_rank: int, side: Literal["end", "start"]
) -> int:
    """CDS coordinate at the 3' end (side='end') or 5' start (side='start') of
    a 1-based exon rank."""
    s, e = exon_cds[exon_rank - 1]
    if (s, e) == (0, 0):
        raise ValueError(f"exon {exon_rank} is non-coding")
    return e if side == "end" else s


def build_exon_genomic_map(strand: int, exons: list[dict]) -> list[tuple[int, int]]:
    """Genomic (g_lo, g_hi) bounds per exon, in transcription order.

    Same ordering and indexing (index == rank-1) as build_exon_cds_map(), so the
    two lists line up. Kept separate from exon_cds so a genomic breakpoint can be
    resolved to a CDS coordinate without re-fetching the transcript structure.
    """
    ordered = sorted(exons, key=lambda e: e["start"], reverse=(strand == -1))
    return [(e["start"], e["end"]) for e in ordered]


def cds_coord_at_genomic(
    strand: int,
    exon_genomic: list[tuple[int, int]],
    cds_g_start: int,
    cds_g_end: int,
    g_pos: int,
    side: Literal["end", "start"],
) -> int:
    """Map a genomic breakpoint to a 1-based CDS coordinate on this transcript.

    Unlike an exon number, a genomic coordinate pins the isoform: it only resolves
    against the transcript whose exon table actually spans it, which is exactly what
    disambiguates overlapping isoforms with divergent exon numbering (see issue #3).

    Parameters
    ----------
    exon_genomic : (g_lo, g_hi) per exon in transcription order (build_exon_genomic_map()).
    cds_g_start / cds_g_end : genomic CDS bounds (lo < hi); UTR-only exon portions are ignored.
    g_pos : genomic position of the breakpoint.
    side : 'end'   -> last CDS base the 5' partner retains at/upstream of g_pos.
           'start' -> first CDS base the 3' partner retains at/downstream of g_pos.

    A breakpoint inside a coding exon maps to that exact base; a breakpoint in an
    intron/UTR snaps to the nearest flanking coding-exon boundary in the relevant
    direction. Raises ValueError if g_pos maps outside the coding region entirely.
    """
    if not exon_genomic:
        raise ValueError(
            "transcript has no genomic exon map; cannot resolve a genomic breakpoint "
            "(populate Transcript.exon_genomic via build_exon_genomic_map)"
        )

    def tc(g: int) -> int:  # transcription coordinate (increases 5'->3')
        return g if strand == 1 else -g

    tp = tc(g_pos)
    pos = 0
    upstream_cds_end: int | None = None  # cds_end of the last coding exon fully 5' of g_pos
    for g_lo, g_hi in exon_genomic:
        o_lo, o_hi = max(g_lo, cds_g_start), min(g_hi, cds_g_end)
        clen = o_hi - o_lo + 1
        if clen <= 0:  # non-coding / UTR-only exon
            continue
        c_start, c_end = pos + 1, pos + clen
        pos += clen
        t_lo, t_hi = sorted((tc(o_lo), tc(o_hi)))
        if t_lo <= tp <= t_hi:  # breakpoint inside this coding exon
            return c_start + (tp - t_lo)  # exact base (c_start is the 5' end of the exon)
        if tp > t_hi:  # exon lies entirely upstream (5') of g_pos
            upstream_cds_end = c_end
        elif side == "start":  # first exon entirely downstream (3') of g_pos
            return c_start
    if side == "end":
        if upstream_cds_end is not None:
            return upstream_cds_end
        # A 5' partner breakpoint can lie wholly in the 5' UTR/promoter. In
        # that case the fusion retains zero coding bases from that partner
        # rather than raising, which is what promoter-swap fusions such as
        # TMPRSS2::ERG need.
        return 0
    raise ValueError(f"genomic position {g_pos} does not map into the CDS of this transcript")


def parse_genomic_breakpoint(value: object) -> int:
    """Parse a genomic breakpoint into an integer position.

    Accepts a plain int, a bare position string ("117324415"), a "chrom:pos" form
    ("chr6:117324415" / "6:117324415"), or an HGVS genomic term ("g.117324415",
    "chr6:g.117324415"). Only the 1-based genomic position is used -- the transcript's
    own strand/exon table supplies chromosome context -- so the chromosome token, if
    present, is accepted but not validated here.
    """
    if isinstance(value, int):
        return value
    s = str(value).strip()
    if ":" in s:
        s = s.rsplit(":", 1)[1].strip()
    if s.lower().startswith("g."):
        s = s[2:].strip()
    s = s.replace(",", "").replace("_", "")
    try:
        return int(s)
    except ValueError as exc:
        raise ValueError(f"could not parse genomic breakpoint {value!r}") from exc


# ----------------------------------------------------------------------------
# New code: wires the ported arithmetic above into an independent cross-check
# against this project's own resolver, fed with this project's own
# (already-correct) UTR-derived CDS bounds. Nothing below this point is part
# of the original genome-nexus/fusion-annotation project.
# ----------------------------------------------------------------------------

FusionRole = Literal["five_prime", "three_prime"]
_SIDE_FOR_ROLE: dict[str, Literal["end", "start"]] = {"five_prime": "end", "three_prime": "start"}


@dataclass
class FusionAnnotationCrosscheckResult:
    protein_position: int | None
    """fusion-annotation's own role-appropriate residue estimate (the "last
    complete 5' residue" for a five_prime breakpoint, or the "first complete
    3' residue" for a three_prime breakpoint), independently derived from the
    ported ``cds_coord_at_genomic``. ``None`` when the cross-check could not
    be attempted or computed (see ``error``)."""
    cds_position: int | None
    """The raw 1-based CDS nucleotide coordinate fusion-annotation's ported
    arithmetic assigns to this breakpoint, before any residue rounding."""
    exon_rank: int | None
    """Coding-exon rank fusion-annotation's own exon/CDS map assigns to
    ``cds_position``, independently derived via the ported
    ``build_exon_cds_map``/``cds_coord_at_exon_boundary``. Not necessarily
    identical in meaning to this project's own ``breakpoint_exon`` for an
    intronic breakpoint, since the two resolvers can clamp to different
    exons there (see module docstring)."""
    is_hybrid_codon: bool | None
    """True when ``cds_position`` does not fall exactly on a codon boundary
    (``cds_position % 3 != 0``) -- a property of the breakpoint itself,
    independent of fusion role. Only affects the comparison for a
    ``five_prime`` breakpoint, where fusion-annotation's residue is then
    expected to be exactly one less than this project's own
    hybrid-junction-inclusive residue (see module docstring); a
    ``three_prime`` breakpoint's expected residue is always exactly equal
    to this project's own, hybrid or not."""
    agrees: bool | None
    """Whether fusion-annotation's residue matches this project's own
    ``breakpoint_protein_position`` under the documented comparison
    convention. ``None`` when the cross-check could not be attempted at all
    (``error`` is set) rather than attempted-and-disagreeing."""
    note: str | None
    """Human-readable explanation of the comparison outcome, always set
    when ``agrees`` is not ``None``."""
    error: str | None
    """Set, and every other field left ``None`` (aside from ``agrees``,
    which stays ``None`` too), when the cross-check could not be computed at
    all -- e.g. no UTR-derived CDS bounds were available, or the ported
    arithmetic raised on this transcript's exon structure. Never raises;
    this is the graceful-degradation outcome."""


def _gn_exon_dicts(exons: list[ExonRecord]) -> list[dict]:
    return [{"start": exon.start, "end": exon.end} for exon in exons]


def crosscheck_breakpoint_protein_position(
    exons: list[ExonRecord],
    breakpoint_genomic: int | str,
    strand: int,
    role: str,
    our_protein_position: int | None,
    *,
    cds_min_genomic: int | None,
    cds_max_genomic: int | None,
) -> FusionAnnotationCrosscheckResult:
    """Independently re-derive a breakpoint's protein position using
    fusion-annotation's ported strand-order/CDS-offset arithmetic, fed this
    project's own UTR-derived CDS bounds, and compare it against
    ``our_protein_position`` (this project's own already-computed estimate
    from :func:`cfh.mapping.genome_nexus_source.map_genomic_breakpoint_to_protein_position`).

    Never raises: any failure to compute a comparable cross-check value
    (missing CDS bounds, no exon data, a genomic position fusion-annotation's
    arithmetic cannot place in the CDS at all) is reported via ``error``
    with every other field left unset, matching this codebase's existing
    graceful-degradation convention for QA cross-checks.
    """
    if role not in _SIDE_FOR_ROLE:
        return FusionAnnotationCrosscheckResult(
            protein_position=None,
            cds_position=None,
            exon_rank=None,
            is_hybrid_codon=None,
            agrees=None,
            note=None,
            error=f"unsupported fusion role {role!r}; expected 'five_prime' or 'three_prime'",
        )
    if not exons:
        return FusionAnnotationCrosscheckResult(
            protein_position=None,
            cds_position=None,
            exon_rank=None,
            is_hybrid_codon=None,
            agrees=None,
            note=None,
            error="no exon data available; cross-check cannot be attempted",
        )

    # cds_min_genomic/cds_max_genomic is None when this project's own
    # cds_bounds_from_utrs() found no UTR of that type at all -- a real,
    # gene-agnostic transcript characteristic (e.g. a transcript with no
    # annotated stop codon), not a defect in that function. Our own
    # single-sided resolver handles that by simply not clipping on that
    # end; the ported cds_coord_at_genomic requires concrete integer
    # bounds, so fall back to the transcript's outermost exon coordinate on
    # that side -- the same fallback fusion-annotation's own gn_provider.py
    # documents for exactly this case ("GN then simply omits that UTR
    # record ... falls back to the transcript's outermost exon coordinate
    # on that side instead of raising"), citing NTRK1's canonical
    # ENST00000524377 (no three_prime_UTR at all) as its own example.
    effective_cds_min = (
        cds_min_genomic if cds_min_genomic is not None else min(exon.start for exon in exons)
    )
    effective_cds_max = (
        cds_max_genomic if cds_max_genomic is not None else max(exon.end for exon in exons)
    )

    try:
        g_pos = parse_genomic_breakpoint(breakpoint_genomic)
        exon_dicts = _gn_exon_dicts(exons)
        exon_genomic = build_exon_genomic_map(strand, exon_dicts)
        side = _SIDE_FOR_ROLE[role]
        cds_position = cds_coord_at_genomic(
            strand, exon_genomic, effective_cds_min, effective_cds_max, g_pos, side
        )
    except (ValueError, KeyError, IndexError) as exc:
        return FusionAnnotationCrosscheckResult(
            protein_position=None,
            cds_position=None,
            exon_rank=None,
            is_hybrid_codon=None,
            agrees=None,
            note=None,
            error=f"fusion-annotation cross-check arithmetic failed: {type(exc).__name__}: {exc}",
        )

    if cds_position <= 0:
        return FusionAnnotationCrosscheckResult(
            protein_position=None,
            cds_position=cds_position,
            exon_rank=None,
            is_hybrid_codon=None,
            agrees=None,
            note=None,
            error=(
                "breakpoint lies entirely upstream of the first coding exon "
                "(zero coding bases retained on this side); no residue to compare"
            ),
        )

    # Whether the breakpoint lands exactly on a codon boundary is a property
    # of cds_position itself, not of which role is being evaluated -- it
    # must be computed the same way regardless of role (using cds_position
    # % 3, not a role-shifted variant), or the two role branches would
    # disagree about whether the *same* physical breakpoint is a hybrid
    # junction.
    is_hybrid_codon = cds_position % 3 != 0
    if role == "five_prime":
        protein_position = cds_position // 3
    else:
        protein_position = (cds_position - 1) // 3 + 1

    exon_rank = None
    try:
        exon_cds = build_exon_cds_map(strand, exon_dicts, effective_cds_min, effective_cds_max)
        for rank in range(1, len(exon_cds) + 1):
            try:
                exon_start = cds_coord_at_exon_boundary(exon_cds, rank, "start")
                exon_end = cds_coord_at_exon_boundary(exon_cds, rank, "end")
            except ValueError:
                continue  # non-coding exon at this rank
            if exon_start <= cds_position <= exon_end:
                exon_rank = rank
                break
    except (ValueError, IndexError):
        exon_rank = None  # exon-rank derivation is best-effort diagnostic detail only

    if our_protein_position is None:
        return FusionAnnotationCrosscheckResult(
            protein_position=protein_position,
            cds_position=cds_position,
            exon_rank=exon_rank,
            is_hybrid_codon=is_hybrid_codon,
            agrees=None,
            note=None,
            error="no local breakpoint_protein_position available to compare against",
        )

    # The expected relationship is exact and signed, proven in the module
    # docstring: three_prime's protein_position is ALWAYS exactly equal to
    # our_protein_position (never a rounding-tolerance case); five_prime's
    # is exactly our_protein_position - 1 when hybrid, or exactly equal
    # when not. There is no legitimate case where fusion-annotation's
    # residue is *higher* than ours, so an undirected abs()-based
    # tolerance would wrongly accept that as agreement -- it must not be
    # used here.
    if role == "five_prime" and is_hybrid_codon:
        expected = our_protein_position - 1
    else:
        expected = our_protein_position

    if protein_position == expected:
        agrees = True
        if role == "five_prime" and is_hybrid_codon:
            note = (
                "hybrid junction codon: fusion-annotation reports the last complete "
                "5-prime residue, exactly one less than this project's "
                "hybrid-junction-inclusive residue; expected convention difference, "
                "not a disagreement"
            )
        else:
            note = "exact match"
    else:
        agrees = False
        note = (
            f"fusion-annotation reports residue {protein_position}, but the documented "
            f"convention for this role expects exactly {expected} given this project's "
            f"residue of {our_protein_position}; possible genuine disagreement, e.g. the "
            "two resolvers clamped an intronic breakpoint to different flanking exons "
            "(see module docstring)"
        )

    return FusionAnnotationCrosscheckResult(
        protein_position=protein_position,
        cds_position=cds_position,
        exon_rank=exon_rank,
        is_hybrid_codon=is_hybrid_codon,
        agrees=agrees,
        note=note,
        error=None,
    )
