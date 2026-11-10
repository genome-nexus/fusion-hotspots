"""Classify a predicted fusion junction against the domain it rebuilds.

The call comes from sequence, not from the fractional spans in
:mod:`cfh.mapping.domain_reconstitution`. The chimeric domain is the 5'
parent's residues up to the junction, any hybrid residue, then the 3' parent's
residues from the junction.

A junction duplication or deletion is read from the junction coordinates, not
from where an aligner places a gap in the chimera. The two parent domains are
aligned with each other; for each co-optimal alignment, ``k`` counts the 3'
residues at or before the 5' end's column, and the in-register parent is the
5' domain through the junction followed by the 3' domain from residue ``k``.
The chimera's length minus the in-register parent's is the indel: a restart
before ``k`` repeats residues, a restart after it skips them. Every co-optimal
alignment is considered (forward and backward dynamic programming), so a
paralog gap that can sit equally well on either side of the junction, as in a
repeat, cannot flip the call through a tie-break: if the alignments disagree
the call is ``ambiguous_junction_register``. ``clean_reconstitution`` needs no
indel under any of them and a chimera the length of the parent it is reported
against. Paralog indels away from the junction do not count. The chimera is
also aligned (BLOSUM62, affine gaps) against each parent domain; the parent
reported is the same-length one if any, else the best-scoring one.

There is no Pfam HMM here: the parents stand in for the family, so the call
says how the junction compares with these two parents only.

Classes, in the order they are decided:

``splice_site_unresolved``
    An intronic break falls inside the splice donor or acceptor the product
    would use, so whether any product is made is unresolved; ``reason`` says
    which site and how many intronic bases the break keeps.
``no_transcript_model``
    No credible chimeric transcript product can be predicted; ``reason`` says why.
``not_split_on_both``
    A product exists, but the family's domain is not cut in both partners.
``out_of_frame``
    The junction shifts the reading frame.
``premature_stop_at_junction``
    In frame, but the junction-spanning codon is a stop.
``ambiguous_junction_register``
    In frame and cut in both partners, but equally good alignments of the two
    parent domains (a paralog gap that can sit on either side of the junction,
    e.g. inside a repeat) disagree on whether residues are repeated or
    skipped, or the chimera is in register but matches neither parent's
    length. Never clean; ``reason`` lists the competing indel sizes.
``tandem_duplication_or_insertion``
    The chimeric domain has extra residues at the junction.
``gap_or_deletion``
    The chimeric domain is missing residues at the junction.
``clean_reconstitution``
    In frame, both halves present, no residue repeated or skipped at the
    junction under every co-optimal parent alignment, and the same length as
    the parent it is reported against.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

from cfh.mapping.junction_translation import (
    DOMAIN_TRUNCATED,
    FRAME_OUT,
    MECHANISM_SPLICE_SITE_UNRESOLVED,
    DomainSpan,
    JunctionPrediction,
    TranscriptModel,
    translate,
)

CLEAN_RECONSTITUTION = "clean_reconstitution"
TANDEM_DUPLICATION_OR_INSERTION = "tandem_duplication_or_insertion"
GAP_OR_DELETION = "gap_or_deletion"
OUT_OF_FRAME = "out_of_frame"
NO_TRANSCRIPT_MODEL = "no_transcript_model"
NOT_SPLIT_ON_BOTH = "not_split_on_both"
PREMATURE_STOP_AT_JUNCTION = "premature_stop_at_junction"
SPLICE_SITE_UNRESOLVED = "splice_site_unresolved"
AMBIGUOUS_JUNCTION_REGISTER = "ambiguous_junction_register"

RECONSTITUTION_CLASSES = (
    CLEAN_RECONSTITUTION,
    TANDEM_DUPLICATION_OR_INSERTION,
    GAP_OR_DELETION,
    OUT_OF_FRAME,
    NO_TRANSCRIPT_MODEL,
    NOT_SPLIT_ON_BOTH,
    PREMATURE_STOP_AT_JUNCTION,
    SPLICE_SITE_UNRESOLVED,
    AMBIGUOUS_JUNCTION_REGISTER,
)

GAP_OPEN = -10
GAP_EXTEND = -1
_AMINO_ACIDS = "ARNDCQEGHILKMFPSTWYV"
_BLOSUM62_ROWS = """
 4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0
-1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3
-2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3
-2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3
 0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
-1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2
-1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2
 0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3
-2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3
-1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3
-1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1
-1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2
-1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1
-2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1
-1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2
 1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2
 0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0
-3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3
-2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1
 0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4
"""
_BLOSUM62 = {
    (a, b): int(value)
    for a, row in zip(_AMINO_ACIDS, _BLOSUM62_ROWS.split("\n")[1:-1])
    for b, value in zip(_AMINO_ACIDS, row.split())
}


_LOW = -(10**9)
_Matrix = list[list[int]]


def _forward(first: str, second: str) -> tuple[_Matrix, _Matrix, _Matrix]:
    """Affine-gap forward matrices ``(match, gap_second, gap_first)``.

    ``gap_second[i][j]`` ends with ``first[i-1]`` against a gap,
    ``gap_first[i][j]`` with ``second[j-1]`` against a gap.
    """
    n, m = len(first), len(second)
    low = _LOW
    match = [[low] * (m + 1) for _ in range(n + 1)]
    gap_second = [[low] * (m + 1) for _ in range(n + 1)]  # first[i-1] against a gap
    gap_first = [[low] * (m + 1) for _ in range(n + 1)]  # second[j-1] against a gap
    match[0][0] = 0
    for i in range(1, n + 1):
        gap_second[i][0] = GAP_OPEN + GAP_EXTEND * (i - 1)
    for j in range(1, m + 1):
        gap_first[0][j] = GAP_OPEN + GAP_EXTEND * (j - 1)
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            pair = _BLOSUM62.get((first[i - 1], second[j - 1]), -4)
            match[i][j] = pair + max(
                match[i - 1][j - 1], gap_second[i - 1][j - 1], gap_first[i - 1][j - 1]
            )
            gap_second[i][j] = max(
                match[i - 1][j] + GAP_OPEN,
                gap_second[i - 1][j] + GAP_EXTEND,
                gap_first[i - 1][j] + GAP_OPEN,
            )
            gap_first[i][j] = max(
                match[i][j - 1] + GAP_OPEN,
                gap_first[i][j - 1] + GAP_EXTEND,
                gap_second[i][j - 1] + GAP_OPEN,
            )
    return match, gap_second, gap_first


def _backward(first: str, second: str) -> tuple[_Matrix, _Matrix, _Matrix]:
    """Best score still to come from each cell and state, to the end at ``(n, m)``."""
    n, m = len(first), len(second)
    match = [[_LOW] * (m + 1) for _ in range(n + 1)]
    gap_second = [[_LOW] * (m + 1) for _ in range(n + 1)]
    gap_first = [[_LOW] * (m + 1) for _ in range(n + 1)]
    for i in range(n, -1, -1):
        for j in range(m, -1, -1):
            if i == n and j == m:
                match[i][j] = gap_second[i][j] = gap_first[i][j] = 0
                continue
            diagonal = (
                _BLOSUM62.get((first[i], second[j]), -4) + match[i + 1][j + 1]
                if i < n and j < m
                else _LOW
            )
            down = gap_second[i + 1][j] if i < n else _LOW
            right = gap_first[i][j + 1] if j < m else _LOW
            match[i][j] = max(diagonal, down + GAP_OPEN, right + GAP_OPEN)
            gap_second[i][j] = max(diagonal, down + GAP_EXTEND, right + GAP_OPEN)
            gap_first[i][j] = max(diagonal, down + GAP_OPEN, right + GAP_EXTEND)
    return match, gap_second, gap_first


def co_optimal_registers(first: str, second: str, consumed: int) -> list[int]:
    """Every ``j`` such that some co-optimal global alignment of ``first`` and
    ``second`` places ``first[consumed - 1]`` in a column with exactly ``j``
    residues of ``second`` at or before it.

    Equivalent gap placements (for instance anywhere along a poly-alanine
    run) are all optimal, so a quantity read from one traceback can depend on
    an arbitrary tie-break; this returns the whole set instead.
    """
    forward = _forward(first, second)
    backward = _backward(first, second)
    n, m = len(first), len(second)
    best = max(forward[0][n][m], forward[1][n][m], forward[2][n][m])
    found = []
    for j in range(m + 1):
        # The step that consumes first[consumed - 1] ends in the match or the
        # gap_second state at row ``consumed``.
        if any(
            forward[state][consumed][j] > _LOW // 2
            and forward[state][consumed][j] + backward[state][consumed][j] == best
            for state in (0, 1)
        ):
            found.append(j)
    return found


def align_global(first: str, second: str) -> tuple[str, str, int]:
    """Affine-gap global alignment; returns aligned strings and the score."""
    n, m = len(first), len(second)
    match, gap_second, gap_first = _forward(first, second)

    i, j = n, m
    aligned_first: list[str] = []
    aligned_second: list[str] = []
    best = max(match[n][m], gap_second[n][m], gap_first[n][m])
    state = "match" if best == match[n][m] else "second" if best == gap_second[n][m] else "first"
    while i > 0 or j > 0:
        if state == "match":
            pair = _BLOSUM62.get((first[i - 1], second[j - 1]), -4)
            current = match[i][j] - pair
            aligned_first.append(first[i - 1])
            aligned_second.append(second[j - 1])
            i, j = i - 1, j - 1
            if current == match[i][j]:
                state = "match"
            elif current == gap_second[i][j]:
                state = "second"
            else:
                state = "first"
        elif state == "second":
            current = gap_second[i][j]
            aligned_first.append(first[i - 1])
            aligned_second.append("-")
            i -= 1
            if current == gap_second[i][j] + GAP_EXTEND:
                state = "second"
            elif current == match[i][j] + GAP_OPEN:
                state = "match"
            else:
                state = "first"
        else:
            current = gap_first[i][j]
            aligned_first.append("-")
            aligned_second.append(second[j - 1])
            j -= 1
            if current == gap_first[i][j] + GAP_EXTEND:
                state = "first"
            elif current == match[i][j] + GAP_OPEN:
                state = "match"
            else:
                state = "second"
    return "".join(reversed(aligned_first)), "".join(reversed(aligned_second)), best


@dataclass(frozen=True)
class JunctionIndel:
    length: int
    """Residues: positive for an insertion in the chimera, negative for a deletion."""
    residues: str
    is_tandem: bool


@dataclass(frozen=True)
class ReconstitutionCall:
    classification: str
    family: str
    reason: str
    five_prime_domain: DomainSpan | None = None
    three_prime_domain: DomainSpan | None = None
    chimeric_domain: str | None = None
    alignment_parent: str | None = None
    alignment_identity: float | None = None
    junction_indel: JunctionIndel | None = None


def _is_tandem(chimera: str, residues: str, junction_index: int) -> bool:
    """True when ``residues`` occurs twice in a row near the junction."""
    size = len(residues)
    if not size:
        return False
    for begin in range(max(0, junction_index - 2 * size), junction_index + 1):
        if chimera[begin : begin + size] == residues == chimera[begin + size : begin + 2 * size]:
            return True
    return False


@dataclass(frozen=True)
class _Register:
    """The junction read against one co-optimal register of the two parents."""

    indel: JunctionIndel | None
    reference: str
    """The in-register parent: the 5' domain through the junction, then the 3'
    domain after it."""


def _junction_registers(
    five_domain: str,
    three_domain: str,
    five_end: int,
    three_start: int,
    chimera: str,
    junction_index: int,
) -> list[_Register]:
    """The junction indel under every co-optimal register of the two parents.

    ``five_end`` is the last 5' domain residue the chimera holds (zero-based,
    counting a hybrid residue as the 5' partner's), ``three_start`` the first
    3' domain residue it holds. For each co-optimal alignment of the two parent
    domains, ``k`` is the number of 3' residues at or before the column of the
    5' end. The in-register parent is then ``five_domain[:five_end + 1] +
    three_domain[k:]``: a restart before ``k`` repeats ``three_domain
    [three_start:k]``, a restart after it skips ``three_domain[k:three_start]``.

    The chimera's own length minus the in-register parent's length is the
    indel, so where an aligner puts a gap inside the chimera never matters.
    Every distinct ``k`` is returned, so a caller can tell when equivalent
    placements of a paralog gap disagree.
    """
    registers = []
    for k in co_optimal_registers(five_domain, three_domain, five_end + 1):
        reference = five_domain[: five_end + 1] + three_domain[k:]
        length = len(chimera) - len(reference)
        if length > 0:
            residues = three_domain[three_start:k]
            indel: JunctionIndel | None = JunctionIndel(
                length, residues, _is_tandem(chimera, residues, junction_index)
            )
        elif length < 0:
            indel = JunctionIndel(length, three_domain[k:three_start], False)
        else:
            indel = None
        assert length == k - three_start  # by construction
        registers.append(_Register(indel, reference))
    return registers


def _identity(first: str, second: str) -> float:
    aligned = sum(a != "-" and b != "-" for a, b in zip(first, second))
    same = sum(a == b and a != "-" for a, b in zip(first, second))
    return round(same / aligned, 3) if aligned else 0.0


def _truncated_copy(prediction: JunctionPrediction, partner: str, family: str) -> DomainSpan | None:
    for span in prediction.domain_spans:
        if span.partner == partner and span.family == family and span.status == DOMAIN_TRUNCATED:
            return span
    return None


def _describe_not_split(prediction: JunctionPrediction, family: str) -> str:
    parts = []
    for partner, gene in (
        ("five_prime", prediction.five_prime_gene),
        ("three_prime", prediction.three_prime_gene),
    ):
        states = [
            span.status
            for span in prediction.domain_spans
            if span.partner == partner and span.family == family
        ]
        parts.append(f"{gene}: {', '.join(states) or 'no copy'}")
    return f"{family} not cut in both partners ({'; '.join(parts)})"


def classify_family(
    prediction: JunctionPrediction,
    five_prime_transcript: TranscriptModel | None,
    three_prime_transcript: TranscriptModel | None,
    family: str,
) -> ReconstitutionCall:
    """Classify one Pfam family shared by the two partners."""
    if prediction.mechanism == MECHANISM_SPLICE_SITE_UNRESOLVED:
        return ReconstitutionCall(
            SPLICE_SITE_UNRESOLVED,
            family,
            prediction.no_product_reason or "an intronic break is inside a splice site",
        )
    if not prediction.has_product:
        return ReconstitutionCall(
            NO_TRANSCRIPT_MODEL, family, prediction.no_product_reason or "no product predicted"
        )
    assert five_prime_transcript is not None and three_prime_transcript is not None

    five_span = _truncated_copy(prediction, "five_prime", family)
    three_span = _truncated_copy(prediction, "three_prime", family)
    if five_span is None or three_span is None:
        return ReconstitutionCall(
            NOT_SPLIT_ON_BOTH,
            family,
            _describe_not_split(prediction, family),
            five_prime_domain=five_span,
            three_prime_domain=three_span,
        )
    split = ReconstitutionCall(
        NOT_SPLIT_ON_BOTH,  # placeholder; every return below sets the class
        family,
        "",
        five_prime_domain=five_span,
        three_prime_domain=three_span,
    )
    if prediction.frame_status == FRAME_OUT:
        return replace(
            split, classification=OUT_OF_FRAME, reason="the junction shifts the reading frame"
        )
    if prediction.premature_stop_at_junction:
        return replace(
            split,
            classification=PREMATURE_STOP_AT_JUNCTION,
            reason=f"junction codon {prediction.hybrid_codon} encodes a stop",
        )

    assert prediction.five_prime_last_residue is not None
    assert prediction.three_prime_first_residue is not None
    five_protein = five_prime_transcript.protein
    three_protein = three_prime_transcript.protein
    hybrid_residue = translate(prediction.hybrid_codon) if prediction.hybrid_codon else ""
    five_piece = (
        five_protein[five_span.start_aa - 1 : prediction.five_prime_last_residue] + hybrid_residue
    )
    three_piece = three_protein[prediction.three_prime_first_residue - 1 : three_span.end_aa]
    chimera = five_piece + three_piece
    junction_index = len(five_piece) - 1

    five_domain = five_protein[five_span.start_aa - 1 : five_span.end_aa]
    three_domain = three_protein[three_span.start_aa - 1 : three_span.end_aa]
    registers = _junction_registers(
        five_domain,
        three_domain,
        five_end=len(five_piece) - 1,
        three_start=prediction.three_prime_first_residue - three_span.start_aa,
        chimera=chimera,
        junction_index=junction_index,
    )

    scored: list[tuple[int, str, str, str, str]] = []
    for parent_gene, parent_domain in (
        (prediction.five_prime_gene, five_domain),
        (prediction.three_prime_gene, three_domain),
    ):
        aligned_chimera, aligned_parent, score = align_global(chimera, parent_domain)
        scored.append((score, parent_gene, parent_domain, aligned_chimera, aligned_parent))
    # Same length as the chimera first, then the better alignment score.
    same_length = [entry for entry in scored if len(entry[2]) == len(chimera)]
    _, parent_gene, _, aligned_chimera, aligned_parent = max(
        same_length or scored, key=lambda entry: entry[0]
    )
    lengths = sorted({register.indel.length if register.indel else 0 for register in registers})
    # The largest indel any co-optimal register implies, for the columns.
    indel = max(
        (register.indel for register in registers),
        key=lambda found: abs(found.length) if found else 0,
    )
    aligned = replace(
        split,
        chimeric_domain=chimera,
        alignment_parent=parent_gene,
        alignment_identity=_identity(aligned_chimera, aligned_parent),
        junction_indel=indel,
    )
    if len(lengths) > 1:
        return replace(
            aligned,
            classification=AMBIGUOUS_JUNCTION_REGISTER,
            reason=(
                "co-optimal alignments of the parent domains disagree on the junction "
                f"register: {', '.join(f'{length:+d}' for length in lengths)} residues "
                f"(chimera {len(chimera)}; parents {len(five_domain)} and {len(three_domain)})"
            ),
        )
    reference = registers[0].reference
    if indel is None:
        if not same_length:
            return replace(
                aligned,
                classification=AMBIGUOUS_JUNCTION_REGISTER,
                reason=(
                    f"in register at the junction, but the {len(chimera)}-residue chimera "
                    f"matches neither parent length ({len(five_domain)}, {len(three_domain)})"
                ),
            )
        return replace(
            aligned,
            classification=CLEAN_RECONSTITUTION,
            reason=(
                f"chimeric {family} domain is in register at the junction under every "
                f"co-optimal parent alignment and has the length of {parent_gene}'s "
                f"({len(chimera)} residues)"
            ),
        )
    if indel.length > 0:
        kind = "tandem duplication" if indel.is_tandem else "insertion"
        return replace(
            aligned,
            classification=TANDEM_DUPLICATION_OR_INSERTION,
            reason=(
                f"{indel.length}-residue {kind} at the junction ({indel.residues}); "
                f"chimera {len(chimera)} vs in-register parent {len(reference)} residues"
            ),
        )
    return replace(
        aligned,
        classification=GAP_OR_DELETION,
        reason=(
            f"{-indel.length}-residue deletion at the junction ({indel.residues}); "
            f"chimera {len(chimera)} vs in-register parent {len(reference)} residues"
        ),
    )
