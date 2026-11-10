"""Classify a predicted fusion junction against the domain it rebuilds.

The call comes from sequence, not from the fractional spans in
:mod:`cfh.mapping.domain_reconstitution`. The chimeric domain (the 5' parent's
residues up to the junction, any hybrid residue, then the 3' parent's residues
from the junction) is aligned globally against each parent domain, and the
best-scoring parent decides. Only indels touching the junction count, because
paralogs often differ by indels elsewhere in a domain.

There is no Pfam HMM here: the parents stand in for the family, so the call
says how the junction compares with these two parents only.

Classes, in the order they are decided:

``no_transcript_model``
    No credible chimeric transcript product can be predicted; ``reason`` says why.
``not_split_on_both``
    A product exists, but the family's domain is not cut in both partners.
``out_of_frame``
    The junction shifts the reading frame.
``premature_stop_at_junction``
    In frame, but the junction-spanning codon is a stop.
``tandem_duplication_or_insertion``
    The chimeric domain has extra residues at the junction.
``gap_or_deletion``
    The chimeric domain is missing residues at the junction.
``clean_reconstitution``
    In frame, both halves present, no junction indel against the parent.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

from cfh.mapping.junction_translation import (
    DOMAIN_TRUNCATED,
    FRAME_OUT,
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

RECONSTITUTION_CLASSES = (
    CLEAN_RECONSTITUTION,
    TANDEM_DUPLICATION_OR_INSERTION,
    GAP_OR_DELETION,
    OUT_OF_FRAME,
    NO_TRANSCRIPT_MODEL,
    NOT_SPLIT_ON_BOTH,
    PREMATURE_STOP_AT_JUNCTION,
)

GAP_OPEN = -10
GAP_EXTEND = -1
JUNCTION_WINDOW = 5
"""Alignment columns either side of the junction in which an indel counts.
A tandem copy can be placed on either side of the junction, so the window
must reach past the junction itself."""

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


def align_global(first: str, second: str) -> tuple[str, str, int]:
    """Affine-gap global alignment; returns aligned strings and the score."""
    n, m = len(first), len(second)
    low = -(10**9)
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


def _column_of_residue(aligned: str, residue_index: int) -> int:
    seen = -1
    for column, char in enumerate(aligned):
        if char != "-":
            seen += 1
            if seen == residue_index:
                return column
    raise ValueError(f"residue {residue_index} is not in the alignment")


def _gap_runs(aligned: str) -> list[tuple[int, int]]:
    runs: list[tuple[int, int]] = []
    start: int | None = None
    for column, char in enumerate(aligned):
        if char == "-" and start is None:
            start = column
        elif char != "-" and start is not None:
            runs.append((start, column))
            start = None
    if start is not None:
        runs.append((start, len(aligned)))
    return runs


def _is_tandem(chimera: str, residues: str, junction_index: int) -> bool:
    """True when ``residues`` occurs twice in a row near the junction."""
    size = len(residues)
    if not size:
        return False
    for begin in range(max(0, junction_index - 2 * size), junction_index + 1):
        if chimera[begin : begin + size] == residues == chimera[begin + size : begin + 2 * size]:
            return True
    return False


def _junction_indel(
    aligned_chimera: str, aligned_parent: str, junction_index: int
) -> JunctionIndel | None:
    """The largest non-terminal indel touching the junction, if any."""
    junction_column = _column_of_residue(aligned_chimera, junction_index)
    low, high = junction_column - JUNCTION_WINDOW, junction_column + JUNCTION_WINDOW + 1
    columns = len(aligned_chimera)
    best: JunctionIndel | None = None
    chimera = aligned_chimera.replace("-", "")
    for sequence, other, sign in (
        (aligned_parent, aligned_chimera, 1),
        (aligned_chimera, aligned_parent, -1),
    ):
        for start, end in _gap_runs(sequence):
            if start == 0 or end == columns or end <= low or start >= high:
                continue
            residues = other[start:end].replace("-", "")
            size = len(residues)
            if best is None or size > abs(best.length):
                tandem = sign == 1 and _is_tandem(chimera, residues, junction_index)
                best = JunctionIndel(sign * size, residues, tandem)
    return best


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

    best: tuple[int, str, str, str] | None = None
    for parent_gene, parent_domain in (
        (prediction.five_prime_gene, five_protein[five_span.start_aa - 1 : five_span.end_aa]),
        (prediction.three_prime_gene, three_protein[three_span.start_aa - 1 : three_span.end_aa]),
    ):
        aligned_chimera, aligned_parent, score = align_global(chimera, parent_domain)
        if best is None or score > best[0]:
            best = (score, parent_gene, aligned_chimera, aligned_parent)
    assert best is not None
    _, parent_gene, aligned_chimera, aligned_parent = best
    indel = _junction_indel(aligned_chimera, aligned_parent, junction_index)
    aligned = replace(
        split,
        chimeric_domain=chimera,
        alignment_parent=parent_gene,
        alignment_identity=_identity(aligned_chimera, aligned_parent),
        junction_indel=indel,
    )
    if indel is None:
        return replace(
            aligned,
            classification=CLEAN_RECONSTITUTION,
            reason=f"chimeric {family} domain aligns to {parent_gene} with no junction indel",
        )
    if indel.length > 0:
        kind = "tandem duplication" if indel.is_tandem else "insertion"
        return replace(
            aligned,
            classification=TANDEM_DUPLICATION_OR_INSERTION,
            reason=(
                f"{indel.length}-residue {kind} at the junction ({indel.residues}) against "
                f"{parent_gene}"
            ),
        )
    return replace(
        aligned,
        classification=GAP_OR_DELETION,
        reason=(
            f"{-indel.length}-residue deletion at the junction ({indel.residues}) against "
            f"{parent_gene}"
        ),
    )
