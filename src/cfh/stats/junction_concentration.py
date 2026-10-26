"""Do a gene's fusion junctions pile up at one exon more than random breakage?

The observed statistic asks which ``(role, junction exon)`` holds the most
fusions relative to what random breakage would put there. The null redraws
every breakpoint uniformly across the gene's genomic span, keeps each
event's role (5' or 3' partner), and maps it to a junction exon with the same
:func:`~cfh.mapping.genome_nexus_source.map_genomic_breakpoint_to_protein_position`
rule that labels the observed events. That rule is constant between exon
edges and intron midpoints, so each exon's chance of being the junction is
the length of DNA that maps to it, and a permutation replicate is one
multinomial draw per role.

Panel caveat: the null assumes a breakpoint anywhere in the gene would be
detected. A hybrid-capture panel such as MSK-IMPACT tiles only selected
introns, so breakpoints outside them are missed and the observed junctions
are concentrated partly by assay design. On panel data the p-value therefore
overstates the evidence for biological concentration; cohorts sequenced
without intron selection (whole genome, RNA) do not have this bias.
"""

from __future__ import annotations

import math
from collections import Counter
from collections.abc import Sequence
from typing import Any

import numpy as np
from scipy.stats import binom

from cfh.mapping.genome_nexus_source import (
    ExonRecord,
    _coding_interval,
    map_genomic_breakpoint_to_protein_position,
)

ROLES = ("five_prime", "three_prime")
DEFAULT_MIN_EVENTS = 5
PANEL_CAPTURE_CAVEAT = (
    "The random-breakage null assumes a breakpoint anywhere in the gene would be "
    "detected. Capture panels such as MSK-IMPACT tile only selected introns, so "
    "concentration on panel data is partly by assay design and this p-value "
    "overstates the evidence for biological concentration. Confirm in a cohort "
    "sequenced without intron selection (whole genome or RNA)."
)
_CHUNK_SIZE = 20_000


def junction_probabilities(
    exons: Sequence[ExonRecord],
    role: str,
    *,
    cds_min_genomic: int | None = None,
    cds_max_genomic: int | None = None,
) -> dict[int, float]:
    """Return ``{junction exon rank: probability}`` for one uniformly random
    breakpoint across the transcript's genomic span, given the gene's role.
    """
    if role not in ROLES:
        raise ValueError(f"role must be one of {ROLES}")
    if not exons:
        raise ValueError("no exons to draw breakpoints across")
    strand = exons[0].strand
    span_start = min(exon.start for exon in exons)
    span_end = max(exon.end for exon in exons)
    intervals = sorted(
        interval
        for exon in exons
        if (interval := _coding_interval(exon, cds_min_genomic, cds_max_genomic)) is not None
    )
    if not intervals:
        raise ValueError("no coding exons to map breakpoints to")

    # The mapped exon can only change at a coding-exon edge or, for the
    # nearest-exon fallback, halfway between two coding exons.
    edges = {span_start, span_end + 1}
    for start, end in intervals:
        edges.update((start, end + 1))
    for (_, left_end), (right_start, _) in zip(intervals, intervals[1:]):
        edges.add((left_end + right_start + 1) // 2)
    edges_sorted = sorted(edge for edge in edges if span_start <= edge <= span_end + 1)

    lengths: Counter[int] = Counter()
    for left, right in zip(edges_sorted, edges_sorted[1:]):
        if right <= left:
            continue
        exon_rank = map_genomic_breakpoint_to_protein_position(
            list(exons),
            (left + right - 1) // 2,
            strand,
            cds_min_genomic=cds_min_genomic,
            cds_max_genomic=cds_max_genomic,
            role=role,
        ).exon_rank
        lengths[exon_rank] += right - left
    total = sum(lengths.values())
    return {rank: length / total for rank, length in sorted(lengths.items())}


def junction_concentration_test(
    junctions: Sequence[tuple[str, int]],
    probabilities_by_role: dict[str, dict[int, float]],
    *,
    seed: int = 42,
    n_permutations: int = 10_000,
) -> dict:
    """Test whether any ``(role, exon)`` junction holds far more of the gene's
    fusions than random breakage predicts.

    ``junctions`` holds one ``(role, junction exon rank)`` per observation.
    Each junction's excess is the binomial upper-tail ``-log10 p`` of its
    count among the role's events, given its random-breakage probability. The
    statistic is the largest excess over all junctions, so the null draws,
    which take the same maximum, correct for having looked at every exon.
    Shares are within the role: "83% of fusions with the gene as 3' partner".
    """
    if n_permutations <= 0:
        raise ValueError("n_permutations must be positive")
    if not junctions:
        raise ValueError("no junctions to test")
    counts = Counter(junctions)
    role_sizes = Counter(role for role, _ in junctions)
    excess_by_role = {
        role: _excess_lookup(probabilities_by_role[role], role_sizes[role]) for role in role_sizes
    }
    exon_index = {
        role: {exon: index for index, exon in enumerate(probabilities_by_role[role])}
        for role in role_sizes
    }

    def excess(role: str, exon: int, count: int) -> float:
        index = exon_index[role].get(exon)
        # An exon random breakage cannot reach makes any count infinitely surprising.
        return math.inf if index is None else float(excess_by_role[role][index, count])

    rows: list[dict[str, Any]] = [
        {
            "role": role,
            "junction_exon": exon,
            "count": count,
            "role_event_count": role_sizes[role],
            "observed_share": count / role_sizes[role],
            "random_breakage_share": probabilities_by_role[role].get(exon, 0.0),
            "excess_neg_log10_p": excess(role, exon, count),
        }
        for (role, exon), count in sorted(counts.items())
    ]
    top = max(rows, key=lambda row: (row["excess_neg_log10_p"], row["count"]))
    observed = float(top["excess_neg_log10_p"])

    rng = np.random.default_rng(seed)
    null_max = np.empty(n_permutations)
    for offset in range(0, n_permutations, _CHUNK_SIZE):
        size = min(_CHUNK_SIZE, n_permutations - offset)
        chunk_max = np.zeros(size)
        for role in sorted(role_sizes):
            probabilities = list(probabilities_by_role[role].values())
            draws = rng.multinomial(role_sizes[role], probabilities, size=size)
            lookup = excess_by_role[role]
            values = lookup[np.arange(lookup.shape[0])[None, :], draws]
            np.maximum(chunk_max, values.max(axis=1), out=chunk_max)
        null_max[offset : offset + size] = chunk_max

    # A small tolerance keeps ties from floating-point noise counting as ties.
    p_value = (1 + int((null_max >= observed - 1e-9).sum())) / (n_permutations + 1)
    return {
        "n_events": len(junctions),
        "top_junction": {key: top[key] for key in top},
        "permutation_empirical_p_value": p_value,
        "junction_counts": sorted(rows, key=lambda row: (-row["count"], row["role"])),
    }


def _excess_lookup(probabilities: dict[int, float], n_events: int) -> np.ndarray:
    """``[exon index, count] -> -log10 P(X >= count)`` for X ~ Binomial(n, p)."""
    counts = np.arange(n_events + 1)
    table = np.empty((len(probabilities), n_events + 1))
    for index, probability in enumerate(probabilities.values()):
        tail = binom.sf(counts - 1, n_events, probability)
        with np.errstate(divide="ignore"):
            table[index] = -np.log10(np.clip(tail, 0.0, 1.0))
    return table
