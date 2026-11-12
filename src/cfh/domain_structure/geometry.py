"""Superposition, RMSD and TM-score on CA coordinates with a given residue mapping.

The structure comparisons in this package always use an explicit residue
correspondence (chimera residue -> parent residue from the parent alignment of
:mod:`cfh.domain_function.metrics`), so the scores say how each *mapped*
residue sits, not how well two shapes can be matched in any register.

* :func:`kabsch` -- least-squares rotation/translation (Kabsch 1976).
* :func:`tm_score_fixed` -- TM-score (Zhang & Skolnick 2004, Proteins 57:702)
  for a fixed correspondence: the rotation maximising the TM sum is searched
  as in the TM-score program (seed fragments of length L, L/2, L/4, ... >= 4,
  each refined by re-superposing the pairs closer than a growing cutoff).
  Normalised by the length passed in (here the chimera's domain length).
* :func:`tm_align` -- sequence-independent TM-align via ``tmtools``, reported
  alongside as a secondary, register-free check.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

Coords = NDArray[np.float64]


def kabsch(mobile: Coords, target: Coords) -> tuple[Coords, Coords]:
    """Rotation ``R`` and translation ``t`` minimising ``|mobile @ R.T + t - target|``."""
    if mobile.shape != target.shape or mobile.shape[0] < 3:
        raise ValueError("kabsch needs two equal-length sets of at least 3 points")
    mc, tc = mobile.mean(axis=0), target.mean(axis=0)
    h = (mobile - mc).T @ (target - tc)
    u, _, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    correction = np.diag([1.0, 1.0, d])
    rotation = vt.T @ correction @ u.T
    return rotation, tc - mc @ rotation.T


def apply(coords: Coords, rotation: Coords, translation: Coords) -> Coords:
    return coords @ rotation.T + translation


def rmsd(a: Coords, b: Coords) -> float:
    return float(np.sqrt(((a - b) ** 2).sum(axis=1).mean()))


def superposed_rmsd(mobile: Coords, target: Coords) -> float:
    """RMSD after optimal superposition of all pairs."""
    rotation, translation = kabsch(mobile, target)
    return rmsd(apply(mobile, rotation, translation), target)


def tm_d0(length: int) -> float:
    """TM-score distance scale; 0.5 floor for short chains as in the TM-score program."""
    if length <= 21:
        return 0.5
    return max(0.5, 1.24 * (length - 15) ** (1.0 / 3.0) - 1.8)


def _tm_sum(distances: NDArray[np.float64], d0: float) -> float:
    return float((1.0 / (1.0 + (distances / d0) ** 2)).sum())


@dataclass(frozen=True)
class TmResult:
    tm_score: float
    """TM-score normalised by ``norm_length``."""
    rmsd: float
    """RMSD of all mapped pairs after their own least-squares superposition."""
    n_aligned: int
    norm_length: int


def tm_score_fixed(mobile: Coords, target: Coords, norm_length: int) -> TmResult:
    """TM-score of mapped pairs ``mobile[i] <-> target[i]`` (see module docstring)."""
    n = mobile.shape[0]
    if n < 3 or norm_length < 1:
        raise ValueError("need at least 3 mapped pairs")
    d0 = tm_d0(norm_length)
    d_search = min(max(4.5, d0), 8.0)
    best = 0.0
    fragment = n
    tried: set[tuple[int, int]] = set()
    while fragment >= min(4, n):
        step = max(1, fragment // 2)
        for start in range(0, n - fragment + 1, step):
            if (start, fragment) in tried:
                continue
            tried.add((start, fragment))
            selected = np.arange(start, start + fragment)
            for _ in range(20):
                if selected.size < 3:
                    break
                rotation, translation = kabsch(mobile[selected], target[selected])
                distances = np.linalg.norm(apply(mobile, rotation, translation) - target, axis=1)
                best = max(best, _tm_sum(distances, d0))
                cutoff = d_search
                new = np.flatnonzero(distances < cutoff)
                while new.size < 3 and cutoff < 30:
                    cutoff += 0.5
                    new = np.flatnonzero(distances < cutoff)
                if np.array_equal(new, selected):
                    break
                selected = new
        if fragment == min(4, n):
            break
        fragment = max(min(4, n), fragment // 2)
    return TmResult(
        tm_score=round(best / norm_length, 4),
        rmsd=round(superposed_rmsd(mobile, target), 3),
        n_aligned=n,
        norm_length=norm_length,
    )


@dataclass(frozen=True)
class TmAlignResult:
    tm_norm_first: float
    tm_norm_second: float
    rmsd: float
    n_aligned: int


def tm_align(first: Coords, second: Coords, seq_first: str, seq_second: str) -> TmAlignResult:
    """Sequence-independent TM-align (``tmtools``); secondary to the mapped scores."""
    from tmtools import tm_align as _tm_align

    result = _tm_align(
        np.ascontiguousarray(first), np.ascontiguousarray(second), seq_first, seq_second
    )
    aligned = sum(
        1 for a, b in zip(result.seqxA, result.seqyA, strict=True) if a != "-" and b != "-"
    )
    return TmAlignResult(
        tm_norm_first=round(float(result.tm_norm_chain1), 4),
        tm_norm_second=round(float(result.tm_norm_chain2), 4),
        rmsd=round(float(result.rmsd), 3),
        n_aligned=aligned,
    )


def finite_or_none(value: float | None) -> float | None:
    return None if value is None or not math.isfinite(value) else value
