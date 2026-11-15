"""Profile-HMM domain integrity with HMMER (pyhmmer) and per-family Pfam HMMs.

Only the families needed are fetched, one HMM each, from the InterPro API
(``/api/entry/pfam/<acc>/?annotation=hmm``) and cached gzipped. Scores come
from HMMER's own pipeline; nothing here re-implements HMM scoring.

What a passing score means: a domain score at or above the family's Pfam
gathering threshold (GA, domain-level value) supports *family membership*.
It is not evidence of a correct fold or of activity. HMMER also absorbs extra
residues in insert states, so a junction insertion can sit in insert states and
still score well; :class:`DomainHit` therefore counts junction residues aligned
to insert states.
"""

from __future__ import annotations

import gzip
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import pyhmmer
import requests
from pyhmmer.plan7 import HMM, HMMFile, Pipeline

INTERPRO_API = "https://www.ebi.ac.uk/interpro/api"
_ALPHABET = pyhmmer.easel.Alphabet.amino()


class HmmUnavailable(Exception):
    """The family's HMM could not be fetched or has no gathering threshold."""


def fetch_hmm(
    accession: str, cache_dir: Path, *, session: requests.Session | None = None, retries: int = 3
) -> Path:
    """Path to the cached gzipped HMM for Pfam ``accession`` (``PF00178``)."""
    accession = accession.split(".")[0]
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir / f"{accession}.hmm.gz"
    if path.exists():
        return path
    session = session or requests.Session()
    url = f"{INTERPRO_API}/entry/pfam/{accession}/?annotation=hmm"
    last: Exception | None = None
    for attempt in range(retries):
        try:
            response = session.get(url, timeout=60)
            if response.status_code == 404:
                raise HmmUnavailable(f"{accession}: InterPro has no HMM (404)")
            response.raise_for_status()
            payload = response.content
            if not gzip.decompress(payload).startswith(b"HMMER3"):
                raise HmmUnavailable(f"{accession}: response is not a HMMER3 file")
            tmp = path.with_suffix(".tmp")
            tmp.write_bytes(payload)
            tmp.replace(path)
            return path
        except HmmUnavailable:
            raise
        except (requests.RequestException, OSError) as error:
            last = error
            time.sleep(2 * (attempt + 1))
    raise HmmUnavailable(f"{accession}: fetch failed: {last}")


def fetch_pfam_release(session: requests.Session | None = None) -> dict[str, str]:
    """Pfam and InterPro release names/dates served by the InterPro API."""
    session = session or requests.Session()
    response = session.get(f"{INTERPRO_API}/", timeout=60)
    response.raise_for_status()
    databases = response.json()["databases"]
    return {
        f"{name}_{key}": str(databases[name][source])
        for name in ("pfam", "interpro")
        for key, source in (("version", "version"), ("release_date", "releaseDate"))
    }


def load_hmm(path: Path) -> HMM:
    with HMMFile(str(path)) as handle:
        hmm = handle.read()
    if hmm is None:
        raise HmmUnavailable(f"{path}: empty HMM file")
    if not hmm.cutoffs.gathering_available():
        raise HmmUnavailable(f"{path}: HMM has no gathering (GA) threshold")
    return hmm


def gathering_thresholds(hmm: HMM) -> tuple[float, float]:
    """``(sequence GA, domain GA)`` in bits."""
    ga = hmm.cutoffs.gathering
    assert ga is not None
    return float(ga[0]), float(ga[1])


@dataclass(frozen=True)
class DomainHit:
    """One HMMER domain alignment, in 1-based coordinates of the scanned sequence."""

    score: float
    i_evalue: float
    c_evalue: float
    env_from: int
    env_to: int
    ali_from: int
    ali_to: int
    hmm_from: int
    hmm_to: int
    hmm_length: int
    insert_positions: tuple[int, ...]
    """Target positions aligned to HMM insert states."""

    @property
    def match_state_coverage(self) -> float:
        """Fraction of the HMM's match states spanned by the alignment."""
        return (self.hmm_to - self.hmm_from + 1) / self.hmm_length


@dataclass(frozen=True)
class SequenceScan:
    sequence_score: float
    """HMMER full-sequence bit score (all domains, null2-corrected)."""
    domains: tuple[DomainHit, ...]


def _insert_positions(target_from: int, hmm_row: str, target_row: str) -> tuple[int, ...]:
    positions = []
    position = target_from
    for hmm_char, target_char in zip(hmm_row, target_row, strict=True):
        if target_char == "-":
            continue
        if hmm_char == ".":
            positions.append(position)
        position += 1
    return tuple(positions)


def scan(hmm: HMM, sequences: Sequence[tuple[str, str]]) -> dict[str, SequenceScan | None]:
    """Score ``(name, protein)`` pairs against ``hmm``; ``None`` = no reported hit.

    Filters are switched off (``--max``) and every hit is reported, so a weak
    chimeric domain still gets a score; E-values use ``Z = domZ = 1`` (a single
    per-sequence comparison), and gathering comparisons use bit scores only.
    """
    digital = [
        pyhmmer.easel.TextSequence(name=name.encode(), sequence=protein).digitize(_ALPHABET)
        for name, protein in sequences
    ]
    pipeline = Pipeline(
        _ALPHABET,
        Z=1,
        domZ=1,
        E=1e12,
        domE=1e12,
        bias_filter=False,
        F1=1.0,
        F2=1.0,
        F3=1.0,
    )
    hits = pipeline.search_hmm(hmm, pyhmmer.easel.DigitalSequenceBlock(_ALPHABET, digital))
    found: dict[str, SequenceScan | None] = {name: None for name, _ in sequences}
    for hit in hits:
        name = hit.name.decode() if isinstance(hit.name, bytes) else str(hit.name)
        domains = []
        for domain in hit.domains:
            alignment = domain.alignment
            domains.append(
                DomainHit(
                    score=float(domain.score),
                    i_evalue=float(domain.i_evalue),
                    c_evalue=float(domain.c_evalue),
                    env_from=int(domain.env_from),
                    env_to=int(domain.env_to),
                    ali_from=int(alignment.target_from),
                    ali_to=int(alignment.target_to),
                    hmm_from=int(alignment.hmm_from),
                    hmm_to=int(alignment.hmm_to),
                    hmm_length=int(alignment.hmm_length),
                    insert_positions=_insert_positions(
                        int(alignment.target_from),
                        str(alignment.hmm_sequence),
                        str(alignment.target_sequence),
                    ),
                )
            )
        found[name] = SequenceScan(float(hit.score), tuple(domains))
    return found


def best_overlapping(
    scan_result: SequenceScan | None, start: int, end: int
) -> tuple[DomainHit | None, int]:
    """The domain hit overlapping ``[start, end]`` most (ties: higher score), and hit count."""
    if scan_result is None:
        return None, 0

    def overlap(hit: DomainHit) -> int:
        return max(0, min(end, hit.ali_to) - max(start, hit.ali_from) + 1)

    candidates = [hit for hit in scan_result.domains if overlap(hit) > 0]
    if not candidates:
        return None, len(scan_result.domains)
    return max(candidates, key=lambda hit: (overlap(hit), hit.score)), len(scan_result.domains)


def residues_in_inserts(hit: DomainHit | None, positions: Iterable[int]) -> int | None:
    """How many of ``positions`` the alignment puts in insert states (None: no hit)."""
    if hit is None:
        return None
    inserts = set(hit.insert_positions)
    return sum(1 for position in positions if position in inserts)


def overlapping_hits(scan_result: SequenceScan | None, start: int, end: int) -> list[DomainHit]:
    """Every domain hit whose alignment overlaps ``[start, end]``."""
    if scan_result is None:
        return []
    return [hit for hit in scan_result.domains if min(end, hit.ali_to) >= max(start, hit.ali_from)]


def union_match_coverage(hits: Sequence[DomainHit]) -> float | None:
    """Fraction of HMM match states covered by the union of ``hits`` (None: no hits)."""
    if not hits:
        return None
    covered: set[int] = set()
    for hit in hits:
        covered.update(range(hit.hmm_from, hit.hmm_to + 1))
    return len(covered) / hits[0].hmm_length
