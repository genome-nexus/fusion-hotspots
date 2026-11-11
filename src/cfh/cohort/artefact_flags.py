"""Artefact flags for domain-screen candidate rows.

Three flags, each a value plus a short reason string. A value is one of
``flagged``, ``pass``, ``not_assessed`` (a check that should run but could
not, with the reason) or ``not_applicable``. A check never passes silently:
missing data is ``not_assessed``.

``paralog_sequence_flag``
    The fused CDS ``JUNCTION_CONTEXT_NT`` bases either side of the junction
    (``junction_context_nt``) is aligned locally to each single sequence that
    could have produced it without a fusion: the CDS of the 5' and 3' parent
    transcripts and of every Ensembl Compara paralogue of either partner.
    Flagged when one sequence covers at least ``PARALOG_MIN_COVERAGE`` of the
    window at ``PARALOG_MIN_IDENTITY`` identity, i.e. a mis-mapping or template
    switch on one transcript could explain the "junction".

``genomic_similarity_flag``
    Reference sequence ``+/-150`` and ``+/-500`` bp around the two genomic
    breakpoints are aligned to each other (both orientations). Flagged at
    ``GENOMIC_MIN_IDENTITY`` identity over ``GENOMIC_MIN_LENGTH`` aligned
    columns. Also flagged when either breakpoint lies in a UCSC segmental
    duplication or RepeatMasker element, when the sample carries at least
    ``SAMPLE_BURDEN_MIN`` distinct shared-family calls, or when the call has
    no split (junction) reads or an imprecise breakpoint. Each part has its own
    column; the flag is ``flagged`` when any part is, ``not_assessed`` when
    none is flagged but a part could not be assessed.

``normal_recurrence_flag``
    The exact junction (gene pair and both breakpoint coordinates, on
    GRCh38) occurs in a GTEx or TCGA normal sample of the TCGA+GTEx table.
    Cell-line-labelled GTEx samples are counted separately and never flag.

Alignments are local (Smith-Waterman), linear gap, match +1, mismatch -1,
gap -2. Identity is matches over aligned columns.
"""

from __future__ import annotations

import csv
import gzip
import json
import struct
import threading
import time
from bisect import bisect_left
from collections import Counter, defaultdict
from collections.abc import Iterable, Iterator, Mapping, Sequence
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Protocol

import numpy as np
import requests

from cfh.mapping.genome_nexus_source import TRANSIENT_REQUEST_ERRORS
from cfh.mapping.junction_translation import normalize_genome_build
from cfh.mapping.transcript_models import check_cache_provenance

FLAGGED = "flagged"
PASS = "pass"
NOT_ASSESSED = "not_assessed"
NOT_APPLICABLE = "not_applicable"

PARALOG_MIN_IDENTITY = 0.90
PARALOG_MIN_COVERAGE = 0.90
GENOMIC_WINDOWS = (150, 500)
GENOMIC_MIN_IDENTITY = 0.80
GENOMIC_MIN_LENGTH = 100
SAMPLE_BURDEN_MIN = 2

MATCH, MISMATCH, GAP = 1, -1, 2

FLAG_COLUMNS = (
    "paralog_sequence_flag",
    "paralog_sequence_reason",
    "paralog_best_target",
    "paralog_best_identity",
    "paralog_best_query_coverage",
    "genomic_similarity_flag",
    "genomic_similarity_reason",
    "genomic_identity_150",
    "genomic_aligned_length_150",
    "genomic_identity_500",
    "genomic_aligned_length_500",
    "breakpoint_repeat_flag",
    "breakpoint_repeat_reason",
    "sample_shared_family_calls",
    "sample_burden_flag",
    "read_evidence_flag",
    "read_evidence_reason",
    "normal_recurrence_flag",
    "normal_recurrence_reason",
    "normal_exact_junction_samples",
    "cell_line_exact_junction_samples",
    "normal_gene_pair_samples",
)

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


@dataclass(frozen=True)
class LocalAlignment:
    score: int
    matches: int
    columns: int
    query_start: int
    query_end: int
    """Half-open query interval covered by the alignment."""

    @property
    def identity(self) -> float:
        return self.matches / self.columns if self.columns else 0.0

    @property
    def query_span(self) -> int:
        return self.query_end - self.query_start


_EMPTY = LocalAlignment(0, 0, 0, 0, 0)


def _codes(sequence: str) -> np.ndarray:
    return np.frombuffer(sequence.upper().encode("ascii"), dtype=np.uint8)


def local_align(query: str, target: str) -> LocalAlignment:
    """Best local alignment of ``query`` and ``target`` (see module docstring).

    Each row of the dynamic programme is vectorised: with a linear gap the
    in-row gap recurrence is a running maximum of ``T[k] + GAP*k``. ``N`` never
    matches.
    """
    if not query or not target:
        return _EMPTY
    q, t = _codes(query), _codes(target)
    n, m = len(q), len(t)
    score = np.zeros((n + 1, m + 1), dtype=np.int32)
    index = np.arange(m + 1, dtype=np.int32) * GAP
    n_code = ord("N")
    for i in range(1, n + 1):
        same = (t == q[i - 1]) & (t != n_code)
        substitution = np.where(same, MATCH, MISMATCH).astype(np.int32)
        previous = score[i - 1]
        candidate = np.zeros(m + 1, dtype=np.int32)
        candidate[1:] = np.maximum(previous[:-1] + substitution, previous[1:] - GAP)
        np.maximum(candidate, 0, out=candidate)
        score[i] = np.maximum.accumulate(candidate + index) - index
    best = int(score.max())
    if best <= 0:
        return _EMPTY
    i, j = (int(v) for v in np.unravel_index(int(score.argmax()), score.shape))
    end_i = i
    matches = columns = 0
    while i > 0 and j > 0 and score[i, j] > 0:
        same = q[i - 1] == t[j - 1] and q[i - 1] != n_code
        if score[i, j] == score[i - 1, j - 1] + (MATCH if same else MISMATCH):
            matches += int(same)
            i, j = i - 1, j - 1
        elif score[i, j] == score[i - 1, j] - GAP:
            i -= 1
        else:
            j -= 1
        columns += 1
    return LocalAlignment(best, matches, columns, i, end_i)


SPACER = "N" * 128
"""Separates concatenated targets. N never matches and a local alignment
scores at most the query length (100), so no alignment crosses a spacer."""

CHUNK_BASES = 400_000


def explaining_min_score(query_length: int) -> int:
    """Lowest score any alignment that explains a ``query_length`` window can have.

    Explaining needs a query span of at least ``ceil(0.9 L)`` (so at least that
    many columns, ``c``) at >= 90% identity; with match +1, mismatch -1 and gap -2
    the worst case is all non-matching columns gaps: ``0.9c - 2 * 0.1c = 0.7c``.
    Rounded down, so the cutoff never exceeds the true bound. Targets scoring
    below it cannot explain the window and need not be traced back.
    """
    span = -(-PARALOG_MIN_COVERAGE_TENTHS * query_length // 10)
    return 7 * span // 10


PARALOG_MIN_COVERAGE_TENTHS = 9
""":data:`PARALOG_MIN_COVERAGE` in tenths, for integer arithmetic."""

EXPLAINING_MIN_SCORE = explaining_min_score(100)
"""The cutoff for a full 100 nt window (63)."""


def _column_maxima(query: np.ndarray, target: np.ndarray) -> np.ndarray:
    """Per target column, the best local score of an alignment ending there."""
    m = len(target)
    index = np.arange(m + 1, dtype=np.int32) * GAP
    previous = np.zeros(m + 1, dtype=np.int32)
    best = np.zeros(m + 1, dtype=np.int32)
    n_code = ord("N")
    valid = target != n_code
    for base in query:
        substitution = np.where((target == base) & valid, MATCH, MISMATCH).astype(np.int32)
        candidate = np.zeros(m + 1, dtype=np.int32)
        candidate[1:] = np.maximum(previous[:-1] + substitution, previous[1:] - GAP)
        np.maximum(candidate, 0, out=candidate)
        previous = np.maximum.accumulate(candidate + index) - index
        np.maximum(best, previous, out=best)
    return best[1:]


def target_scores(query: str, targets: Mapping[str, str]) -> dict[str, int]:
    """Best local alignment score of ``query`` against each target."""
    names = [name for name in sorted(targets) if targets[name]]
    q = _codes(query)
    scores: dict[str, int] = {}
    chunk: list[str] = []
    length = 0

    def flush() -> None:
        if not chunk:
            return
        joined = SPACER.join(targets[name] for name in chunk)
        maxima = _column_maxima(q, _codes(joined))
        offset = 0
        for name in chunk:
            size = len(targets[name])
            scores[name] = int(maxima[offset : offset + size].max())
            offset += size + len(SPACER)
        chunk.clear()

    for name in names:
        if length + len(targets[name]) > CHUNK_BASES:
            flush()
            length = 0
        chunk.append(name)
        length += len(targets[name]) + len(SPACER)
    flush()
    return scores


@dataclass(frozen=True)
class SequenceMatch:
    target: str
    alignment: LocalAlignment
    query_length: int

    @property
    def coverage(self) -> float:
        return self.alignment.query_span / self.query_length if self.query_length else 0.0

    @property
    def explains(self) -> bool:
        return (
            self.coverage >= PARALOG_MIN_COVERAGE
            and self.alignment.identity >= PARALOG_MIN_IDENTITY
        )


def best_single_sequence_match(query: str, targets: Mapping[str, str]) -> SequenceMatch | None:
    """The target that best explains all of ``query`` (explaining first, then score).

    All targets are scored at once; the best-scoring local alignment of each
    target scoring at least :func:`explaining_min_score` of the query length, and of the overall
    best target, is traced back. ``None`` when no target aligns at all.
    Memoised on the query and the whole target set: many rows share a junction.
    """
    return _best_single_sequence_match(query, tuple(sorted(targets.items())))


@lru_cache(maxsize=8192)
def _best_single_sequence_match(
    query: str, target_items: tuple[tuple[str, str], ...]
) -> SequenceMatch | None:
    targets = dict(target_items)
    scores = target_scores(query, targets)
    if not scores or max(scores.values()) <= 0:
        return None
    top = max(scores, key=lambda name: (scores[name], name))
    cutoff = explaining_min_score(len(query))
    traced = {name for name, score in scores.items() if score >= cutoff} | {top}
    best: SequenceMatch | None = None
    for name in sorted(traced):
        match = SequenceMatch(name, local_align(query, targets[name]), len(query))
        key = (match.explains, match.alignment.query_span * match.alignment.identity)
        if best is None or key > (
            best.explains,
            best.alignment.query_span * best.alignment.identity,
        ):
            best = match
    return best


@dataclass(frozen=True)
class GenomicSimilarity:
    window: int
    alignment: LocalAlignment
    orientation: str

    @property
    def high(self) -> bool:
        return (
            self.alignment.columns >= GENOMIC_MIN_LENGTH
            and self.alignment.identity >= GENOMIC_MIN_IDENTITY
        )


def breakpoint_windows(
    sequences: SequenceSource, chromosome: str, position: int
) -> dict[int, str | None]:
    """Reference ``+/-w`` bp around ``position`` for every ``w`` in :data:`GENOMIC_WINDOWS`.

    One request for the widest window; narrower windows are sliced from it
    around the breakpoint base (windows are clipped at the chromosome start).
    """
    widest = max(GENOMIC_WINDOWS)
    start = max(1, position - widest)
    sequence = sequences.region(chromosome, start, position + widest)
    if not sequence:
        return dict.fromkeys(GENOMIC_WINDOWS)
    centre = position - start
    return {
        window: sequence[max(0, centre - window) : centre + window + 1]
        for window in GENOMIC_WINDOWS
    }


@lru_cache(maxsize=8192)
def breakpoint_similarity(five_window: str, three_window: str, window: int) -> GenomicSimilarity:
    """Best local alignment of two breakpoint windows in either orientation."""
    same = local_align(five_window, three_window)
    opposite = local_align(five_window, reverse_complement(three_window))
    if opposite.score > same.score:
        return GenomicSimilarity(window, opposite, "opposite_strand")
    return GenomicSimilarity(window, same, "same_strand")


# --- data sources -----------------------------------------------------------------

ENSEMBL_REST = {"GRCh37": "https://grch37.rest.ensembl.org", "GRCh38": "https://rest.ensembl.org"}
UCSC_ASSEMBLY = {"GRCh37": "hg19", "GRCh38": "hg38"}
UCSC_DOWNLOAD = "https://hgdownload.soe.ucsc.edu/goldenPath/{assembly}/database/{table}.txt.gz"
_RETRY = {429, 500, 502, 503, 504}


def chromosome_name(value: Any) -> str | None:
    """``chr7``/``7`` -> ``7``; ``chrM`` -> ``MT``."""
    if value is None or str(value).strip() == "":
        return None
    name = str(value).strip()
    name = name[3:] if name.lower().startswith("chr") else name
    return "MT" if name.upper() == "M" else name


@dataclass(frozen=True)
class ParalogLookup:
    """Paralogue CDS for one gene, and whether the lookup is complete.

    ``complete`` is False when the gene could not be resolved, the lookup
    failed (e.g. HTTP 400/404), or a paralogue's CDS is missing; an empty
    ``sequences`` with ``complete`` True is a verified "no paralogues".
    """

    sequences: dict[str, str]
    complete: bool
    missing: tuple[str, ...] = ()


class SequenceSource(Protocol):
    genome_build: str

    def region(self, chromosome: str, start: int, end: int) -> str | None: ...

    def paralog_lookup(self, gene: str, transcript_id: str | None = None) -> ParalogLookup: ...

    def transcript_cds(self, transcript_id: str) -> str | None: ...


class EnsemblReleaseFiles:
    """Paralogue and transcript CDS from one Ensembl release's published files.

    The same data :class:`EnsemblSequenceSource` fetches over REST, read from
    ``Homo_sapiens.<assembly>.cds.all.fa.gz``, ``...pep.all.fa.gz`` (protein to
    transcript) and the human ``Compara.<release>.protein_default.homologies.tsv.gz``.
    Paralogues are human-human pairs whose homology type contains ``paralog``;
    each pair's CDS is that of the Compara member protein's transcript, as the
    REST ``sequence=cdna`` view returns it. A symbol maps to primary-assembly
    genes only (patch and haplotype scaffolds are ignored); a caller that
    knows the transcript passes it so its gene is used instead of the symbol.
    """

    def __init__(self, cds_fasta: str | Path, pep_fasta: str | Path, homologies: str | Path):
        self.cds: dict[str, str] = {}
        self.gene_of_transcript: dict[str, str] = {}
        self.genes_of_symbol: dict[str, set[str]] = defaultdict(set)
        for header, sequence in _fasta(cds_fasta):
            fields = header.split()
            transcript = _unversioned(fields[0])
            if "_PAR_Y" in fields[0]:
                continue
            self.cds[transcript] = sequence
            tags = dict(field.split(":", 1) for field in fields[2:] if ":" in field)
            gene = _unversioned(tags.get("gene", ""))
            self.gene_of_transcript[transcript] = gene
            if (
                fields[1] == "cds"
                and fields[2].startswith("chromosome:")
                and tags.get("gene_symbol")
            ):
                self.genes_of_symbol[tags["gene_symbol"]].add(gene)
        self.transcript_of_protein: dict[str, str] = {}
        for header, _ in _fasta(pep_fasta):
            fields = header.split()
            tags = dict(field.split(":", 1) for field in fields[2:] if ":" in field)
            if "_PAR_Y" not in fields[0] and tags.get("transcript"):
                self.transcript_of_protein[_unversioned(fields[0])] = _unversioned(
                    tags["transcript"]
                )
        self.paralogues: dict[str, dict[str, str]] = defaultdict(dict)
        self.member_protein: dict[str, str] = {}
        with gzip.open(homologies, "rt", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if (
                    row["species"] != "homo_sapiens"
                    or row["homology_species"] != "homo_sapiens"
                    or "paralog" not in row["homology_type"]
                ):
                    continue
                gene, other = row["gene_stable_id"], row["homology_gene_stable_id"]
                self.member_protein[gene] = row["protein_stable_id"]
                self.member_protein[other] = row["homology_protein_stable_id"]
                self.paralogues[gene][other] = row["homology_protein_stable_id"]
                self.paralogues[other][gene] = row["protein_stable_id"]

    def _protein_cds(self, protein: str) -> str | None:
        transcript = self.transcript_of_protein.get(_unversioned(protein))
        return self.cds.get(transcript) if transcript else None

    def gene_ids(self, gene: str, transcript_id: str | None = None) -> set[str]:
        if transcript_id:
            known = self.gene_of_transcript.get(_unversioned(transcript_id))
            if known:
                return {known}
        return set(self.genes_of_symbol.get(gene, ()))

    def paralog_lookup(self, gene: str, transcript_id: str | None = None) -> ParalogLookup:
        """Same keys as :meth:`EnsemblSequenceSource.paralog_lookup`.

        Incomplete when the gene does not resolve to an Ensembl gene, or when a
        Compara member protein's CDS is absent from the release files.
        """
        gene_ids = self.gene_ids(gene, transcript_id)
        if not gene_ids:
            return ParalogLookup({}, False, (f"{gene} not in the Ensembl release files",))
        result: dict[str, str] = {}
        missing: list[str] = []
        for gene_id in sorted(gene_ids):
            own = self.member_protein.get(gene_id)
            if own:
                own_cds = self._protein_cds(own)
                if own_cds:
                    result[f"{gene}:compara_canonical:{gene_id}"] = own_cds
                else:
                    missing.append(f"no CDS for {gene} Compara member {own}")
            for other, protein in sorted(self.paralogues.get(gene_id, {}).items()):
                sequence = self._protein_cds(protein)
                if sequence:
                    result[f"{gene}:paralog:{other}"] = sequence
                else:
                    missing.append(f"no CDS for {gene} paralogue {other} ({protein})")
        return ParalogLookup(result, not missing, tuple(missing))

    def paralog_cds(self, gene: str, transcript_id: str | None = None) -> dict[str, str]:
        return self.paralog_lookup(gene, transcript_id).sequences

    def transcript_cds(self, transcript_id: str) -> str | None:
        return self.cds.get(_unversioned(transcript_id))


def _unversioned(identifier: str) -> str:
    return identifier.split(".")[0]


def _fasta(path: str | Path) -> Iterator[tuple[str, str]]:
    header: str | None = None
    parts: list[str] = []
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts)
                header, parts = line[1:].strip(), []
            else:
                parts.append(line.strip())
    if header is not None:
        yield header, "".join(parts)


class TwoBitReference:
    """Read-only UCSC ``.2bit`` genome: forward-strand bases of a region.

    Implements the published format (header, sequence index, then per
    sequence its size, N blocks, soft-mask blocks and 2-bit packed bases).
    Soft-masking is ignored, so bases are upper case like Ensembl's; N blocks
    read as ``N``. Chromosomes are matched as ``chrN`` or ``N`` (``MT`` ->
    ``chrM``).
    """

    SIGNATURE = 0x1A412743

    def __init__(self, path: str | Path):
        self.path = Path(path)
        self._handle = self.path.open("rb")
        self._lock = threading.Lock()
        header = self._handle.read(16)
        for order in ("<", ">"):
            signature, version, count, _ = struct.unpack(order + "IIII", header)
            if signature == self.SIGNATURE:
                break
        else:
            raise ValueError(f"{self.path} is not a 2bit file")
        if version not in (0, 1):
            raise ValueError(f"{self.path}: unsupported 2bit version {version}")
        self._order = order
        offset_format = order + ("Q" if version == 1 else "I")
        offset_size = struct.calcsize(offset_format)
        self._offsets: dict[str, int] = {}
        for _ in range(count):
            (size,) = struct.unpack("B", self._handle.read(1))
            name = self._handle.read(size).decode("ascii")
            (self._offsets[name],) = struct.unpack(offset_format, self._handle.read(offset_size))
        self._layout: dict[str, tuple[int, list[tuple[int, int]], int]] = {}

    def _name(self, chromosome: str) -> str | None:
        bare = chromosome_name(chromosome)
        if bare is None:
            return None
        for name in (f"chr{'M' if bare == 'MT' else bare}", bare):
            if name in self._offsets:
                return name
        return None

    def _read_ints(self, count: int) -> list[int]:
        return list(struct.unpack(f"{self._order}{count}I", self._handle.read(4 * count)))

    def _sequence_layout(self, name: str) -> tuple[int, list[tuple[int, int]], int]:
        """``(size, N blocks as (start, end), offset of the packed bases)``."""
        if name not in self._layout:
            self._handle.seek(self._offsets[name])
            size, n_count = self._read_ints(2)
            n_starts, n_sizes = self._read_ints(n_count), self._read_ints(n_count)
            (mask_count,) = self._read_ints(1)
            self._handle.seek(8 * mask_count + 4, 1)  # mask starts + sizes, reserved
            blocks = [(start, start + length) for start, length in zip(n_starts, n_sizes)]
            self._layout[name] = (size, blocks, self._handle.tell())
        return self._layout[name]

    def region(self, chromosome: str, start: int, end: int) -> str | None:
        """Bases ``start..end`` (1-based, inclusive), clipped to the sequence."""
        name = self._name(chromosome)
        if name is None:
            return None
        with self._lock:
            size, n_blocks, packed = self._sequence_layout(name)
            begin, stop = max(0, start - 1), min(size, end)
            if begin >= stop:
                return None
            self._handle.seek(packed + begin // 4)
            data = self._handle.read((stop - 1) // 4 - begin // 4 + 1)
        bases = "".join(_TWO_BIT_BYTE[byte] for byte in data)
        first = begin - (begin // 4) * 4
        sequence = list(bases[first : first + stop - begin])
        for block_start, block_end in n_blocks:
            for index in range(max(block_start, begin), min(block_end, stop)):
                sequence[index - begin] = "N"
        return "".join(sequence)


_TWO_BIT_BYTE = [
    "".join("TCAG"[(byte >> shift) & 3] for shift in (6, 4, 2, 0)) for byte in range(256)
]


PARALOG_CACHE_PREFIX = "paralog_lookup_v2_"


def parse_homology_response(gene: str, data: Any) -> dict[str, Any]:
    """Paralogue CDS from an Ensembl REST ``/homology/symbol`` payload, fail-closed.

    Complete only when the payload has a non-empty ``data`` list whose every
    entry names the resolved gene (``id``) and carries an explicit
    ``homologies`` list, and every homology has a source and target with an id
    and a sequence. An empty ``homologies`` list on a resolved gene is the one
    verified "no paralogues". Anything else (no payload, ``{}``, ``data: []``,
    an entry without ``homologies``, a member without a sequence) is incomplete.
    """
    missing: list[str] = []
    result: dict[str, str] = {}
    entries = data.get("data") if isinstance(data, dict) else None
    if data is None:
        missing.append(f"Ensembl homology lookup for {gene} returned 400/404")
    elif not isinstance(entries, list) or not entries:
        missing.append(f"Ensembl homology response for {gene} resolved no gene")
    else:
        for entry in entries:
            if not isinstance(entry, dict) or not entry.get("id"):
                missing.append(f"Ensembl homology entry for {gene} names no gene")
                continue
            homologies = entry.get("homologies")
            if not isinstance(homologies, list):
                missing.append(f"Ensembl homology entry {entry['id']} for {gene} has no list")
                continue
            for homology in homologies:
                if not isinstance(homology, dict):
                    missing.append(f"malformed homology for {gene}")
                    continue
                for role, label in (("source", "compara_canonical"), ("target", "paralog")):
                    member = homology.get(role)
                    if isinstance(member, dict) and member.get("id") and member.get("seq"):
                        result[f"{gene}:{label}:{member['id']}"] = member["seq"]
                    else:
                        member_id = member.get("id") if isinstance(member, dict) else None
                        missing.append(f"no sequence for {gene} {label} {member_id}")
    return {"sequences": result, "complete": not missing, "missing": missing}


class EnsemblSequenceSource:
    """Reference regions, paralogue CDS and liftover from Ensembl REST, cached per build."""

    def __init__(
        self,
        genome_build: str,
        cache_dir: str | Path,
        *,
        session: requests.Session | None = None,
        max_retries: int = 5,
        backoff_seconds: float = 1.0,
        timeout: float = 120,
        reference: TwoBitReference | None = None,
        release_files: EnsemblReleaseFiles | None = None,
    ):
        """``reference`` serves :meth:`region` and ``release_files`` serves
        :meth:`paralog_cds` / :meth:`transcript_cds` locally instead of REST."""
        build = normalize_genome_build(genome_build)
        if build not in ENSEMBL_REST:
            raise ValueError(f"no Ensembl REST host for genome build {genome_build!r}")
        self.genome_build = build
        self.host = ENSEMBL_REST[build]
        self.cache_dir = Path(cache_dir) / "ensembl_sequence" / build
        check_cache_provenance(self.cache_dir, build, self.host)
        self.session = session or requests.Session()
        self.max_retries = max_retries
        self.backoff_seconds = backoff_seconds
        self.timeout = timeout
        self.reference = reference
        self.release_files = release_files
        self._memory: dict[str, Any] = {}

    def _cached(self, key: str, fetch: Any) -> Any:
        if key in self._memory:
            return self._memory[key]
        path = self.cache_dir / f"{key}.json"
        if path.exists():
            value = json.loads(path.read_text())["value"]
        else:
            value = fetch()
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps({"value": value}))
        self._memory[key] = value
        return value

    def _get(self, path: str, params: Mapping[str, str] | None = None) -> Any | None:
        for attempt in range(self.max_retries + 1):
            try:
                response = self.session.get(
                    f"{self.host}{path}",
                    params=params,
                    headers={"Content-Type": "application/json"},
                    timeout=self.timeout,
                )
            except TRANSIENT_REQUEST_ERRORS:
                if attempt == self.max_retries:
                    raise
                time.sleep(self.backoff_seconds * 2**attempt)
                continue
            if response.status_code not in _RETRY or attempt == self.max_retries:
                break
            wait = float(response.headers.get("Retry-After", 0) or 0)
            time.sleep(max(wait, self.backoff_seconds * 2**attempt))
        if response.status_code in (400, 404):
            return None
        response.raise_for_status()
        return response.json()

    def region(self, chromosome: str, start: int, end: int) -> str | None:
        """Forward-strand reference bases ``start..end`` (1-based, inclusive)."""
        if self.reference is not None:
            return self.reference.region(chromosome, max(1, start), end)
        chromosome = chromosome_name(chromosome) or ""
        start = max(1, start)

        def fetch() -> str | None:
            data = self._get(
                f"/sequence/region/human/{chromosome}:{start}..{end}:1",
                {"coord_system_version": self.genome_build},
            )
            return data.get("seq") if data else None

        return self._cached(f"region_{chromosome}_{start}_{end}", fetch)

    def paralog_lookup(self, gene: str, transcript_id: str | None = None) -> ParalogLookup:
        """``{"<gene>:paralog:<ENSG>": cds, "<gene>:compara_canonical:<ENSG>": cds}``.

        REST looks the gene up by symbol; release files use ``transcript_id``'s
        gene when given. A REST 400/404 (unknown symbol) is an incomplete
        lookup, not an empty one; so is a homology entry without a sequence.
        """
        if self.release_files is not None:
            return self.release_files.paralog_lookup(gene, transcript_id)

        def fetch() -> dict[str, Any]:
            data = self._get(
                f"/homology/symbol/human/{gene}",
                {"type": "paralogues", "sequence": "cdna", "aligned": "0"},
            )
            return parse_homology_response(gene, data)

        # Versioned key: entries cached under earlier completeness rules
        # (``paralogs_<gene>``, ``paralog_lookup_<gene>``) are never reused.
        value = self._cached(f"{PARALOG_CACHE_PREFIX}{gene}", fetch)
        return ParalogLookup(value["sequences"], value["complete"], tuple(value["missing"]))

    def paralog_cds(self, gene: str, transcript_id: str | None = None) -> dict[str, str]:
        return self.paralog_lookup(gene, transcript_id).sequences

    def transcript_cds(self, transcript_id: str) -> str | None:
        if self.release_files is not None:
            return self.release_files.transcript_cds(transcript_id)

        def fetch() -> str | None:
            data = self._get(f"/sequence/id/{transcript_id}", {"type": "cds"})
            return data.get("seq") if data else None

        return self._cached(f"cds_{transcript_id}", fetch)

    def liftover_to_grch38(self, chromosome: str, position: int) -> tuple[str, int] | None:
        """One GRCh37 base on GRCh38, or ``None`` unless it maps to exactly one base."""
        if self.genome_build != "GRCh37":
            raise ValueError("liftover is only defined from a GRCh37 source")
        chromosome = chromosome_name(chromosome) or ""

        def fetch() -> list[Any] | None:
            data = self._get(f"/map/human/GRCh37/{chromosome}:{position}..{position}:1/GRCh38")
            mappings = (data or {}).get("mappings", [])
            if len(mappings) != 1:
                return None
            mapped = mappings[0]["mapped"]
            if mapped["start"] != mapped["end"]:
                return None
            return [str(mapped["seq_region_name"]), int(mapped["start"])]

        value = self._cached(f"liftover_{chromosome}_{position}", fetch)
        return (value[0], value[1]) if value else None


@dataclass(frozen=True)
class TrackHit:
    track: str
    name: str
    detail: str


class RepeatTracks:
    """UCSC ``genomicSuperDups`` and ``rmsk`` overlaps at chosen positions.

    The full tables are streamed once and only records covering a queried
    base are kept. UCSC coordinates are 0-based half-open, so base ``p``
    (1-based) is covered when ``start < p <= end``.
    """

    def __init__(self, genome_build: str, superdups: Path, rmsk: Path):
        self.genome_build = genome_build
        self.paths = {"segdup": superdups, "rmsk": rmsk}

    def overlaps(
        self, positions: Iterable[tuple[str, int]]
    ) -> dict[tuple[str, int], list[TrackHit]]:
        wanted: dict[str, list[int]] = defaultdict(list)
        for chromosome, position in positions:
            name = chromosome_name(chromosome)
            if name is not None:
                wanted[name].append(position)
        for name in wanted:
            wanted[name] = sorted(set(wanted[name]))
        hits: dict[tuple[str, int], list[TrackHit]] = defaultdict(list)
        for track, path in self.paths.items():
            for chromosome, start, end, hit in _track_records(track, path):
                points = wanted.get(chromosome, [])
                if not points:
                    continue
                index = bisect_left(points, start + 1)
                while index < len(points) and points[index] <= end:
                    hits[(chromosome, points[index])].append(hit)
                    index += 1
        return hits


def _track_records(track: str, path: Path) -> Iterator[tuple[str, int, int, TrackHit]]:
    with gzip.open(path, "rt", newline="") as handle:
        for fields in csv.reader(handle, delimiter="\t"):
            if track == "segdup":
                chromosome, start, end = fields[1], int(fields[2]), int(fields[3])
                detail = (
                    f"other={fields[7]}:{int(fields[8]) + 1}-{fields[9]} "
                    f"fracMatch={float(fields[26]):.3f}"
                )
                hit = TrackHit("segdup", fields[4], detail)
            else:
                chromosome, start, end = fields[5], int(fields[6]), int(fields[7])
                hit = TrackHit("rmsk", fields[10], f"{fields[11]}/{fields[12]}")
            name = chromosome_name(chromosome)
            if name is not None:
                yield name, start, end, hit


# --- normal-recurrence index --------------------------------------------------------

JunctionKey = tuple[str, int, str, int]


@dataclass
class JunctionObservation:
    gene_5p: str
    gene_3p: str
    sample_id: str
    tissue_class: str
    tissue_label: str


def tcga_gtex_index(
    path: str | Path,
    keys: set[JunctionKey],
    gene_pairs: set[tuple[str, str]],
) -> tuple[dict[JunctionKey, list[JunctionObservation]], dict[tuple[str, str], Counter[str]]]:
    """Exact-junction observations and per-gene-pair tissue-class counts.

    ``path`` is the normalized TCGA+GTEx TSV. Only the requested junction keys
    and gene pairs are kept.
    """
    exact: dict[JunctionKey, list[JunctionObservation]] = defaultdict(list)
    pairs: dict[tuple[str, str], Counter[str]] = defaultdict(Counter)
    with Path(path).open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            pair = (row["gene_5p"], row["gene_3p"])
            if pair in gene_pairs:
                pairs[pair][row["tissue_class"]] += 1
            five, three = chromosome_name(row["chrom_5p"]), chromosome_name(row["chrom_3p"])
            if not five or not three or not row["pos_5p"] or not row["pos_3p"]:
                continue
            key = (five, int(row["pos_5p"]), three, int(row["pos_3p"]))
            if key in keys:
                label = json.loads(row["annotation_flags"]).get("tissue_subtype", "")
                exact[key].append(
                    JunctionObservation(
                        row["gene_5p"], row["gene_3p"], row["sample_id"], row["tissue_class"], label
                    )
                )
    return exact, pairs


# --- per-row flags ------------------------------------------------------------------


def _int(value: Any) -> int | None:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def paralog_flag(row: Mapping[str, Any], sequences: SequenceSource) -> dict[str, Any]:
    window = row.get("junction_context_nt") or ""
    if not window:
        return {
            "paralog_sequence_flag": NOT_ASSESSED,
            "paralog_sequence_reason": "no fused CDS at the junction "
            f"({row.get('reconstitution_class')}); no junction sequence to test",
        }
    targets: dict[str, str] = {}
    gaps: list[str] = []  # missing data: the test is incomplete
    notes: list[str] = []  # verified facts, e.g. a gene with no paralogues
    for side in ("five_prime", "three_prime"):
        gene = row[f"{side}_gene"]
        transcript = row.get(f"{side}_ensembl_transcript")
        if transcript:
            cds = sequences.transcript_cds(transcript)
            if cds:
                targets[f"{gene}:parent:{transcript}"] = cds
            else:
                gaps.append(f"no CDS for {gene} parent {transcript}")
        else:
            gaps.append(f"no Ensembl transcript for {gene}")
        lookup = sequences.paralog_lookup(gene, transcript or None)
        if not lookup.complete:
            gaps.extend(lookup.missing or (f"paralogue lookup for {gene} incomplete",))
        elif not any(":paralog:" in name for name in lookup.sequences):
            notes.append(f"{gene} has no Compara paralogues")
        targets.update(lookup.sequences)
    gap_note = f"; incomplete: {'; '.join(gaps)}" if gaps else ""
    note = f" ({'; '.join(notes)})" if notes else ""
    if not targets:
        return {
            "paralog_sequence_flag": NOT_ASSESSED,
            "paralog_sequence_reason": "no parent or paralogue CDS available" + note + gap_note,
        }
    match = best_single_sequence_match(window, targets)
    if match is not None and match.explains:
        flag = FLAGGED  # an available sequence already explains the junction
    elif gaps:
        flag = NOT_ASSESSED  # a missing sequence might have explained it
    else:
        flag = PASS
    if match is None:
        return {
            "paralog_sequence_flag": flag,
            "paralog_sequence_reason": f"no local alignment to any of {len(targets)} sequences"
            + note
            + gap_note,
        }
    summary = (
        f"{match.target} {match.alignment.identity:.1%} identity over "
        f"{match.alignment.query_span}/{len(window)} nt"
    )
    return {
        "paralog_sequence_flag": flag,
        "paralog_sequence_reason": (
            f"junction window explained by one sequence: {summary}"
            if match.explains
            else f"best single sequence ({len(targets)} tested): {summary}{note}{gap_note}"
        ),
        "paralog_best_target": match.target,
        "paralog_best_identity": round(match.alignment.identity, 4),
        "paralog_best_query_coverage": round(match.coverage, 4),
    }


def genomic_flags(
    row: Mapping[str, Any],
    sequences: SequenceSource,
    repeat_hits: Mapping[tuple[str, int], list[TrackHit]] | None,
    sample_calls: int,
    repeat_unavailable_reason: str | None = None,
) -> dict[str, Any]:
    out: dict[str, Any] = {}
    parts: list[str] = []
    states: list[str] = []
    five_chrom = chromosome_name(row.get("five_prime_chromosome"))
    three_chrom = chromosome_name(row.get("three_prime_chromosome"))
    five_pos = _int(row.get("five_prime_breakpoint"))
    three_pos = _int(row.get("three_prime_breakpoint"))
    located = None not in (five_chrom, three_chrom, five_pos, three_pos)

    if not located:
        states.append(NOT_ASSESSED)
        parts.append("similarity not assessed: a breakpoint chromosome or position is missing")
    else:
        assert five_chrom and three_chrom and five_pos and three_pos
        high = []
        five_windows = breakpoint_windows(sequences, five_chrom, five_pos)
        three_windows = breakpoint_windows(sequences, three_chrom, three_pos)
        for window in GENOMIC_WINDOWS:
            five_seq, three_seq = five_windows[window], three_windows[window]
            if not five_seq or not three_seq:
                states.append(NOT_ASSESSED)
                parts.append(f"+/-{window} bp reference unavailable")
                continue
            similarity = breakpoint_similarity(five_seq, three_seq, window)
            out[f"genomic_identity_{window}"] = round(similarity.alignment.identity, 4)
            out[f"genomic_aligned_length_{window}"] = similarity.alignment.columns
            if similarity.high:
                high.append(
                    f"+/-{window} bp {similarity.alignment.identity:.1%} over "
                    f"{similarity.alignment.columns} bp ({similarity.orientation})"
                )
        if high:
            states.append(FLAGGED)
            parts.append("breakpoint regions similar: " + "; ".join(high))
        elif not any(s == NOT_ASSESSED for s in states):
            states.append(PASS)

    if repeat_hits is None or not located:
        out["breakpoint_repeat_flag"] = NOT_ASSESSED
        out["breakpoint_repeat_reason"] = repeat_unavailable_reason or "breakpoint not located"
        states.append(NOT_ASSESSED)
        parts.append("segdup/repeat not assessed")
    else:
        assert five_chrom and three_chrom and five_pos and three_pos
        found = []
        for label, key in (("5'", (five_chrom, five_pos)), ("3'", (three_chrom, three_pos))):
            for hit in repeat_hits.get(key, []):
                text = f"{label} in {hit.track} {hit.name} ({hit.detail})"
                if text not in found:
                    found.append(text)
        out["breakpoint_repeat_flag"] = FLAGGED if found else PASS
        out["breakpoint_repeat_reason"] = "; ".join(found) or "neither breakpoint base is in a "
        if not found:
            out["breakpoint_repeat_reason"] += "segmental duplication or RepeatMasker element"
        if found:
            states.append(FLAGGED)
            parts.append("; ".join(found))

    out["sample_shared_family_calls"] = sample_calls
    out["sample_burden_flag"] = FLAGGED if sample_calls >= SAMPLE_BURDEN_MIN else PASS
    if sample_calls >= SAMPLE_BURDEN_MIN:
        states.append(FLAGGED)
        parts.append(f"sample has {sample_calls} distinct shared-family calls")

    split = _int(row.get("tumor_split_reads"))
    precision = str(row.get("breakpoint_precision") or "")
    weak = []
    if split == 0:
        weak.append(f"0 split/junction reads ({_int(row.get('tumor_paired_end_reads'))} pairs)")
    if precision == "IMPRECISE":
        weak.append("imprecise breakpoint")
    out["read_evidence_flag"] = FLAGGED if weak else (PASS if split is not None else NOT_ASSESSED)
    out["read_evidence_reason"] = (
        "; ".join(weak) if weak else (f"{split} split reads" if split is not None else "no counts")
    )
    if weak:
        states.append(FLAGGED)
        parts.append("; ".join(weak))
    elif split is None:
        states.append(NOT_ASSESSED)
        parts.append("read support not assessed: no split-read count")

    if FLAGGED in states:
        flag = FLAGGED
    elif NOT_ASSESSED in states:
        flag = NOT_ASSESSED
    else:
        flag = PASS
    out["genomic_similarity_flag"] = flag
    out["genomic_similarity_reason"] = "; ".join(parts) or (
        f"breakpoint regions <{GENOMIC_MIN_IDENTITY:.0%} over >={GENOMIC_MIN_LENGTH} bp; "
        "no segdup/repeat; single call in sample; split reads present"
    )
    return out


def normal_flag(
    row: Mapping[str, Any],
    key: JunctionKey | None,
    key_note: str | None,
    exact: Mapping[JunctionKey, list[JunctionObservation]],
    pair_counts: Mapping[tuple[str, str], Counter[str]],
) -> dict[str, Any]:
    pair = (row["five_prime_gene"], row["three_prime_gene"])
    counts = pair_counts.get(pair, Counter())
    pair_note = (
        f"gene pair in TCGA+GTEx: {counts.get('normal', 0)} normal, "
        f"{counts.get('cell_line', 0)} cell-line, {counts.get('tumor', 0)} tumor rows"
    )
    out: dict[str, Any] = {"normal_gene_pair_samples": counts.get("normal", 0)}
    own_normal = (
        str(row.get("Open_source_table") or "").startswith("TCGA_n_GTEx")
        and row.get("Open_tissue_class") == "normal"
    )
    if key is None:
        out["normal_recurrence_flag"] = FLAGGED if own_normal else NOT_ASSESSED
        out["normal_recurrence_reason"] = (
            "call is itself from a TCGA/GTEx normal sample; " if own_normal else ""
        ) + f"exact junction not assessed: {key_note}; {pair_note}"
        return out
    observations = [o for o in exact.get(key, []) if o.sample_id != row.get("sample_id")]
    same_genes = [o for o in observations if (o.gene_5p, o.gene_3p) == pair]
    other_genes = [o for o in observations if (o.gene_5p, o.gene_3p) != pair]
    normals = [o for o in same_genes if o.tissue_class == "normal"]
    cells = [o for o in same_genes if o.tissue_class == "cell_line"]
    out["normal_exact_junction_samples"] = len(normals)
    out["cell_line_exact_junction_samples"] = len(cells)
    reasons = []
    if own_normal:
        reasons.append("call is itself from a TCGA/GTEx normal sample")
    if normals:
        tissues = Counter(o.tissue_label or "TCGA normal" for o in normals)
        reasons.append(
            f"exact junction in {len(normals)} other normal sample(s): "
            + ", ".join(f"{t} x{n}" for t, n in tissues.most_common(4))
        )
    if cells:
        reasons.append(f"exact junction in {len(cells)} GTEx cell-line sample(s) (not counted)")
    other_normals = [o for o in other_genes if o.tissue_class == "normal"]
    if other_normals:
        names = sorted({f"{o.gene_5p}::{o.gene_3p}" for o in other_normals})
        reasons.append(
            f"same coordinates in {len(other_normals)} normal sample(s) under other gene "
            f"symbols {', '.join(names[:3])} (not counted)"
        )
    if key_note:
        reasons.append(key_note)
    reasons.append(pair_note)
    out["normal_recurrence_flag"] = FLAGGED if normals or own_normal else PASS
    if not normals and not own_normal:
        reasons.insert(0, "exact junction absent from all TCGA/GTEx normal samples")
    out["normal_recurrence_reason"] = "; ".join(reasons)
    return out


def sample_call_counts(rows: Sequence[Mapping[str, Any]]) -> Counter[tuple[str, str]]:
    """Distinct shared-family events (pair + breakpoints) per (cohort, sample)."""
    events = {
        (
            str(row.get("cohort")),
            str(row.get("sample_id")),
            row["five_prime_gene"],
            row["three_prime_gene"],
            row.get("five_prime_breakpoint"),
            row.get("three_prime_breakpoint"),
        )
        for row in rows
        if row.get("reconstitution_class") != "no_shared_family"
    }
    return Counter((event[0], event[1]) for event in events)


Liftover = Any
"""Callable ``(chromosome, position) -> (chromosome, position) | None`` onto GRCh38."""


def junction_key(
    row: Mapping[str, Any], liftover: Liftover | None
) -> tuple[JunctionKey | None, str | None]:
    """The row's GRCh38 junction key, lifting GRCh37 breakpoints over base by base."""
    five = (
        chromosome_name(row.get("five_prime_chromosome")),
        _int(row.get("five_prime_breakpoint")),
    )
    three = (
        chromosome_name(row.get("three_prime_chromosome")),
        _int(row.get("three_prime_breakpoint")),
    )
    if None in (*five, *three):
        return None, "a breakpoint chromosome or position is missing"
    build = normalize_genome_build(row.get("genome_build"))
    if build == "GRCh38":
        return (five[0], five[1], three[0], three[1]), None  # type: ignore[return-value]
    if build != "GRCh37" or liftover is None:
        return None, f"no liftover from {build} to GRCh38"
    lifted_five, lifted_three = liftover(*five), liftover(*three)
    if lifted_five is None or lifted_three is None:
        return None, "a breakpoint does not lift over to a single GRCh38 base"
    note = (
        f"GRCh37 breakpoints lifted to GRCh38 {lifted_five[0]}:{lifted_five[1]} / "
        f"{lifted_three[0]}:{lifted_three[1]}"
    )
    return (lifted_five[0], lifted_five[1], lifted_three[0], lifted_three[1]), note


UNFLAGGED_CLASSES = frozenset({"no_transcript_model"})
"""Shared-family classes that never reach junction translation and are not flagged."""


def prefetch_sequences(
    rows: Sequence[Mapping[str, Any]],
    sequences: Mapping[str, SequenceSource],
    *,
    workers: int,
    progress: Any = None,
) -> None:
    """Warm each source's cache with every fetch the rows' flags will make."""
    tasks: set[tuple[str, str, tuple[Any, ...]]] = set()
    for row in rows:
        build = normalize_genome_build(row.get("genome_build")) or ""
        if build not in sequences:
            continue
        for side in ("five_prime", "three_prime"):
            if row.get("junction_context_nt"):
                tasks.add(
                    (
                        build,
                        "paralog_lookup",
                        (row[f"{side}_gene"], row.get(f"{side}_ensembl_transcript") or None),
                    )
                )
                transcript = row.get(f"{side}_ensembl_transcript")
                if transcript:
                    tasks.add((build, "transcript_cds", (transcript,)))
            chromosome = chromosome_name(row.get(f"{side}_chromosome"))
            position = _int(row.get(f"{side}_breakpoint"))
            if chromosome and position:
                tasks.add((build, "windows", (chromosome, position)))

    def run(task: tuple[str, str, tuple[Any, ...]]) -> None:
        build, kind, args = task
        source = sequences[build]
        if kind == "windows":
            breakpoint_windows(source, *args)
        else:
            getattr(source, kind)(*args)

    ordered = sorted(tasks, key=repr)
    if progress is not None:
        progress(f"prefetching {len(ordered)} sequence requests with {workers} workers")
    if workers <= 1:
        for task in ordered:
            run(task)
        return
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for done, _ in enumerate(pool.map(run, ordered), start=1):
            if progress is not None and done % 500 == 0:
                progress(f"prefetched {done}/{len(ordered)}")


def annotate_rows(
    rows: Sequence[dict[str, Any]],
    *,
    sequences: Mapping[str, SequenceSource],
    repeat_tracks: Mapping[str, RepeatTracks | None],
    repeat_unavailable: Mapping[str, str],
    tcga_gtex_table: str | Path,
    liftover: Liftover | None = None,
    progress: Any = None,
    prefetch_workers: int = 4,
) -> list[dict[str, Any]]:
    """Return ``rows`` with every :data:`FLAG_COLUMNS` column filled.

    Rows of :data:`UNFLAGGED_CLASSES` (no transcript model, so no junction
    translation) and ``no_shared_family`` rows get ``not_applicable`` with the
    reason; they still count toward a sample's shared-family burden.
    ``prefetch_workers`` threads fetch reference windows and CDS before the
    rows are flagged (each fetch is cached).

    ``sequences`` and ``repeat_tracks`` are keyed by each row's
    ``genome_build``. A build without repeat tracks gets ``not_assessed`` with
    ``repeat_unavailable[build]`` as the reason. ``no_shared_family`` rows get
    ``not_applicable``.
    """
    candidates = [row for row in rows if row.get("reconstitution_class") != "no_shared_family"]
    flagged_rows = [
        row for row in candidates if row.get("reconstitution_class") not in UNFLAGGED_CLASSES
    ]
    keys = {id(row): junction_key(row, liftover) for row in flagged_rows}
    exact, pairs = tcga_gtex_index(
        tcga_gtex_table,
        {key for key, _ in keys.values() if key is not None},
        {(row["five_prime_gene"], row["three_prime_gene"]) for row in flagged_rows},
    )
    prefetch_sequences(flagged_rows, sequences, workers=prefetch_workers, progress=progress)
    hits: dict[str, Mapping[tuple[str, int], list[TrackHit]]] = {}
    for build, tracks in repeat_tracks.items():
        if tracks is None:
            continue
        positions = [
            (str(row.get(f"{side}_chromosome")), int(row[f"{side}_breakpoint"]))
            for row in flagged_rows
            if normalize_genome_build(row.get("genome_build")) == build
            for side in ("five_prime", "three_prime")
            if row.get(f"{side}_chromosome") and _int(row.get(f"{side}_breakpoint"))
        ]
        hits[build] = tracks.overlaps(positions)
    burden = sample_call_counts(candidates)

    annotated = []
    for index, row in enumerate(rows, start=1):
        out = dict(row)
        for column in FLAG_COLUMNS:
            out.setdefault(column, None)
        row_class = row.get("reconstitution_class")
        if row_class == "no_shared_family" or row_class in UNFLAGGED_CLASSES:
            reason = f"row did not reach junction translation ({row_class})"
            for flag in (
                "paralog_sequence_flag",
                "genomic_similarity_flag",
                "normal_recurrence_flag",
                "breakpoint_repeat_flag",
                "read_evidence_flag",
            ):
                out[flag] = NOT_APPLICABLE
                out[flag.replace("_flag", "_reason")] = reason
            out["sample_burden_flag"] = NOT_APPLICABLE
            annotated.append(out)
            continue
        build = normalize_genome_build(row.get("genome_build")) or ""
        source = sequences[build]
        out.update(paralog_flag(row, source))
        out.update(
            genomic_flags(
                row,
                source,
                hits.get(build),
                burden[(str(row.get("cohort")), str(row.get("sample_id")))],
                repeat_unavailable.get(build),
            )
        )
        key, note = keys[id(row)]
        out.update(normal_flag(row, key, note, exact, pairs))
        annotated.append(out)
        if progress is not None and index % 25 == 0:
            progress(f"flagged {index}/{len(rows)} rows")
    return annotated
