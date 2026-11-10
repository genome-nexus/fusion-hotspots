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
import time
from bisect import bisect_left
from collections import Counter, defaultdict
from collections.abc import Iterable, Iterator, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Protocol

import numpy as np
import requests

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

EXPLAINING_MIN_SCORE = 63
"""Lowest score of an alignment that explains a 100 nt window: >=90 columns at
>=90% identity scores at least 0.9c - 2 * 0.1c = 0.7c >= 63. Only targets whose
best local score reaches it (or the overall best target) are traced back."""


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
    target scoring at least :data:`EXPLAINING_MIN_SCORE`, and of the overall
    best target, is traced back. ``None`` when no target aligns at all.
    """
    scores = target_scores(query, targets)
    if not scores or max(scores.values()) <= 0:
        return None
    top = max(scores, key=lambda name: (scores[name], name))
    traced = {name for name, score in scores.items() if score >= EXPLAINING_MIN_SCORE} | {top}
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


class SequenceSource(Protocol):
    genome_build: str

    def region(self, chromosome: str, start: int, end: int) -> str | None: ...

    def paralog_cds(self, gene: str) -> dict[str, str]: ...

    def transcript_cds(self, transcript_id: str) -> str | None: ...


class EnsemblSequenceSource:
    """Reference regions, paralogue CDS and liftover from Ensembl REST, cached per build."""

    def __init__(
        self,
        genome_build: str,
        cache_dir: str | Path,
        *,
        session: requests.Session | None = None,
        max_retries: int = 4,
    ):
        build = normalize_genome_build(genome_build)
        if build not in ENSEMBL_REST:
            raise ValueError(f"no Ensembl REST host for genome build {genome_build!r}")
        self.genome_build = build
        self.host = ENSEMBL_REST[build]
        self.cache_dir = Path(cache_dir) / "ensembl_sequence" / build
        check_cache_provenance(self.cache_dir, build, self.host)
        self.session = session or requests.Session()
        self.max_retries = max_retries
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
            response = self.session.get(
                f"{self.host}{path}",
                params=params,
                headers={"Content-Type": "application/json"},
                timeout=60,
            )
            if response.status_code not in _RETRY or attempt == self.max_retries:
                break
            wait = float(response.headers.get("Retry-After", 0) or 0)
            time.sleep(max(wait, 0.5 * 2**attempt))
        if response.status_code in (400, 404):
            return None
        response.raise_for_status()
        return response.json()

    def region(self, chromosome: str, start: int, end: int) -> str | None:
        """Forward-strand reference bases ``start..end`` (1-based, inclusive)."""
        chromosome = chromosome_name(chromosome) or ""
        start = max(1, start)

        def fetch() -> str | None:
            data = self._get(
                f"/sequence/region/human/{chromosome}:{start}..{end}:1",
                {"coord_system_version": self.genome_build},
            )
            return data.get("seq") if data else None

        return self._cached(f"region_{chromosome}_{start}_{end}", fetch)

    def paralog_cds(self, gene: str) -> dict[str, str]:
        """``{"paralog:<ENSG> (<perc_id>%)": cds}`` for Compara paralogues of ``gene``."""

        def fetch() -> dict[str, str]:
            data = self._get(
                f"/homology/symbol/human/{gene}",
                {"type": "paralogues", "sequence": "cdna", "aligned": "0"},
            )
            result: dict[str, str] = {}
            for entry in (data or {}).get("data", []):
                for homology in entry.get("homologies", []):
                    source, target = homology.get("source", {}), homology.get("target", {})
                    if source.get("seq"):
                        result[f"{gene}:compara_canonical:{source.get('id')}"] = source["seq"]
                    if target.get("seq"):
                        result[f"{gene}:paralog:{target.get('id')}"] = target["seq"]
            return result

        return self._cached(f"paralogs_{gene}", fetch)

    def transcript_cds(self, transcript_id: str) -> str | None:
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
    for side in ("five_prime", "three_prime"):
        gene = row[f"{side}_gene"]
        transcript = row.get(f"{side}_ensembl_transcript")
        if transcript:
            cds = sequences.transcript_cds(transcript)
            if cds:
                targets[f"{gene}:parent:{transcript}"] = cds
        targets.update(sequences.paralog_cds(gene))
    if not targets:
        return {
            "paralog_sequence_flag": NOT_ASSESSED,
            "paralog_sequence_reason": "no parent or paralogue CDS available",
        }
    match = best_single_sequence_match(window, targets)
    if match is None:
        return {
            "paralog_sequence_flag": PASS,
            "paralog_sequence_reason": f"no local alignment to any of {len(targets)} sequences",
        }
    summary = (
        f"{match.target} {match.alignment.identity:.1%} identity over "
        f"{match.alignment.query_span}/{len(window)} nt"
    )
    return {
        "paralog_sequence_flag": FLAGGED if match.explains else PASS,
        "paralog_sequence_reason": (
            f"junction window explained by one sequence: {summary}"
            if match.explains
            else f"best single sequence ({len(targets)} tested): {summary}"
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
        for window in GENOMIC_WINDOWS:
            five_seq = sequences.region(five_chrom, five_pos - window, five_pos + window)
            three_seq = sequences.region(three_chrom, three_pos - window, three_pos + window)
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


def annotate_rows(
    rows: Sequence[dict[str, Any]],
    *,
    sequences: Mapping[str, SequenceSource],
    repeat_tracks: Mapping[str, RepeatTracks | None],
    repeat_unavailable: Mapping[str, str],
    tcga_gtex_table: str | Path,
    liftover: Liftover | None = None,
    progress: Any = None,
) -> list[dict[str, Any]]:
    """Return ``rows`` with every :data:`FLAG_COLUMNS` column filled.

    ``sequences`` and ``repeat_tracks`` are keyed by each row's
    ``genome_build``. A build without repeat tracks gets ``not_assessed`` with
    ``repeat_unavailable[build]`` as the reason. ``no_shared_family`` rows get
    ``not_applicable``.
    """
    candidates = [row for row in rows if row.get("reconstitution_class") != "no_shared_family"]
    keys = {id(row): junction_key(row, liftover) for row in candidates}
    exact, pairs = tcga_gtex_index(
        tcga_gtex_table,
        {key for key, _ in keys.values() if key is not None},
        {(row["five_prime_gene"], row["three_prime_gene"]) for row in candidates},
    )
    hits: dict[str, Mapping[tuple[str, int], list[TrackHit]]] = {}
    for build, tracks in repeat_tracks.items():
        if tracks is None:
            continue
        positions = [
            (str(row.get(f"{side}_chromosome")), int(row[f"{side}_breakpoint"]))
            for row in candidates
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
        if row.get("reconstitution_class") == "no_shared_family":
            for flag in (
                "paralog_sequence_flag",
                "genomic_similarity_flag",
                "normal_recurrence_flag",
            ):
                out[flag] = NOT_APPLICABLE
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
