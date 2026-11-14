"""What to score for each record, and the client that drives :mod:`.worker`.

For every record three sequences are scored with identical settings, the same
windows PR #160 predicted: the chimeric domain +/- ``FLANK`` residues
(:attr:`ChimeraRecord.chimeric_domain_flanked`) and each parent's Pfam domain
+/- ``FLANK`` (:meth:`ParentDomain.window`). Flanks are context only: the
masked positions are the domain residues, the junction window, and every
parent residue a chimera domain residue was copied from (for the paired
comparison in :mod:`.metrics`).
"""

from __future__ import annotations

import json
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping, Protocol

from cfh.domain_function.record import ChimeraRecord
from cfh.domain_plm.cache import ModelSpec, ScoreCache, valid_logp

CANONICAL = frozenset("ACDEFGHIKLMNPQRSTVWY")
MAX_RESIDUES = 1022
"""ESM-2 context limit used here (1024 tokens including BOS/EOS); longer sequences are scored
in per-position centred windows (see :func:`cfh.domain_plm.worker.context_windows`)."""
HARD_LIMIT = 4000
"""Refuse sequences longer than this (none of the scored windows comes close)."""
JUNCTION_HALF_WINDOW = 10
"""Residues on each side of the crossover in the junction window (2 x 10 residues)."""

ROLES = ("chimera", "parent5", "parent3")


@dataclass(frozen=True)
class Side:
    """One scored sequence of a record, in its protein's own 1-based coordinates."""

    sequence: str
    start: int
    """Protein position of ``sequence[0]``."""
    domain: tuple[int, int]
    junction: tuple[int, int]
    """Junction window (clipped to the scored sequence)."""
    paired: tuple[int, ...] = ()
    """Extra protein positions needed for the chimera-vs-parent pairing."""

    @property
    def end(self) -> int:
        return self.start + len(self.sequence) - 1

    def contains(self, position: int) -> bool:
        return self.start <= position <= self.end

    def seq_pos(self, position: int) -> int:
        return position - self.start + 1

    def span(self, bounds: tuple[int, int]) -> list[int]:
        return [p for p in range(bounds[0], bounds[1] + 1) if self.contains(p)]

    def needed(self) -> list[int]:
        """1-based sequence positions to mask and score."""
        wanted = set(self.span(self.domain)) | set(self.span(self.junction))
        wanted |= {p for p in self.paired if self.contains(p)}
        return sorted(self.seq_pos(p) for p in wanted)


@dataclass(frozen=True)
class RecordPlan:
    record: ChimeraRecord
    chimera: Side
    parent5: Side
    parent3: Side
    pairs5: tuple[tuple[int, int], ...]
    """(chimera position, 5' parent position) for 5'-derived chimera domain residues."""
    pairs3: tuple[tuple[int, int], ...]
    problem: str = ""
    """Non-empty: the record cannot be scored (reason); it is reported as not assessed."""

    def windowed_roles(self) -> list[str]:
        """Roles whose sequence is longer than the model context (scored in windows)."""
        return [role for role, side in self.sides().items() if len(side.sequence) > MAX_RESIDUES]

    def sides(self) -> dict[str, Side]:
        return {"chimera": self.chimera, "parent5": self.parent5, "parent3": self.parent3}


def _clip(lo: int, hi: int, start: int, end: int) -> tuple[int, int]:
    return max(lo, start), min(hi, end)


def plan_record(record: ChimeraRecord, w: int = JUNCTION_HALF_WINDOW) -> RecordPlan:
    c_start, c_end = record.window()
    f_start, f_end = record.five.window()
    t_start, t_end = record.three.window()
    j = record.junction_index
    j5 = record.five.junction_residue
    j3 = record.three.junction_residue
    pairs5: list[tuple[int, int]] = []
    pairs3: list[tuple[int, int]] = []
    for position in range(record.domain_start, record.domain_end + 1):
        origin, parent_position = record.origin(position)
        if origin == "five" and parent_position is not None and f_start <= parent_position <= f_end:
            pairs5.append((position, parent_position))
        elif (
            origin == "three"
            and parent_position is not None
            and t_start <= parent_position <= t_end
        ):
            pairs3.append((position, parent_position))
    chimera = Side(
        record.chimeric_domain_flanked,
        c_start,
        (record.domain_start, record.domain_end),
        _clip(j - w + 1, j + w, c_start, c_end),
    )
    parent5 = Side(
        record.five.protein[f_start - 1 : f_end],
        f_start,
        (record.five.domain_start, record.five.domain_end),
        _clip(j5 - w + 1, j5 + w, f_start, f_end),
        tuple(p for _, p in pairs5),
    )
    parent3 = Side(
        record.three.protein[t_start - 1 : t_end],
        t_start,
        (record.three.domain_start, record.three.domain_end),
        _clip(j3 - w, j3 + w - 1, t_start, t_end),
        tuple(p for _, p in pairs3),
    )
    problem = ""
    for role, side in (("chimera", chimera), ("parent5", parent5), ("parent3", parent3)):
        if not side.sequence:
            problem = f"{role}: empty sequence"
        elif not side.start <= side.domain[0] <= side.domain[1] <= side.end:
            problem = f"{role}: domain span is outside the scored sequence"
        elif len(side.sequence) > HARD_LIMIT:
            problem = f"{role}: {len(side.sequence)} residues exceeds limit {HARD_LIMIT}"
        elif set(side.sequence) - CANONICAL:
            bad = "".join(sorted(set(side.sequence) - CANONICAL))
            problem = f"{role}: non-canonical residue(s) {bad}"
        if problem:
            break
    if (
        not problem
        and chimera.sequence[
            chimera.seq_pos(record.domain_start) - 1 : chimera.seq_pos(record.domain_end)
        ]
        != record.chimeric_domain
    ):
        problem = "chimera: domain is not a slice of the scored window"
    return RecordPlan(record, chimera, parent5, parent3, tuple(pairs5), tuple(pairs3), problem)


def check_against_metrics(plan: RecordPlan, row: Mapping[str, Any] | None) -> str:
    """Cross-check a record against its PR #155 metrics row; '' or the mismatch."""
    if row is None:
        return "record not in the PR #155 metrics table"
    record = plan.record
    expected = {
        "sequence_id": record.sequence_id,
        "fusion": record.fusion,
        "chimeric_domain_span": f"{record.domain_start}-{record.domain_end}",
        "five_prime_domain": f"{record.five.domain_start}-{record.five.domain_end}",
        "three_prime_domain": f"{record.three.domain_start}-{record.three.domain_end}",
        "junction_index": str(record.junction_index),
        "chimeric_domain_length": str(len(record.chimeric_domain)),
    }
    for column, value in expected.items():
        got = str(row.get(column, "")).strip()
        if got.endswith(".0"):
            got = got[:-2]
        if got != value:
            return f"sequence mismatch vs PR #155 metrics: {column} {got!r} != {value!r}"
    return ""


class Scorer(Protocol):
    spec: ModelSpec

    def score(self, sequence: str, positions: list[int]) -> dict[int, float]: ...


class ScoringError(RuntimeError):
    pass


@dataclass
class WorkerScorer:
    """Runs :mod:`.worker` under ``python`` and exchanges JSON lines with it."""

    python: str
    model: str
    device: str = "mps"
    batch_size: int = 64
    env: dict[str, str] | None = None
    script: Path | None = None
    """Worker script; defaults to :mod:`.worker` (tests substitute a fake)."""
    spec: ModelSpec = field(init=False)
    header: dict[str, Any] = field(init=False)

    def __post_init__(self) -> None:
        worker = self.script or Path(__file__).with_name("worker.py")
        self._proc = subprocess.Popen(
            [self.python, "-I", str(worker), self.model, self.device, str(self.batch_size)],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            text=True,
            env=self.env,
        )
        line = self._readline()
        try:
            self.header = json.loads(line)["header"]
        except (json.JSONDecodeError, KeyError, TypeError) as error:
            raise ScoringError(f"PLM worker sent no valid header: {line[:200]!r}") from error
        self.spec = ModelSpec.from_header(self.header)
        self._next = 0

    def _readline(self) -> str:
        assert self._proc.stdout is not None
        line = self._proc.stdout.readline()
        if not line:
            raise ScoringError(f"PLM worker exited (status {self._proc.poll()})")
        return line

    def score(self, sequence: str, positions: list[int]) -> dict[int, float]:
        assert self._proc.stdin is not None
        self._next += 1
        rid = str(self._next)
        request = {"id": rid, "sequence": sequence, "positions": positions}
        try:
            self._proc.stdin.write(json.dumps(request) + "\n")
            self._proc.stdin.flush()
        except OSError as error:
            raise ScoringError("PLM worker input pipe closed") from error
        try:
            reply = json.loads(self._readline())
        except json.JSONDecodeError as error:
            raise ScoringError("PLM worker sent an unparsable reply") from error
        if not isinstance(reply, dict):
            raise ScoringError("PLM worker reply is not an object")
        if reply.get("id") != rid:
            raise ScoringError("PLM worker reply out of order")
        if "error" in reply:
            raise ScoringError(reply["error"])
        if reply.get("sequence") != sequence:
            raise ScoringError("PLM worker scored a different sequence")
        try:
            return {int(k): v for k, v in reply["logp"].items()}
        except (KeyError, TypeError, ValueError, AttributeError) as error:
            raise ScoringError("PLM worker reply has no valid logp map") from error

    def close(self) -> None:
        if self._proc.stdin:
            self._proc.stdin.close()
        # MPS teardown can outlast a minute even after scoring has finished.
        # Keep waiting for graceful EOF shutdown: aborting here loses the
        # driver's completion status or prevents cached analysis from writing.
        self._proc.wait()


def score_sequences(
    sequences: Iterable[tuple[str, list[int]]],
    scorer: Scorer,
    cache: ScoreCache,
    log: Callable[[str], None] = print,
) -> dict[str, str]:
    """Score the missing positions of each sequence and cache them; returns failures.

    A reply that is missing a requested position or holds a non-finite or
    positive value is a failure for that sequence: nothing from it is cached.
    """
    failures: dict[str, str] = {}
    for sequence, positions in sequences:
        have = cache.get(sequence)
        missing = sorted(set(positions) - set(have))
        if not missing:
            continue
        try:
            got = scorer.score(sequence, missing)
        except ScoringError as error:
            failures[sequence] = str(error)
            log(f"FAILED len={len(sequence)}: {error}")
            continue
        if set(got) != set(missing):
            failures[sequence] = "worker returned wrong score positions or length"
            log(f"FAILED len={len(sequence)}: {failures[sequence]}")
            continue
        bad = [p for p in missing if not valid_logp(got[p])]
        if bad:
            failures[sequence] = f"invalid or missing scores at {len(bad)} position(s)"
            log(f"FAILED len={len(sequence)}: {failures[sequence]}")
            continue
        cache.put(sequence, {p: float(got[p]) for p in missing})
    return failures
