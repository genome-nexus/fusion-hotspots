"""Driver stages: scoring order, resumable scoring, metrics rows and outputs."""

from __future__ import annotations

import csv
import gzip
import random
import time
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping

from cfh.domain_function.record import ChimeraRecord
from cfh.domain_plm.cache import ScoreCache
from cfh.domain_plm.metrics import METRIC_COLUMNS, not_assessed, record_metrics
from cfh.domain_plm.scoring import RecordPlan, Scorer, score_sequences

Log = Callable[[str], None]

NAMED_FUSIONS = (
    "FLI1::ERG",
    "STAT5B::STAT5A",
    "GPC6::GPC5",
    "GPC4::GPC3",
    "GPC3::GPC4",
    "LPO::MPO",
    "PIM3::ULK1",
)
NAMED_DNA_FUSIONS = ("MAP3K3::STRADA", "GNA14::GNAQ", "ABCA6::ABCA10", "AKR1C2::AKR1C3")
ORDER_SEED = 20261113
"""Seed of the random order in which the non-priority records are scored."""


def read_tsv(path: Path) -> list[dict[str, str]]:
    opener: Any = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def scoring_order(
    plans: list[RecordPlan], shortlist_ids: list[str], named: Iterable[str], extra_ids: set[str]
) -> list[tuple[str, RecordPlan]]:
    """(priority group, plan) in scoring order.

    1. the PR #155 shortlist (in shortlist rank order); 2. every record of a
    named fusion; 3. the extra (PCAWG DNA) records; 4. every other record in a
    seeded random order, so whatever is scored when time runs out is a random
    sample of the rest rather than a biased slice.
    """
    by_id = {plan.record.record_id: plan for plan in plans}
    named = set(named)
    out: list[tuple[str, RecordPlan]] = []
    seen: set[str] = set()

    def add(group: str, record_id: str) -> None:
        if record_id in by_id and record_id not in seen:
            seen.add(record_id)
            out.append((group, by_id[record_id]))

    for record_id in shortlist_ids:
        add("shortlist", record_id)
    for plan in sorted(plans, key=lambda p: (p.record.fusion, p.record.record_id)):
        if plan.record.fusion in named and plan.record.record_id not in extra_ids:
            add("named", plan.record.record_id)
    for record_id in sorted(extra_ids):
        add("pcawg_dna", record_id)
    rest = sorted(rid for rid in by_id if rid not in seen)
    random.Random(ORDER_SEED).shuffle(rest)
    for record_id in rest:
        add("random_order", record_id)
    return out


def score_in_order(
    ordered: list[tuple[str, RecordPlan]],
    scorer: Scorer,
    cache: ScoreCache,
    log: Log,
    deadline: float | None = None,
) -> tuple[int, dict[str, str]]:
    """Score records in order (cached per sequence after each one); stop at ``deadline``."""
    failures: dict[str, str] = {}
    done = 0
    started = time.time()
    for index, (group, plan) in enumerate(ordered, start=1):
        if deadline is not None and time.time() > deadline:
            log(f"deadline reached after {done} records; stopping (resume picks up here)")
            break
        if plan.problem:
            continue
        work = [(side.sequence, side.needed()) for side in plan.sides().values()]
        failures.update(score_sequences(work, scorer, cache, log))
        done += 1
        if index % 25 == 0 or group != "random_order":
            rate = (time.time() - started) / max(done, 1)
            log(f"[{index}/{len(ordered)}] {group} {plan.record.fusion}: {rate:.1f} s/record")
    return done, failures


def collect(
    plan: RecordPlan, cache: ScoreCache | None
) -> tuple[dict[str, dict[int, float]], dict[str, str]]:
    """Cached scores per role, restricted to the positions the plan needs."""
    scores: dict[str, dict[int, float]] = {}
    for role, side in plan.sides().items():
        cached = cache.get(side.sequence) if cache is not None else {}
        scores[role] = {p: cached[p] for p in side.needed() if p in cached}
    return scores, {}


def metrics_row(
    plan: RecordPlan,
    cache: ScoreCache | None,
    failures: Mapping[str, str],
    mismatch: str,
    model_problem: str,
) -> dict[str, Any]:
    """Every fail-safe rule is applied here, before any metric is computed."""
    if model_problem:
        return not_assessed(model_problem)
    if mismatch:
        return not_assessed(mismatch)
    if plan.problem:
        return not_assessed(plan.problem)
    scores, _ = collect(plan, cache)
    for role, side in plan.sides().items():
        if side.sequence in failures:
            return not_assessed(f"{role}: scoring failed ({failures[side.sequence]})")
        if len(scores[role]) < len(side.needed()):
            return not_assessed(f"{role}: not scored (coverage limit or stale cache)")
    return record_metrics(plan, scores, dict(failures))


__all__ = [
    "METRIC_COLUMNS",
    "NAMED_DNA_FUSIONS",
    "NAMED_FUSIONS",
    "collect",
    "metrics_row",
    "read_tsv",
    "score_in_order",
    "scoring_order",
]


def records_by_id(records: Iterable[ChimeraRecord]) -> dict[str, ChimeraRecord]:
    return {record.record_id: record for record in records}
