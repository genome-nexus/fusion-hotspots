"""Length-normalised PLM comparisons of a chimeric domain with its two parents.

All values are masked pseudo-log-likelihoods: ``log p(residue | everything else,
residue masked)`` in nats, one per residue. Fusions with indels change the
domain length, so summed PLLs are not comparable across chimera and parents.
Every comparison here is length-robust:

* ``*_domain_pll`` -- **mean** per-residue PLL over the Pfam domain;
  ``delta_vs_5p`` / ``delta_vs_3p`` are chimera minus parent means;
* ``context_shift_5p`` / ``_3p`` -- **paired** mean difference over the same
  residues (chimera residue copied from the parent, scored in each context):
  inserted/hybrid residues have no partner and are excluded, deleted residues
  are absent, so length drops out entirely;
* junction window -- 2 x ``JUNCTION_HALF_WINDOW`` residues around the
  crossover. ``*_junction_dip`` is the window mean minus the same sequence's
  domain mean (self-normalised); ``dip_excess_vs_*`` subtracts the parent's
  dip at its own crossover residues, so a junction that falls in a naturally
  low-scoring loop is not counted against the chimera.

Fail closed: if any required score is missing or invalid the record is ``not
assessed`` with a reason and every metric is left empty.
"""

from __future__ import annotations

import math
from typing import Any

from cfh.domain_plm.cache import valid_logp
from cfh.domain_plm.scoring import RecordPlan, Side

METRIC_COLUMNS = (
    "plm_status",
    "plm_not_assessed_reason",
    "chimera_domain_pll",
    "parent5_domain_pll",
    "parent3_domain_pll",
    "delta_vs_5p",
    "delta_vs_3p",
    "delta_vs_better_parent",
    "context_shift_5p",
    "n_paired_5p",
    "context_shift_3p",
    "n_paired_3p",
    "chimera_junction_window",
    "chimera_junction_pll",
    "chimera_junction_min",
    "chimera_junction_dip",
    "parent5_junction_dip",
    "parent3_junction_dip",
    "dip_excess_vs_5p",
    "dip_excess_vs_3p",
    "dip_excess_worst",
)


def _mean(values: list[float]) -> float:
    return sum(values) / len(values)


class _Missing(Exception):
    pass


def _values(side: Side, scores: dict[int, float], positions: list[int]) -> list[float]:
    out = []
    for position in positions:
        value = scores.get(side.seq_pos(position))
        if value is None or not valid_logp(value):
            raise _Missing
        out.append(value)
    if not out:
        raise _Missing
    return out


def not_assessed(reason: str) -> dict[str, Any]:
    row: dict[str, Any] = {column: "" for column in METRIC_COLUMNS}
    row["plm_status"] = "not_assessed"
    row["plm_not_assessed_reason"] = reason
    return row


def record_metrics(
    plan: RecordPlan, scores: dict[str, dict[int, float]], failures: dict[str, str]
) -> dict[str, Any]:
    """Metrics for one record; ``scores`` maps role -> sequence position -> logp."""
    if plan.problem:
        return not_assessed(plan.problem)
    sides = plan.sides()
    for role, side in sides.items():
        if side.sequence in failures:
            return not_assessed(f"{role}: scoring failed ({failures[side.sequence]})")
    try:
        means: dict[str, float] = {}
        dips: dict[str, float] = {}
        junction: dict[str, list[float]] = {}
        for role, side in sides.items():
            domain = _values(side, scores.get(role, {}), side.span(side.domain))
            window = _values(side, scores.get(role, {}), side.span(side.junction))
            means[role] = _mean(domain)
            junction[role] = window
            dips[role] = _mean(window) - means[role]
        shifts: dict[str, float | None] = {}
        counts: dict[str, int] = {}
        for role, pairs in (("parent5", plan.pairs5), ("parent3", plan.pairs3)):
            counts[role] = len(pairs)
            shifts[role] = None
            if pairs:
                chim = _values(plan.chimera, scores.get("chimera", {}), [c for c, _ in pairs])
                par = _values(sides[role], scores.get(role, {}), [p for _, p in pairs])
                shifts[role] = _mean([a - b for a, b in zip(chim, par)])
    except _Missing:
        return not_assessed("missing or invalid PLM score for a required position")
    row: dict[str, Any] = {
        "plm_status": "scored",
        "plm_not_assessed_reason": "",
        "chimera_domain_pll": means["chimera"],
        "parent5_domain_pll": means["parent5"],
        "parent3_domain_pll": means["parent3"],
        "delta_vs_5p": means["chimera"] - means["parent5"],
        "delta_vs_3p": means["chimera"] - means["parent3"],
        "delta_vs_better_parent": means["chimera"] - max(means["parent5"], means["parent3"]),
        "context_shift_5p": "" if shifts["parent5"] is None else shifts["parent5"],
        "n_paired_5p": counts["parent5"],
        "context_shift_3p": "" if shifts["parent3"] is None else shifts["parent3"],
        "n_paired_3p": counts["parent3"],
        "chimera_junction_window": f"{plan.chimera.junction[0]}-{plan.chimera.junction[1]}",
        "chimera_junction_pll": _mean(junction["chimera"]),
        "chimera_junction_min": min(junction["chimera"]),
        "chimera_junction_dip": dips["chimera"],
        "parent5_junction_dip": dips["parent5"],
        "parent3_junction_dip": dips["parent3"],
        "dip_excess_vs_5p": dips["chimera"] - dips["parent5"],
        "dip_excess_vs_3p": dips["chimera"] - dips["parent3"],
        "dip_excess_worst": min(dips["chimera"] - dips["parent5"], dips["chimera"] - dips["parent3"]),
    }
    for key, value in row.items():
        if isinstance(value, float) and not math.isfinite(value):
            return not_assessed(f"non-finite metric {key}")
    return row


def junction_profile(plan: RecordPlan, scores: dict[int, float]) -> list[tuple[int, str, float]]:
    """(chimera position, residue, logp) over the junction window, for the report."""
    side = plan.chimera
    out = []
    for position in side.span(side.junction):
        value = scores.get(side.seq_pos(position))
        if value is None or not valid_logp(value):
            return []
        out.append((position, side.sequence[side.seq_pos(position) - 1], value))
    return out
