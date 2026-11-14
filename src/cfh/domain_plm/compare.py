"""A descriptive PLM rank label, and where PLM, PR #155 tier and PR #160 fold flag disagree.

There is no validated way to combine these measures for chimeric domains, so
no combined score and no functional cutoff is defined. To list disagreements,
the PLM side needs *some* categorical reading; it is purely a rank within this
run, not a calibrated threshold:

``plm_rank`` = ``bottom_decile`` when the record's ``delta_vs_better_parent``
or ``dip_excess_worst`` is in the lowest ``RANK_FRACTION`` of the PR #155
records scored in this run (ties included), otherwise ``not_bottom_decile``;
``not_assessed`` whenever the PLM metrics are not assessed. A ``not_assessed``
measure never takes part in a disagreement (it is listed as not comparable).
"""

from __future__ import annotations

import math
from typing import Any

RANK_FRACTION = 0.10
RANK_MEASURES = ("delta_vs_better_parent", "dip_excess_worst")

DISAGREEMENT_RULES = {
    "plm_low_but_tier_likely_intact": "PLM bottom decile; PR #155 tier likely_intact",
    "plm_not_low_but_tier_likely_broken": "PLM not bottom decile; PR #155 tier likely_broken",
    "plm_low_but_fold_supported": "PLM bottom decile; PR #160 fold_supported",
    "plm_not_low_but_fold_doubtful": "PLM not bottom decile; PR #160 fold_doubtful",
    "tier_likely_intact_but_fold_doubtful": "PR #155 likely_intact; PR #160 fold_doubtful",
    "tier_likely_broken_but_fold_supported": "PR #155 likely_broken; PR #160 fold_supported",
}


def _number(value: Any) -> float | None:
    if value in ("", None) or isinstance(value, bool):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def percentile(value: float, reference: list[float]) -> float:
    """Fraction of ``reference`` values <= ``value``."""
    return sum(1 for ref in reference if ref <= value) / len(reference)


def add_ranks(rows: list[dict[str, Any]], reference_source: str = "pr155") -> None:
    """Add ``rank_<measure>`` percentiles and ``plm_rank`` to every row in place.

    The reference distribution is the scored rows of ``reference_source``;
    other rows (PCAWG DNA records) are ranked against it but do not shape it.
    """
    reference = {
        measure: sorted(
            n
            for row in rows
            if row.get("source") == reference_source and row.get("plm_status") == "scored"
            for n in [_number(row.get(measure))]
            if n is not None
        )
        for measure in RANK_MEASURES
    }
    for row in rows:
        values = {m: _number(row.get(m)) for m in RANK_MEASURES}
        if (
            row.get("plm_status") != "scored"
            or any(v is None for v in values.values())
            or any(not reference[m] for m in RANK_MEASURES)
        ):
            for measure in RANK_MEASURES:
                row[f"rank_{measure}"] = ""
            row["plm_rank"] = "not_assessed"
            continue
        low = False
        for measure, value in values.items():
            assert value is not None
            pct = percentile(value, reference[measure])
            row[f"rank_{measure}"] = round(pct, 4)
            low = low or pct <= RANK_FRACTION
        row["plm_rank"] = "bottom_decile" if low else "not_bottom_decile"


def disagreements(row: dict[str, Any]) -> list[str]:
    """Disagreement labels for one row (see :data:`DISAGREEMENT_RULES`)."""
    plm = row.get("plm_rank", "not_assessed")
    tier = row.get("pr155_tier", "")
    fold = row.get("pr160_fold_flag", "")
    out = []
    if plm == "bottom_decile" and tier == "likely_intact":
        out.append("plm_low_but_tier_likely_intact")
    if plm == "not_bottom_decile" and tier == "likely_broken":
        out.append("plm_not_low_but_tier_likely_broken")
    if plm == "bottom_decile" and fold == "fold_supported":
        out.append("plm_low_but_fold_supported")
    if plm == "not_bottom_decile" and fold == "fold_doubtful":
        out.append("plm_not_low_but_fold_doubtful")
    if tier == "likely_intact" and fold == "fold_doubtful":
        out.append("tier_likely_intact_but_fold_doubtful")
    if tier == "likely_broken" and fold == "fold_supported":
        out.append("tier_likely_broken_but_fold_supported")
    return out
