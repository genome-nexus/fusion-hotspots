"""HEURISTIC structure flag: a screen-out, not a confirmation.

``fold_doubtful`` -- a doubt rule fired; ``fold_supported`` -- every input
was assessed and nothing fired; ``not_assessed`` -- anything needed is missing
(failed prediction, MSA error, missing model, parse or mapping failure, or
parents that are not themselves confidently predicted under these settings).

Asymmetry. For near-identical paralogue chimeras an AlphaFold-family model is
fed essentially the parents' MSA and will usually predict a confident,
parent-like fold whatever the junction does, so ``fold_supported`` is *weak*
evidence that the chimera folds. A failed or clearly worse fold is
informative, which is why the flag is used to screen records out.

Thresholds (documented heuristics, not calibrated values):

* pLDDT 70 -- AlphaFold's "confident" band boundary; < 50 is "very low"
  (Jumper et al. 2021; Tunyasuvunakool et al. 2021).
* TM-score 0.5 -- same-fold threshold (Xu & Zhang 2010, Bioinformatics 26:889).
* 10 / 15 pLDDT drops and +5 A inter-half PAE -- about one confidence band /
  a clear shift relative to the like-for-like parent baselines.
"""

from __future__ import annotations

from typing import Any

FOLD_SUPPORTED = "fold_supported"
FOLD_DOUBTFUL = "fold_doubtful"
NOT_ASSESSED = "not_assessed"

CONFIDENT_PLDDT = 70.0
VERY_LOW_PLDDT = 50.0
SAME_FOLD_TM = 0.5
DOMAIN_PLDDT_DROP = -10.0
JUNCTION_PLDDT_DROP = -15.0
INTERHALF_PAE_RISE = 5.0
INTERHALF_PAE_FLOOR = 10.0

RULES = {
    "D1_domain_plddt_low": f"chimera domain pLDDT < {CONFIDENT_PLDDT:g} (a parent is >= it)",
    "D2_domain_plddt_drop": f"chimera minus better parent domain pLDDT <= {DOMAIN_PLDDT_DROP:g}",
    "D3_junction_very_low": f"chimera junction-window min pLDDT < {VERY_LOW_PLDDT:g} while the "
    f"better parent's junction-window mean is >= {CONFIDENT_PLDDT:g}",
    "D4_junction_plddt_drop": f"chimera minus better parent junction-window mean pLDDT <= "
    f"{JUNCTION_PLDDT_DROP:g}",
    "D5_halves_not_packed": f"chimera inter-half PAE >= {INTERHALF_PAE_FLOOR:g} A and >= "
    f"{INTERHALF_PAE_RISE:g} A above the lower parent inter-half PAE",
    "D6_fold_differs": f"mapped TM-score (normalised by the parent domain) < {SAME_FOLD_TM:g} "
    "against both parent predictions",
    "N0_prediction_missing": "a chimera or parent prediction failed / is missing",
    "N1_parents_not_confident": f"neither parent's domain pLDDT is >= {CONFIDENT_PLDDT:g}: "
    "no confident reference fold to compare with",
    "N2_metric_missing": "not doubtful, but a metric needed for fold_supported is missing",
    "S1_all_pass": "every input assessed, no doubt rule fired",
}

REQUIRED_FOR_SUPPORTED = (
    "chimera_domain_plddt",
    "parent5_domain_plddt",
    "parent3_domain_plddt",
    "delta_domain_plddt_vs_better_parent",
    "delta_junction_plddt_mean_vs_better_parent",
    "chimera_junction_plddt_min",
    "chimera_interhalf_pae",
    "vs_5p_tm",
    "vs_3p_tm",
)


def _num(row: dict[str, Any], key: str) -> float | None:
    value = row.get(key)
    return float(value) if isinstance(value, (int, float)) else None


def assign_flag(row: dict[str, Any]) -> tuple[str, list[str], list[str]]:
    """(flag, fired rules, missing inputs) for one metrics row."""
    missing_predictions = [
        role
        for role in ("chimera", "parent5", "parent3")
        if row.get(f"{role}_prediction_status") != "ok"
    ]
    if missing_predictions:
        reasons = [f"{r}: {row.get(f'{r}_prediction_status')}" for r in missing_predictions]
        return NOT_ASSESSED, ["N0_prediction_missing"], reasons

    parent_domain = [_num(row, "parent5_domain_plddt"), _num(row, "parent3_domain_plddt")]
    known = [v for v in parent_domain if v is not None]
    if known and max(known) < CONFIDENT_PLDDT:
        return NOT_ASSESSED, ["N1_parents_not_confident"], []

    fired: list[str] = []
    chimera_domain = _num(row, "chimera_domain_plddt")
    if chimera_domain is not None and known and chimera_domain < CONFIDENT_PLDDT:
        fired.append("D1_domain_plddt_low")
    delta_domain = _num(row, "delta_domain_plddt_vs_better_parent")
    if delta_domain is not None and delta_domain <= DOMAIN_PLDDT_DROP:
        fired.append("D2_domain_plddt_drop")
    junction_min = _num(row, "chimera_junction_plddt_min")
    parent_junction = [
        v
        for v in (
            _num(row, "parent5_junction_plddt_mean"),
            _num(row, "parent3_junction_plddt_mean"),
        )
        if v is not None
    ]
    if (
        junction_min is not None
        and parent_junction
        and junction_min < VERY_LOW_PLDDT
        and max(parent_junction) >= CONFIDENT_PLDDT
    ):
        fired.append("D3_junction_very_low")
    delta_junction = _num(row, "delta_junction_plddt_mean_vs_better_parent")
    if delta_junction is not None and delta_junction <= JUNCTION_PLDDT_DROP:
        fired.append("D4_junction_plddt_drop")
    chimera_pae = _num(row, "chimera_interhalf_pae")
    parent_pae = [
        v
        for v in (_num(row, "parent5_interhalf_pae"), _num(row, "parent3_interhalf_pae"))
        if v is not None
    ]
    if (
        chimera_pae is not None
        and parent_pae
        and chimera_pae >= INTERHALF_PAE_FLOOR
        and chimera_pae - min(parent_pae) >= INTERHALF_PAE_RISE
    ):
        fired.append("D5_halves_not_packed")
    tms = [_num(row, "vs_5p_tm"), _num(row, "vs_3p_tm")]
    if all(v is not None for v in tms) and max(v for v in tms if v is not None) < SAME_FOLD_TM:
        fired.append("D6_fold_differs")
    if fired:
        return FOLD_DOUBTFUL, fired, []

    missing = [key for key in REQUIRED_FOR_SUPPORTED if _num(row, key) is None]
    if not parent_pae:
        missing.append("parent_interhalf_pae")
    for short in ("5p", "3p"):
        if row.get(f"vs_{short}_status") != "ok":
            missing.append(f"vs_{short}: {row.get(f'vs_{short}_status')}")
    if missing:
        return NOT_ASSESSED, ["N2_metric_missing"], missing
    return FOLD_SUPPORTED, ["S1_all_pass"], []
