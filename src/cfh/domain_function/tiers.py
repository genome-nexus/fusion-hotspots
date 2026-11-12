"""HEURISTIC rule-based tier for a chimeric domain.

No calibrated intact/broken cutoff exists for these metrics (calibration would
need assayed chimeras), so the tier is a transparent summary of the reported
metrics, not a prediction. Every rule is listed in :data:`RULES`; a record's
``tier_rules`` column names the rules that fired. Missing inputs never count in
favour of ``likely_intact``: a record that is not ``likely_broken`` and lacks
any required input is ``not_assessed``, and each missing input is named in
``tier_missing_inputs``.

Inputs required for ``likely_intact`` (:data:`REQUIRED_INPUTS`): the family HMM
scores for chimera and both parents; an assessed key-residue annotation for both
parents (UniProt entry present, transcript identity >= 90 %, no feature with
unusable coordinates near the domain) with every annotated residue resolved; and
a validated AlphaFold DB structure with DSSP for both parents
(``structure_status_* == "ok"``), even when the junction has no indel -- the
structure also feeds the reported contact metrics, and a missing or rejected
model must not read as a pass.

The SCHEMA-like broken-contact count E is reported but not used for tiering
(there is no universal E threshold); it appears only in disagreement flags.
"""

from __future__ import annotations

from typing import Any

from cfh.domain_function.metrics import KEY_ANNOTATION_OK, NOT_ASSESSED, STRUCTURE_OK

TIER_LABEL = "HEURISTIC"
LIKELY_INTACT = "likely_intact"
UNCERTAIN = "uncertain"
LIKELY_BROKEN = "likely_broken"
TIER_ORDER = (LIKELY_INTACT, UNCERTAIN, NOT_ASSESSED, LIKELY_BROKEN)

MIN_COVERAGE_BROKEN = 0.5
MIN_COVERAGE_INTACT = 0.8
MAX_BIT_DROP_INTACT = 10.0

RULES: dict[str, str] = {
    "B1_below_GA": (
        "likely_broken: the chimera's best domain hit scores below the family's Pfam domain "
        "GA while both parents' domains score at or above GA (or the chimera has no hit at all)."
    ),
    "B2_low_coverage": (
        f"likely_broken: the chimera's hits in the domain together span < "
        f"{MIN_COVERAGE_BROKEN:.0%} of the HMM's match states while both parents pass GA."
    ),
    "B3_point_feature_lost": (
        "likely_broken: a UniProt Active site / Binding site / Site residue of either parent "
        "is absent from the chimera and no chimera residue occupies its aligned column "
        "(junction-deleted column, or a column where the other parent has a gap)."
    ),
    "U1_parent_below_GA": (
        "uncertain: a parent's own domain scores below GA, so GA is not informative for this pair."
    ),
    "U2_bit_drop": (
        f"uncertain: chimera domain bits are more than {MAX_BIT_DROP_INTACT:g} below the weaker "
        "parent's."
    ),
    "U3_partial_coverage": (
        f"uncertain: chimera match-state coverage below {MIN_COVERAGE_INTACT:.0%}."
    ),
    "U4_point_feature_altered": (
        "uncertain: a point key residue differs from both parents' aligned residues."
    ),
    "U5_indel_in_ss": "uncertain: a junction indel residue is helix/strand in a parent model.",
    "U6_indel_near_key": (
        "uncertain: the junction has an indel and lies within 5 residues of a key feature."
    ),
    "U7_junction_outside_hit": (
        "uncertain: the chimera's best hit does not span the junction (it may be one parent's "
        "fragment)."
    ),
    "U8_domain_split": (
        "uncertain: HMMER splits the chimeric domain into more hits than either parent's domain."
    ),
    "U9_span_feature_lost_or_altered": (
        "uncertain: a UniProt DNA-binding / Motif / Region feature has a residue that is lost "
        "or altered (not merely swapped for the other parent's aligned residue)."
    ),
    "N1_missing": (
        "not_assessed: an input needed to call likely_intact is missing, rejected or "
        "incomplete (listed in tier_missing_inputs)."
    ),
    "I1_all_pass": (
        "likely_intact: none of the above fired and every input was assessed (passes GA, "
        f"coverage >= {MIN_COVERAGE_INTACT:.0%}, bit drop <= {MAX_BIT_DROP_INTACT:g}, no key "
        "residue or span feature lost/altered or unresolved, any indel outside parent "
        "helix/strand and away from key residues, both parent structures validated)."
    ),
}


REQUIRED_INPUTS = (
    "hmm (GA, parent and chimera hits)",
    "key_annotation_5p / key_annotation_3p in " + "/".join(KEY_ANNOTATION_OK),
    "no unresolved key residue (key_residue_mapping)",
    "indel secondary structure when the junction has an indel",
    f"structure_status_5p / structure_status_3p == {STRUCTURE_OK!r}",
)


def _missing(value: Any) -> bool:
    return value is None or value == NOT_ASSESSED or value == ""


def _number(value: Any) -> float | None:
    """A numeric metric, or None when it is missing / not assessed."""
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    return float(value)


def assign_tier(metrics: dict[str, Any]) -> tuple[str, list[str], list[str]]:
    """``(tier, fired_rules, missing_inputs)`` from one record's reported metrics."""
    fired: list[str] = []
    missing: list[str] = []
    ga = _number(metrics.get("ga_domain"))
    chimera_bits = _number(metrics.get("chimera_domain_bits"))
    five_bits = _number(metrics.get("parent5_domain_bits"))
    three_bits = _number(metrics.get("parent3_domain_bits"))
    coverage = _number(metrics.get("chimera_match_coverage"))
    union = _number(metrics.get("chimera_match_coverage_union"))
    broken_coverage = union if union is not None else coverage
    hmm_ok = ga is not None and five_bits is not None and three_bits is not None
    if not hmm_ok:
        missing.append("hmm")
    parents_pass = (
        ga is not None
        and five_bits is not None
        and three_bits is not None
        and five_bits >= ga
        and three_bits >= ga
    )

    if parents_pass and (chimera_bits is None or (ga is not None and chimera_bits < ga)):
        fired.append("B1_below_GA")
    if parents_pass and broken_coverage is not None and broken_coverage < MIN_COVERAGE_BROKEN:
        fired.append("B2_low_coverage")
    if metrics.get("key_point_lost"):
        fired.append("B3_point_feature_lost")
    if any(rule.startswith("B") for rule in fired):
        return LIKELY_BROKEN, fired, missing

    if hmm_ok and not parents_pass:
        fired.append("U1_parent_below_GA")
    if hmm_ok and chimera_bits is None:
        missing.append("chimera_hit")
    if five_bits is not None and three_bits is not None and chimera_bits is not None:
        if chimera_bits < min(five_bits, three_bits) - MAX_BIT_DROP_INTACT:
            fired.append("U2_bit_drop")
        if coverage is None:
            missing.append("coverage")
        elif coverage < MIN_COVERAGE_INTACT:
            fired.append("U3_partial_coverage")
        if metrics.get("chimera_hit_spans_junction") is False:
            fired.append("U7_junction_outside_hit")
    split = _number(metrics.get("chimera_hits_in_domain"))
    parent_hits = [_number(metrics.get(f"parent{p}_hits_in_domain")) for p in ("5", "3")]
    if split is not None and all(h is not None for h in parent_hits):
        if split > max(h for h in parent_hits if h is not None):
            fired.append("U8_domain_split")
    if metrics.get("key_point_altered"):
        fired.append("U4_point_feature_altered")
    if metrics.get("key_span_lost_or_altered"):
        fired.append("U9_span_feature_lost_or_altered")
    for prefix in ("5p", "3p"):
        status = metrics.get(f"key_annotation_status_{prefix}")
        if status not in KEY_ANNOTATION_OK:
            missing.append(f"key_annotation_{prefix}:{status or NOT_ASSESSED}")
    if metrics.get("key_features_assessed") is not True and not any(
        item.startswith("key_annotation_") for item in missing
    ):
        missing.append("key_residues")
    unresolved = metrics.get("key_unresolved")
    if isinstance(unresolved, int) and unresolved > 0:
        missing.append("key_residue_mapping")
    for prefix in ("5p", "3p"):
        status = metrics.get(f"structure_status_{prefix}")
        if status != STRUCTURE_OK:
            missing.append(f"structure_{prefix}:{status or NOT_ASSESSED}")
    indel_in_ss = metrics.get("indel_in_ss")
    if indel_in_ss is True:
        fired.append("U5_indel_in_ss")
    elif _missing(indel_in_ss):
        missing.append("indel_secondary_structure")
    has_indel = metrics.get("indel_ss") not in ("no_indel", None)
    if has_indel and metrics.get("junction_near_key") is True:
        fired.append("U6_indel_near_key")

    if fired:
        return UNCERTAIN, fired, missing
    if missing:
        return NOT_ASSESSED, ["N1_missing"], missing
    return LIKELY_INTACT, ["I1_all_pass"], missing


def disagreements(metrics: dict[str, Any], tier: str) -> list[str]:
    """Pairs of metrics pointing different ways; reported, never resolved."""
    flags = []
    ga = _number(metrics.get("ga_domain"))
    bits = _number(metrics.get("chimera_domain_bits"))
    passes = ga is not None and bits is not None and bits >= ga
    fails = (
        ga is not None
        and (bits is None or bits < ga)
        and _number(metrics.get("parent5_domain_bits")) is not None
    )
    if passes and metrics.get("key_point_lost"):
        flags.append("passes_GA_but_key_point_lost")
    if passes and metrics.get("indel_in_ss") is True:
        flags.append("passes_GA_but_indel_in_helix_or_strand")
    junction_inserts = metrics.get("junction_residues_in_inserts")
    if passes and isinstance(junction_inserts, int) and junction_inserts > 0:
        flags.append("passes_GA_with_junction_residues_in_insert_states")
    if (
        fails
        and metrics.get("key_features_assessed") is True
        and not (metrics.get("key_point_lost") or metrics.get("key_point_altered"))
    ):
        flags.append("below_GA_but_key_points_retained")
    core = [metrics.get(f"schema_E_core_{p}") for p in ("5p", "3p")]
    if tier == LIKELY_INTACT and any(isinstance(c, int) and c > 0 for c in core):
        flags.append("likely_intact_but_broken_core_contacts")
    energy = [metrics.get(f"schema_E_{p}") for p in ("5p", "3p")]
    if fails and all(e == 0 for e in energy):
        flags.append("below_GA_but_no_broken_contacts")
    return flags
