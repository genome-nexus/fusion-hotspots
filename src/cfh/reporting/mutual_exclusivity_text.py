"""Templated (never free-text/LLM-generated) report sentences for the
``mutation_cooccurrence`` algorithm's findings.

Follows the same report-generation discipline as
:func:`cfh.reporting.domain_names.domain_interpretation_sentence`: every
rendered claim is a deterministic string built from real numbers already
present in an ``AlgorithmResult``, and the absence of a computed finding is
stated explicitly rather than silently omitted.
"""

from __future__ import annotations

from typing import Any


def mutual_exclusivity_report_lines(
    result: dict[str, Any] | None, gene_symbol: str | None
) -> list[str]:
    """Render one markdown bullet per configured comparator target.

    ``result`` is a ``mutation_cooccurrence`` :class:`~cfh.model.algorithm_result.AlgorithmResult`
    (or its ``.model_dump()`` dict), or ``None`` if that algorithm was not
    part of this run at all (in which case nothing is rendered -- a run
    that never requested this algorithm has nothing to state an absence
    about). When the algorithm DID run but produced a no-op (unconfigured
    gene, or configured without the live-fetched comparator data), its own
    ``Warnings`` message -- already a real, specific statement of why --
    is rendered as the absence line, never a generic placeholder.
    """
    if result is None:
        return []
    algorithm = (
        result.get("Algorithm") if isinstance(result, dict) else getattr(result, "Algorithm", None)
    )
    if algorithm != "mutation_cooccurrence":
        return []

    summary = (result.get("Summary") if isinstance(result, dict) else result.Summary) or {}
    warnings = (result.get("Warnings") if isinstance(result, dict) else result.Warnings) or []
    targets = summary.get("targets") or []
    gene = gene_symbol or "This gene"

    if not targets:
        reason = warnings[0] if warnings else f"{gene} has no comparator alteration configured."
        return [f"- Mutation/CNA co-occurrence: not computed ({reason})"]

    lines = []
    for row in targets:
        p_value = row.get("p_value")
        odds_ratio = row.get("odds_ratio")
        direction = row.get("direction")
        table = row.get("contingency_table")
        label = row.get("label") or row.get("comparator_gene")
        direction_phrase = {
            "co_occurring": "tend to co-occur with",
            "mutually_exclusive": "tend toward mutual exclusivity with",
            "indeterminate": "show no directional tendency relative to",
        }.get(direction, "were compared against")
        p_value_is_numeric = isinstance(p_value, (int, float))
        significance = (
            "significantly" if p_value_is_numeric and p_value < 0.05 else "not significantly"
        )
        odds_ratio_text = (
            f"{odds_ratio:.3g}" if isinstance(odds_ratio, (int, float)) else "unavailable"
        )
        p_value_text = f"{p_value:.3g}" if p_value_is_numeric else "unavailable"
        lines.append(
            f"- {gene} fusions are {significance} associated with {label} across "
            f"{row.get('cohort_sample_count', 'the')} cohort samples and {direction_phrase} it "
            f"(2x2 table {table}; odds ratio {odds_ratio_text}, Fisher's exact two-sided "
            f"p={p_value_text})."
        )
    return lines
