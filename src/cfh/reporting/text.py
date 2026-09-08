"""Deterministic, ``results.json``-derived text templates for the PDF report.

Every sentence produced here is computed directly from numeric fields
already present in a run's ``results.json`` payload -- the same dict written
by ``cfh.real_benchmark.write_outputs`` and read back via ``json.load``.
Nothing here calls an LLM and nothing here is a fixed generic sentence
independent of the run's numbers: a needed value that is absent (no
``benchmark_reference`` configured, an algorithm that was skipped or
inapplicable, ...) causes the corresponding sentence/clause to be omitted or
replaced with an explicit "not configured"/"unavailable" statement, never an
invented value.

Every function here is a pure function of its ``payload``/``result`` dict
argument, so the produced text is byte-for-byte reproducible for a given
``results.json``.
"""

from __future__ import annotations

import math
from typing import Any

from cfh.reporting.domain_names import (
    configured_domain_names,
    domain_interpretation_sentence,
    domain_label_for_accession,
    format_domain_names,
)
from cfh.reporting.exon_labels import exon_label_for_protein_position

ALPHA = 0.05

# The order results-summary paragraphs are rendered in when the underlying
# algorithm is present in a run's ``algorithm_results``. Any algorithm run
# that isn't in this list (including a not-yet-existing one, e.g. a future
# ``composite_score``) still gets a paragraph -- via ``_generic_paragraph``
# unless a dedicated renderer is registered for it -- appended after the
# canonically-ordered ones, in the order it appears in ``algorithm_results``.
CANONICAL_ALGORITHM_ORDER = [
    "frequency",
    "domain_retention",
    "domain_disruption",
    "cutpoint_detection",
    "confidence_stats",
    "composite_score",
]

_HEADINGS = {
    "frequency": "Fusion partner frequency",
    "domain_retention": "Domain retention",
    "domain_disruption": "Domain disruption",
    "cutpoint_detection": "Cutpoint detection",
    "confidence_stats": "Corroborating confidence statistics",
    "composite_score": "Composite evidence score",
}


def _finite(value: Any) -> float | None:
    """Coerce ``value`` to a finite float, or ``None`` if that isn't possible."""
    if isinstance(value, bool):
        return None
    if isinstance(value, (int, float)) and math.isfinite(value):
        return float(value)
    return None


def format_stat(value: Any) -> str:
    """Render a numeric statistic (e.g. a p-value) with fixed precision."""
    finite = _finite(value)
    return "unavailable" if finite is None else f"{finite:.6g}"


def format_percent(value: Any, digits: int = 1) -> str:
    """Render a value already on a 0-100 percent scale, or 'unavailable'."""
    finite = _finite(value)
    return "unavailable" if finite is None else f"{finite:.{digits}f}%"


def significance_clause(p_value: Any, alpha: float = ALPHA) -> str | None:
    """Return the fixed 'statistically significant at alpha=X' clause, or None."""
    finite = _finite(p_value)
    if finite is None:
        return None
    verdict = "statistically significant" if finite < alpha else "not statistically significant"
    return f"{verdict} at alpha={alpha:g}"


def _configured_domains(payload: dict, field: str) -> list[dict]:
    summary = payload.get("summary") or {}
    domains = summary.get(field) or []
    return domains if isinstance(domains, list) else []


def _retention_domain_names(payload: dict) -> list[str]:
    domains = _configured_domains(payload, "configured_key_domains")
    if not domains:
        domains = _configured_domains(payload, "key_domains")
    names = configured_domain_names(domains)
    if names:
        return names
    summary = payload.get("summary") or {}
    label = domain_label_for_accession(
        (payload.get("gene_track") or {}).get("domains") or summary.get("key_domains") or [],
        summary.get("domain_accession"),
    )
    return [label] if label else []


def _find_algorithm(payload: dict, name: str) -> dict | None:
    for result in payload.get("algorithm_results") or []:
        if result.get("Algorithm") == name:
            return result
    return None


def _pretty_name(algorithm: str) -> str:
    return _HEADINGS.get(algorithm, algorithm.replace("_", " ").capitalize())


def render_abstract(payload: dict) -> str:
    """Render the 3-5 sentence ABSTRACT section for one run's ``results.json``."""
    gene = payload.get("gene_symbol") or "the gene"
    study = payload.get("study_id") or "the configured study"
    summary = payload.get("summary") or {}
    total = summary.get("total_fusions")

    sentences: list[str] = []

    if total is not None:
        sentences.append(
            f"This report analyzes {gene} gene fusions from the {study} study, covering "
            f"{total} fusion event{'' if total == 1 else 's'}."
        )
    else:
        sentences.append(f"This report analyzes {gene} gene fusions from the {study} study.")

    in_frame_count = summary.get("in_frame_count")
    if total and in_frame_count is not None:
        sentences.append(
            f"{in_frame_count}/{total} fusions ({format_percent(summary.get('in_frame_percent'))}) "
            "were in-frame."
        )

    domain_names = _retention_domain_names(payload)
    domain = format_domain_names(domain_names)
    retained_count = summary.get("kinase_retained_count")
    if total and domain and retained_count is not None:
        sentences.append(
            f"{retained_count}/{total} fusions "
            f"({format_percent(summary.get('kinase_retained_percent'))}) retained the "
            f"{domain}{'' if domain.lower().endswith('domain') else ' domain'}."
        )

    domain_retention = _find_algorithm(payload, "domain_retention")
    fisher_p = (
        (domain_retention.get("Summary") or {}).get("fisher_p_value") if domain_retention else None
    )
    if _finite(fisher_p) is not None:
        sig = significance_clause(fisher_p)
        sentence = (
            f"Domain retention was tested with Fisher's exact test (p={format_stat(fisher_p)}"
        )
        sentence += f", {sig})." if sig else ")."
        sentences.append(sentence)
    else:
        sentences.append("No domain-retention statistical test could be computed for this run.")

    reference = payload.get("reference")
    if reference:
        sentences.append(
            f"Compared to the literature benchmark ({reference.get('citation')}: "
            f"{format_percent(reference.get('in_frame_percent'))} in-frame, "
            f"{format_percent(reference.get('domain_retained_percent'))} domain-retained), "
            f"this run observed {format_percent(summary.get('in_frame_percent'))} in-frame and "
            f"{format_percent(summary.get('kinase_retained_percent'))} domain retention."
        )
    else:
        sentences.append(f"No literature benchmark is configured for {gene}.")

    return " ".join(sentences)


def _frequency_paragraph(result: dict, payload: dict) -> str:
    del payload
    summary = result.get("Summary") or {}
    analyzed = summary.get("analyzed_event_count")
    partners = summary.get("unique_partner_gene_count")
    sentences = [
        f"Fusion-partner frequency was tabulated across "
        f"{analyzed if analyzed is not None else 'an unknown number of'} analyzed fusion "
        f"event{'' if analyzed == 1 else 's'}, identifying "
        f"{partners if partners is not None else 'an unknown number of'} distinct fusion "
        f"partner gene{'' if partners == 1 else 's'}."
    ]
    counts = (result.get("Tables") or {}).get("Partner_gene_counts") or []
    if counts:
        top = sorted(
            counts, key=lambda row: (-(row.get("Event_count") or 0), str(row.get("Partner_gene")))
        )[0]
        top_count = top.get("Event_count")
        sentences.append(
            f"The most frequent partner was {top.get('Partner_gene')}, observed in "
            f"{top_count} event{'' if top_count == 1 else 's'}."
        )
    return " ".join(sentences)


def _domain_retention_paragraph(result: dict, payload: dict) -> str:
    summary = result.get("Summary") or {}
    top = payload.get("summary") or {}
    warnings = result.get("Warnings") or []
    fisher_p = summary.get("fisher_p_value")
    if _finite(fisher_p) is None:
        reason = warnings[0] if warnings else "insufficient mapped in-frame domain-status data"
        return f"Domain-retention statistics were not computed for this run: {reason}"

    configured_domains = _configured_domains(payload, "configured_key_domains")
    if not configured_domains:
        configured_domains = _configured_domains(payload, "key_domains")
    domain_names = _retention_domain_names(payload)
    domain = format_domain_names(domain_names)
    domain_subject = domain or "The configured domain"
    if domain and not domain.lower().endswith("domain"):
        domain_subject += " domain"
    in_frame_count = top.get("in_frame_count")
    in_frame_retained = top.get("in_frame_kinase_retained_count")
    sig = significance_clause(fisher_p)
    sentence = (
        f"{domain_subject} retention was tested with Fisher's exact test "
        "comparing in-frame "
        "fusions against all others; "
    )
    if in_frame_count and in_frame_retained is not None:
        retained_pct = 100 * in_frame_retained / in_frame_count
        sentence += (
            f"{in_frame_retained}/{in_frame_count} ({retained_pct:.1f}%) of in-frame fusions "
            f"retained the domain, p={format_stat(fisher_p)}"
        )
    else:
        sentence += f"p={format_stat(fisher_p)}"
    sentence += f" ({sig})." if sig else "."
    sentences = [sentence]

    interpretation = domain_interpretation_sentence(
        configured_domains,
        fisher_p_value=summary.get("fisher_p_value"),
        fisher_odds_ratio=summary.get("fisher_odds_ratio"),
        effect="retention",
    )
    if interpretation:
        sentences.append(interpretation)

    perm_p = summary.get("permutation_empirical_p_value")
    if _finite(perm_p) is not None:
        sentences.append(
            "A breakpoint-position permutation test produced a corroborating empirical "
            f"p-value of {format_stat(perm_p)}."
        )
    return " ".join(sentences)


def _domain_disruption_paragraph(result: dict, payload: dict) -> str:
    summary = result.get("Summary") or {}
    warnings = result.get("Warnings") or []
    fisher_p = summary.get("fisher_p_value")
    if _finite(fisher_p) is None:
        reason = (
            warnings[0] if warnings else "no disruption-required domain is configured for this gene"
        )
        return f"Domain-disruption analysis was skipped: {reason}"

    table = (result.get("Tables") or {}).get("frame_domain_contingency_table") or [[0, 0], [0, 0]]
    disrupted_in_frame = table[0][0]
    other_in_frame = table[1][0]
    in_frame_total = disrupted_in_frame + other_in_frame
    sig = significance_clause(fisher_p)
    domain_names = configured_domain_names(
        _configured_domains(payload, "configured_disruption_required_domains")
    )
    domain = format_domain_names(domain_names)
    subject = f"Disruption of the {domain}" if domain else "Domain disruption"
    sentence = (
        f"{subject} was tested with "
        "Fisher's exact test comparing in-frame fusions against all others; "
    )
    if in_frame_total:
        pct = 100 * disrupted_in_frame / in_frame_total
        sentence += (
            f"{disrupted_in_frame}/{in_frame_total} ({pct:.1f}%) of in-frame fusions disrupted "
            f"the domain, p={format_stat(fisher_p)}"
        )
    else:
        sentence += f"p={format_stat(fisher_p)}"
    sentence += f" ({sig})." if sig else "."
    sentences = [sentence]

    interpretation = domain_interpretation_sentence(
        _configured_domains(payload, "configured_disruption_required_domains"),
        fisher_p_value=summary.get("fisher_p_value"),
        fisher_odds_ratio=summary.get("fisher_odds_ratio"),
        effect="disruption",
    )
    if interpretation:
        sentences.append(interpretation)

    perm_p = summary.get("permutation_empirical_p_value")
    if _finite(perm_p) is not None:
        sentences.append(
            "A breakpoint-position permutation test produced a corroborating empirical "
            f"p-value of {format_stat(perm_p)}."
        )
    return " ".join(sentences)


def _cutpoint_detection_paragraph(result: dict, payload: dict) -> str:
    summary = result.get("Summary") or {}
    if not summary.get("determinable"):
        reason = summary.get("reason") or "insufficient data"
        return (
            f"Cutpoint detection could not determine a breakpoint boundary for this run: {reason}"
        )

    n = summary.get("n_events_analyzed")
    cutpoint = summary.get("inferred_cutpoint_aa")
    corrected_p = summary.get("corrected_p_value")
    exon_label = exon_label_for_protein_position(
        (payload.get("gene_track") or {}).get("exon_boundaries_aa"), cutpoint
    )
    sig = significance_clause(corrected_p)
    sentence = (
        f"Cutpoint detection scanned {n if n is not None else 'an unknown number of'} mapped "
        "breakpoints for the protein position that best separates domain-retained from "
        f"lost/disrupted fusions; the inferred cutpoint was position {cutpoint} aa "
        f"({exon_label}; permutation-corrected p={format_stat(corrected_p)}"
    )
    sentence += f", {sig})." if sig else ")."
    sentences = [sentence]

    boundary = summary.get("known_domain_boundary_comparison")
    if boundary:
        sentences.append(
            f"This is {boundary.get('distance_aa')} aa from the nearest configured domain "
            f"boundary at {boundary.get('nearest_known_domain_boundary_aa')} aa."
        )
    return " ".join(sentences)


def _confidence_stats_paragraph(result: dict, payload: dict) -> str:
    del payload
    summary = result.get("Summary") or {}
    warnings = result.get("Warnings") or []
    group_field = summary.get("group_field")
    if not group_field:
        reason = warnings[0] if warnings else "no comparison groups were configured for this run"
        return f"Corroborating confidence statistics were not computed for this run: {reason}"

    a_label = summary.get("group_a_label")
    b_label = summary.get("group_b_label")
    n_a = summary.get("n_group_a")
    n_b = summary.get("n_group_b")
    sentences = [
        f'Corroborating confidence statistics compared {group_field} groups "{a_label}" '
        f'(n={n_a}) and "{b_label}" (n={n_b}).'
    ]

    mle = summary.get("mle")
    if mle:
        groups = mle.get("groups") or {}
        group_texts = []
        for label in (a_label, b_label):
            stats = groups.get(str(label))
            if not stats:
                continue
            point = stats.get("point_estimate")
            ci_low = stats.get("ci_low")
            ci_high = stats.get("ci_high")
            point_pct = format_percent(point * 100 if _finite(point) is not None else None)
            ci_low_pct = format_percent(ci_low * 100 if _finite(ci_low) is not None else None)
            ci_high_pct = format_percent(ci_high * 100 if _finite(ci_high) is not None else None)
            group_texts.append(
                f"{label}: {stats.get('successes')}/{stats.get('n')} ({point_pct}, "
                f"95% CI {ci_low_pct}-{ci_high_pct})"
            )
        if group_texts:
            sentences.append(f"For {mle.get('outcome_field')}, " + "; ".join(group_texts) + ".")

    ttest = summary.get("ttest")
    if ttest:
        p_value = ttest.get("p_value")
        sig = significance_clause(p_value)
        sentence = (
            f"Welch's t-test comparing {ttest.get('numeric_field')} between groups gave means "
            f"of {format_stat(ttest.get('mean_a'))} vs {format_stat(ttest.get('mean_b'))}, "
            f"p={format_stat(p_value)}"
        )
        sentence += f" ({sig})." if sig else "."
        sentences.append(sentence)

    return " ".join(sentences)


def _generic_paragraph(result: dict, payload: dict) -> str:
    """Fallback for an algorithm with no dedicated renderer (e.g. a future
    ``composite_score`` before it ships a known schema, or any other
    forward/unknown registered algorithm).

    Only ever describes the shape of the result already present in
    ``results.json`` -- counts of rows in its own ``Tables`` entries, its own
    ``Warnings`` -- never anything gene- or algorithm-specific that isn't
    literally there.
    """
    del payload
    name = result.get("Algorithm") or "unknown_algorithm"
    warnings = result.get("Warnings") or []
    if warnings:
        return f"{_pretty_name(name)} reported: {warnings[0]}"

    tables = result.get("Tables") or {}
    row_counts = [
        (table_name, len(rows)) for table_name, rows in tables.items() if isinstance(rows, list)
    ]
    if row_counts:
        parts = ", ".join(f"{count} row(s) in {table_name}" for table_name, count in row_counts)
        return (
            f"{_pretty_name(name)} produced results for this run ({parts}); "
            "see the accompanying table."
        )
    return (
        f"{_pretty_name(name)} produced results for this run; "
        "see the accompanying data for details."
    )


_PARAGRAPH_BUILDERS = {
    "frequency": _frequency_paragraph,
    "domain_retention": _domain_retention_paragraph,
    "domain_disruption": _domain_disruption_paragraph,
    "cutpoint_detection": _cutpoint_detection_paragraph,
    "confidence_stats": _confidence_stats_paragraph,
}


def render_results_summary(payload: dict) -> list[dict[str, str]]:
    """Render one paragraph per algorithm present in ``payload["algorithm_results"]``.

    Returns an ordered list of ``{"algorithm", "heading", "paragraph"}``
    dicts: algorithms from :data:`CANONICAL_ALGORITHM_ORDER` come first (only
    if actually present in this run), followed by any other algorithm this
    run happens to carry, in the order ``algorithm_results`` lists them.
    """
    results_by_name: dict[str, dict] = {}
    encounter_order: list[str] = []
    for result in payload.get("algorithm_results") or []:
        name = result.get("Algorithm")
        if name is None:
            continue
        results_by_name[name] = result
        if name not in encounter_order:
            encounter_order.append(name)

    ordered_names = [name for name in CANONICAL_ALGORITHM_ORDER if name in results_by_name]
    ordered_names += [name for name in encounter_order if name not in ordered_names]

    sections = []
    for name in ordered_names:
        result = results_by_name[name]
        builder = _PARAGRAPH_BUILDERS.get(name, _generic_paragraph)
        sections.append(
            {
                "algorithm": name,
                "heading": _HEADINGS.get(name, name.replace("_", " ").capitalize()),
                "paragraph": builder(result, payload),
            }
        )
    return sections
