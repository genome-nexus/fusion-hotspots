"""Deterministic, ``summary.json``-derived text templates for the cross-gene
manuscript-style synthesis report (``paper.md``/``paper.pdf``).

Every sentence produced here is computed directly from numeric/list fields
already present in a cohort scan's ``summary.json`` payload -- the same dict
shape written by ``cfh.cohort.outputs.write_cohort_scan_outputs`` (the
``genes`` rows from ``build_summary_rows`` and the ``honorable_mentions``
list from ``build_honorable_mentions``), plus a small number of additive
top-level fields (``significance_level``, ``generated_at``, ``algorithms_run``).
Nothing here calls an LLM, imports a live/mutable registry, or is a fixed
generic sentence independent of the run's numbers: re-rendering the SAME
``summary.json`` must always produce byte-for-byte the same text, even after
the code that produced that ``summary.json`` has since changed. A needed
value that is absent causes the corresponding sentence/clause to be omitted
or replaced with an explicit "unavailable"/"unknown number of" statement,
never an invented value -- and a field that is missing entirely from
``payload`` (an older/different-schema run artifact) is never treated the
same as a field that is present with a confirmed empty/zero value, since
those are different facts ("we don't know" vs. "we know it's zero"). This
follows the same discipline already established by
:func:`cfh.reporting.text.render_abstract`, which this module extends
rather than reimplements (the shared numeric formatting helpers are reused
directly from there).

Counts that describe different pipeline stages are kept genuinely distinct
and are never conflated even when, in a common/uncapped run, they happen to
share the same value: the number of genes that passed the recurrence gate
(``genes_after_gating``, unaffected by a ``max_genes`` cap) is not the same
as the number of genes this run actually attempted (``len(payload["genes"])``,
which IS bounded by ``max_genes``), which is not the same as the number that
were successfully analyzed (``status == "ok"``, which can be smaller still --
a gene can have a resolved config yet fail during analysis), which is not
the same as the number that produced at least one computable p-value and
so entered the Benjamini-Hochberg correction (``min_fdr_adjusted_q_value``
is not ``None``). See :func:`_scan_counts`.

Every function here is a pure function of its arguments, so the produced
text is byte-for-byte reproducible for a given ``summary.json``.
"""

from __future__ import annotations

from cfh.reporting.domain_names import domain_interpretation_sentence
from cfh.reporting.text import (
    format_percent,
    format_stat,
    significance_clause,
)


def _plural(count: int | None, singular: str = "", plural: str = "s") -> str:
    return singular if count == 1 else plural


def _count_display(value: int | None) -> str:
    return "an unknown number of" if value is None else str(value)


def _level_display(significance_level: float | None) -> str:
    return "the configured" if significance_level is None else f"{significance_level:g}"


def _present_list(payload: dict, key: str) -> tuple[bool, list]:
    """Distinguish ``payload[key]`` being present (even as a confirmed empty
    list) from ``key`` being missing/``None`` entirely -- the latter means
    "this run's summary.json carries no data for this field", not "zero",
    and callers must render those two states differently (omission-over-
    invention: never say "no gene reached significance" when the truth is
    "we don't know, this field wasn't recorded")."""
    if key not in payload or payload.get(key) is None:
        return False, []
    return True, payload[key]


def _scan_counts(payload: dict) -> dict[str, int | bool]:
    """Derive every gene count this module reports from the SAME source --
    ``payload["genes"]``, one row per gene this run actually attempted --
    rather than from ``genes_after_gating`` or any other aggregate that can
    legitimately differ from it (see the module docstring). Returns:

    * ``n_scanned`` -- genes this run attempted (``len(payload["genes"])``).
    * ``n_analyzed`` -- of those, how many completed analysis successfully
      (``status == "ok"``); this can be smaller than the hand-curated +
      auto-configured count, since a gene with a resolved config can still
      fail during analysis.
    * ``n_with_q`` -- of those, how many produced at least one computable
      p-value and so actually entered the Benjamini-Hochberg correction
      (``min_fdr_adjusted_q_value is not None``); this is the true
      denominator for any statement about what the FDR correction was
      "applied across", and is never the same as ``n_scanned`` in a run
      where any gene had no computable p-value at all.
    * ``capped`` -- whether ``n_scanned`` is smaller than
      ``genes_after_gating`` (e.g. a ``max_genes``-capped run), meaning not
      every gene that passed the recurrence gate was actually attempted.
    """
    rows = payload.get("genes") or []
    n_scanned = len(rows)
    n_analyzed = sum(1 for row in rows if row.get("status") == "ok")
    n_with_q = sum(1 for row in rows if row.get("min_fdr_adjusted_q_value") is not None)
    genes_after_gating = payload.get("genes_after_gating")
    capped = genes_after_gating is not None and n_scanned < genes_after_gating
    return {
        "n_scanned": n_scanned,
        "n_analyzed": n_analyzed,
        "n_with_q": n_with_q,
        "capped": capped,
    }


def render_manuscript_title(payload: dict) -> str:
    """Render the manuscript's title from the scan's study id."""
    study = payload.get("study_id") or "the configured study"
    return f"Genome-wide fusion-hotspot analysis of {study}"


def render_manuscript_abstract(payload: dict) -> str:
    """Render the 3-7 sentence ABSTRACT section for one cohort scan's
    ``summary.json``-shaped payload."""
    study = payload.get("study_id") or "the configured study"
    total_before = payload.get("total_genes_before_gating")
    total_after = payload.get("genes_after_gating")
    min_patients = payload.get("min_distinct_patients")
    curated = payload.get("curated_gene_count")
    auto = payload.get("auto_config_gene_count")
    unresolved = payload.get("unresolved_gene_count")
    counts = _scan_counts(payload)

    sentences: list[str] = []

    if total_before is not None and total_after is not None and min_patients is not None:
        sentences.append(
            f"This manuscript synthesizes a genome-wide fusion-hotspot cohort scan of {study}: "
            f"{total_before} gene{_plural(total_before)} carried at least one structural-variant "
            f"record in the cohort, of which {total_after} passed the >= {min_patients}-distinct-"
            "patient recurrence gate."
        )
    else:
        sentences.append(
            f"This manuscript synthesizes a genome-wide fusion-hotspot cohort scan of {study}."
        )

    # Attempted/analyzed counts are derived from ``payload["genes"]`` itself
    # (see :func:`_scan_counts`), never from ``genes_after_gating`` -- that
    # gate count is unaffected by a ``max_genes`` cap and can legitimately
    # exceed the number of genes this run actually attempted, and even a
    # gene with a resolved (curated/auto) config can still fail during
    # analysis, so "curated + auto-configured" is not itself a safe stand-in
    # for "successfully analyzed" either.
    if counts["n_scanned"]:
        if counts["capped"]:
            lead = (
                f"{counts['n_scanned']} of those gated genes were attempted in this run "
                "(capped by this run's configured gene limit)"
            )
        else:
            lead = (
                f"All {counts['n_scanned']} gated gene{_plural(counts['n_scanned'])} were attempted"
            )
        sentences.append(
            f"{lead}: {_count_display(curated)} using hand-curated gene configs and "
            f"{_count_display(auto)} auto-configured, with {_count_display(unresolved)} "
            f"gated-in gene(s) left unresolvable; {counts['n_analyzed']} of the "
            f"{counts['n_scanned']} attempted gene{_plural(counts['n_scanned'])} were "
            "successfully analyzed with the full registered algorithm suite."
        )

    significance_level = payload.get("significance_level")
    level_display = _level_display(significance_level)
    has_significant_field, significant_genes = _present_list(payload, "significant_genes")
    rows_by_gene = {row["gene_symbol"]: row for row in payload.get("genes") or []}

    if significant_genes:

        def _q_sort_key(gene: str) -> float:
            q_value = rows_by_gene.get(gene, {}).get("min_fdr_adjusted_q_value")
            return q_value if q_value is not None else 1.0

        clauses = []
        for gene in sorted(significant_genes, key=_q_sort_key):
            row = rows_by_gene.get(gene, {})
            n_events = row.get("n_events_analyzed")
            clauses.append(
                f"{gene} ({_count_display(n_events)} event{_plural(n_events)}, "
                f"{format_percent(row.get('in_frame_percent'))} in-frame, "
                f"{format_percent(row.get('domain_retention_percent'))} domain-retained, "
                f"Fisher p={format_stat(row.get('fisher_p_value'))}, "
                f"q={format_stat(row.get('min_fdr_adjusted_q_value'))})"
            )
        sentences.append(
            f"{len(significant_genes)} gene{_plural(len(significant_genes))} reached genome-wide "
            f"Benjamini-Hochberg FDR significance (q < {level_display}): "
            + "; ".join(clauses)
            + "."
        )
    elif has_significant_field:
        sentences.append(
            "No gene reached genome-wide Benjamini-Hochberg FDR significance "
            f"(q < {level_display}) in this scan."
        )
    else:
        sentences.append(
            "FDR-significance data is unavailable for this run (this run's summary.json "
            "carries no `significant_genes` field)."
        )

    has_honorable_mentions_field, honorable_mentions = _present_list(payload, "honorable_mentions")
    if honorable_mentions:
        verb = "forms" if len(honorable_mentions) == 1 else "form"
        sentences.append(
            f"{len(honorable_mentions)} additional gene{_plural(len(honorable_mentions))} {verb} "
            "a highly ranked non-FDR-significant tier flagged for targeted follow-up (see "
            "Honorable mentions, below)."
        )
    elif has_honorable_mentions_field:
        sentences.append("No honorable-mentions tier was produced for this scan.")
    else:
        sentences.append(
            "Honorable-mentions tier data is unavailable for this run (this run's "
            "summary.json carries no `honorable_mentions` field)."
        )

    return " ".join(sentences)


def render_manuscript_methods(payload: dict) -> str:
    """Render the templated METHODS paragraph: data source, gating
    threshold, the algorithm suite actually recorded as having run for this
    scan (read from ``payload["algorithms_run"]`` -- data this specific run
    actually produced, never a live import of the currently-registered
    algorithm list, so re-rendering the SAME ``summary.json`` later, after
    the registry has changed, or for a programmatic run that used a
    restricted algorithm subset, still produces the correct text for THIS
    run), and the statistical tests applied."""
    study = payload.get("study_id") or "the configured study"
    min_patients = payload.get("min_distinct_patients")
    total_before = payload.get("total_genes_before_gating")
    total_after = payload.get("genes_after_gating")
    curated = payload.get("curated_gene_count")
    auto = payload.get("auto_config_gene_count")
    unresolved = payload.get("unresolved_gene_count")
    counts = _scan_counts(payload)

    sentences: list[str] = []

    if min_patients is not None and total_before is not None and total_after is not None:
        sentences.append(
            f"Structural-variant records were retrieved from the {study} cBioPortal study and "
            f"gated to genes with at least {min_patients} distinct patient(s) carrying a "
            f"structural-variant record ({total_after} of {total_before} genes passed the gate)."
        )
    else:
        sentences.append(
            f"Structural-variant records were retrieved from the {study} cBioPortal study."
        )

    # See the module docstring / :func:`_scan_counts`: "passed the gate"
    # (above), "attempted", "successfully analyzed", and "produced a
    # computable p-value" (below) are four genuinely different counts and
    # are never conflated, even though they coincide in the common,
    # uncapped, all-succeeded case.
    if counts["n_scanned"]:
        if counts["capped"]:
            lead = (
                f"{counts['n_scanned']} of those gated genes were attempted in this run "
                "(capped by this run's configured gene limit)"
            )
        else:
            lead = (
                f"All {counts['n_scanned']} gated gene{_plural(counts['n_scanned'])} were attempted"
            )
        sentences.append(
            f"{lead}: {_count_display(curated)} using hand-curated gene configs and "
            f"{_count_display(auto)} auto-configured from Genome Nexus canonical-transcript/"
            f"Pfam-domain data, with {_count_display(unresolved)} gated-in gene(s) left "
            f"unresolvable; {counts['n_analyzed']} of the {counts['n_scanned']} attempted "
            f"gene{_plural(counts['n_scanned'])} were successfully analyzed."
        )

    has_algorithms_field, algorithm_names = _present_list(payload, "algorithms_run")
    if algorithm_names:
        sentences.append(
            "Each successfully analyzed gene was run through the algorithm suite recorded "
            f"for this scan ({', '.join(algorithm_names)})."
        )
    elif has_algorithms_field:
        sentences.append("No algorithms are recorded as having run for this scan.")
    else:
        sentences.append(
            "Which algorithms ran for this scan is not recorded in this run's summary.json."
        )

    level_display = _level_display(payload.get("significance_level"))
    sentences.append(
        "Domain-retention and domain-disruption significance were assessed per gene with "
        "Fisher's exact test and a breakpoint-position permutation test; the resulting "
        f"p-values across the {counts['n_with_q']} gene{_plural(counts['n_with_q'])} that "
        "produced at least one computable p-value were jointly corrected with "
        f"Benjamini-Hochberg false-discovery-rate correction at q < {level_display}."
    )
    return " ".join(sentences)


def render_manhattan_caption(payload: dict) -> str:
    """Render the templated caption for the embedded genome-wide Manhattan
    figure: how many genes are plotted, the significance threshold, and how
    many genes sit above it -- computed with the same filter
    ``cfh.reporting.manhattan._plottable_points`` applies, so the caption
    never drifts from what the figure actually shows."""
    rows = payload.get("genes") or []
    plottable = [
        row
        for row in rows
        if row.get("min_fdr_adjusted_q_value") is not None
        and row.get("top_composite_score") is not None
    ]
    above = [row for row in plottable if row.get("fdr_significant")]
    level_display = _level_display(payload.get("significance_level"))
    return (
        f"Genome-wide summary of {len(plottable)} scanned gene{_plural(len(plottable))} with an "
        "FDR-adjusted q-value, ranked left-to-right by composite evidence score; the dashed line "
        f"marks the q={level_display} significance threshold, with {len(above)} "
        f"gene{_plural(len(above))} above it."
    )


def render_gene_highlight(
    row: dict,
    *,
    n_fdr_corrected_genes: int | None = None,
    honorable_mention_note: str | None = None,
) -> str:
    """Render the 2-5 sentence templated highlight paragraph for one
    highlighted gene's summary row: gene name, event count, in-frame%,
    domain-retention%, and Fisher/FDR-adjusted p-values.

    The Fisher p-value and the FDR-adjusted q-value answer two genuinely
    different questions, and are worded so a reader cannot mistake one for a
    retraction of the other:

    * the Fisher sentence answers "is this domain-retention finding real for
      *this gene*", computed from this gene's own fusion events alone, with
      no other gene involved;
    * the FDR sentence answers the separate, cohort-scale question "is this
      gene notable *relative to every other gene scanned genome-wide in this
      run*" -- a gene's own within-gene finding can be real even when it
      does not stand out in that larger comparison, so a non-significant
      q-value here is not evidence against the Fisher sentence above.

    ``n_fdr_corrected_genes`` is the number of genes that produced at least
    one computable p-value and so entered the Benjamini-Hochberg correction
    this q-value comes from (the same count reported in
    :func:`render_manuscript_methods`/:func:`render_discussion_bullets` --
    pass ``_scan_counts(payload)["n_with_q"]``); omitted (``None``) rather
    than guessed when the caller doesn't have the full cohort payload handy,
    per the omission-over-invention discipline used throughout this module.
    ``honorable_mention_note``, if given, is the gene's own already-
    established honorable-mentions note text (see
    :func:`cfh.cohort.outputs.build_honorable_mentions`), appended verbatim
    rather than re-derived.
    """
    gene = row.get("gene_symbol") or "This gene"
    n_events = row.get("n_events_analyzed")
    sentences: list[str] = []

    if n_events is not None:
        sentences.append(
            f"{gene} was analyzed across {n_events} fusion event{_plural(n_events)}, "
            f"{format_percent(row.get('in_frame_percent'))} in-frame and "
            f"{format_percent(row.get('domain_retention_percent'))} domain-retained."
        )
    else:
        sentences.append(f"Event counts are unavailable for {gene} in this scan.")

    fisher_p = row.get("fisher_p_value")
    fisher_display = format_stat(fisher_p)
    if fisher_display != "unavailable":
        sig = significance_clause(fisher_p)
        # Framed as a question about this gene alone -- "on {gene}'s own
        # fusion events alone" -- so it stands on its own and carries no
        # cohort-scale/cross-gene language that could be conflated with the
        # separate FDR sentence below.
        sentence = (
            f"Domain-retention Fisher's exact test on {gene}'s own fusion events alone "
            f"gives p={fisher_display}"
        )
        sentence += f" ({sig})." if sig else "."
        sentences.append(sentence)

        interpretation = domain_interpretation_sentence(
            row.get("key_domains") or [],
            fisher_p_value=fisher_p,
            fisher_odds_ratio=row.get("fisher_odds_ratio"),
            effect="retention",
        )
        if interpretation:
            sentences.append(interpretation)

        q_display = format_stat(row.get("min_fdr_adjusted_q_value"))
        if q_display != "unavailable":
            # ``fdr_significant`` is the same precomputed BH-correction verdict
            # used for this gene's badge/tier elsewhere (see
            # ``cfh.cohort.outputs._gene_badges``) -- reusing it here (rather
            # than re-deriving a threshold check on q) guarantees this clause
            # can never drift from, or be confused with, the Fisher sentence's
            # own significance clause above. This sentence deliberately opens
            # with explicit cohort-scale framing ("Among the N genes scanned
            # genome-wide...") -- a different grammatical subject than the
            # Fisher sentence above -- so the two verdicts can never be misread
            # as describing the same comparison: a gene's within-gene Fisher
            # p-value can be significant while its genome-wide BH-adjusted
            # q-value is not (that is the entire point of multiple-testing
            # correction), and neither statement retracts the other.
            fdr_significant = row.get("fdr_significant")
            if fdr_significant is True:
                q_clause = "reaches genome-wide FDR significance"
            elif fdr_significant is False:
                q_clause = "does not reach genome-wide FDR significance"
            else:
                q_clause = None
            gene_count_display = _count_display(n_fdr_corrected_genes)
            q_sentence = (
                f"Among {gene_count_display} gene{_plural(n_fdr_corrected_genes)} scanned "
                f"genome-wide in this cohort run, {gene}'s Benjamini-Hochberg-adjusted "
                f"q={q_display}"
            )
            q_sentence += f" ({q_clause})." if q_clause else "."
            sentences.append(q_sentence)
    else:
        sentences.append(
            f"No domain-retention statistical test could be computed for {gene} in this scan."
        )

    disruption_p = row.get("disruption_fisher_p_value")
    disruption_interpretation = domain_interpretation_sentence(
        row.get("disruption_required_domains") or [],
        fisher_p_value=disruption_p,
        fisher_odds_ratio=row.get("disruption_fisher_odds_ratio"),
        effect="disruption",
    )
    if disruption_interpretation:
        sentences.append(disruption_interpretation)

    if honorable_mention_note:
        sentences.append(honorable_mention_note)

    return " ".join(sentences)


def render_discussion_bullets(payload: dict) -> list[str]:
    """Render the templated INTERPRETATION/DISCUSSION caveats. Each caveat
    that cites a number is gated on that number actually being present in
    ``payload`` (never invented); the final caveat is a standing
    interpretive statement -- true regardless of any specific number in this
    run -- not a claim derived from a numeric field.
    """
    bullets: list[str] = []

    counts = _scan_counts(payload)
    n_with_q = counts["n_with_q"]
    level_display = _level_display(payload.get("significance_level"))
    bullets.append(
        "Cross-gene Benjamini-Hochberg FDR correction was applied jointly across the "
        f"{n_with_q} gene{_plural(n_with_q)} that produced at least one computable p-value "
        "in this scan. This reduces false-positive findings relative to testing each gene in "
        "isolation, but is a conservative correction: a real per-gene effect can fail to "
        f"reach the q < {level_display} threshold once corrected across that full "
        "p-value-bearing gene set."
    )

    curated = payload.get("curated_gene_count")
    auto = payload.get("auto_config_gene_count")
    if auto:
        bullets.append(
            f"{_count_display(curated)} gene(s) used hand-curated gene configurations, while "
            f"{auto} gene(s) were auto-configured from Genome Nexus canonical-transcript/Pfam-"
            "domain data using a kinase/catalytic-keyword heuristic to select the tracked domain; "
            "auto-configured domains have not been manually verified the way hand-curated ones "
            "have."
        )

    generated_at = payload.get("generated_at")
    study = payload.get("study_id") or "the configured study"
    if generated_at:
        bullets.append(
            f"Data were retrieved live from the public {study} cBioPortal study as of "
            f"{generated_at}. Some cBioPortal cohorts (e.g. actively accruing clinical-"
            "sequencing panels such as MSK-IMPACT) are updated periodically, so exact counts "
            "could shift if this scan is re-run later against such a cohort."
        )

    bullets.append(
        "All findings in this manuscript are computational, hypothesis-generating candidate "
        "evidence from bioinformatic analysis of public cohort data, not validated clinical "
        "calls; a gene's presence here is not a therapeutic or diagnostic recommendation."
    )
    return bullets
