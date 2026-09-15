"""Unit tests for the deterministic summary.json -> manuscript text templates
(:mod:`cfh.reporting.manuscript_text`).

Every assertion checks literal, byte-for-byte sentence text produced from a
synthetic ``summary.json``-shaped fixture, per the same never-LLM-generated,
never-a-fixed-generic-sentence discipline already enforced for
``cfh.reporting.text.render_abstract`` in ``tests/test_report_text.py``.
"""

from __future__ import annotations

import json
import re

from cfh.reporting.manuscript_text import (
    render_discussion_bullets,
    render_gene_highlight,
    render_manhattan_caption,
    render_manuscript_abstract,
    render_manuscript_methods,
    render_manuscript_title,
)
from conftest import latest_run_dir

REAL_COHORT_SCAN_SUMMARY_JSON = (
    latest_run_dir("cohort-scan_msk_impact_50k_2026") / "cohort_scan" / "summary.json"
)

# ``genes_after_gating`` (2) intentionally equals ``len(genes)`` below (an
# uncapped run where every gated gene was attempted) -- these are genuinely
# different counts in general (see the module docstring / ``_scan_counts``),
# and ``test_manuscript_abstract_flags_a_max_genes_capped_run`` below covers
# the case where they differ.
PAYLOAD_WITH_SIGNIFICANT_AND_HONORABLE_MENTIONS = {
    "study_id": "msk_impact_50k_2026",
    "min_distinct_patients": 5,
    "total_genes_before_gating": 3919,
    "genes_after_gating": 2,
    "curated_gene_count": 1,
    "auto_config_gene_count": 1,
    "unresolved_gene_count": 0,
    "significant_genes": ["ETV6"],
    "significance_level": 0.05,
    "generated_at": "2026-09-04T14:36:53.295005+00:00",
    "algorithms_run": ["composite_score", "domain_disruption", "domain_retention"],
    "genes": [
        {
            "gene_symbol": "ETV6",
            "config_source": "auto",
            "status": "ok",
            "n_events_analyzed": 90,
            "in_frame_percent": 71.11111111111111,
            "domain_retention_percent": 75.55555555555556,
            "fisher_p_value": 5.123255667283157e-06,
            "fisher_odds_ratio": 4.2,
            "key_domains": [{"name": "Sterile alpha motif (SAM) domain"}],
            "min_fdr_adjusted_q_value": 0.004334274294521551,
            "fdr_significant": True,
            "top_composite_score": 0.41577810494717626,
        },
        {
            "gene_symbol": "RET",
            "config_source": "curated",
            "status": "ok",
            "n_events_analyzed": 194,
            "in_frame_percent": 75.25773195876289,
            "domain_retention_percent": 92.26804123711341,
            "fisher_p_value": 0.00041966557652448966,
            "fisher_odds_ratio": 3.1,
            "key_domains": [{"name": "Protein kinase domain"}],
            "min_fdr_adjusted_q_value": 0.11976161828924137,
            "fdr_significant": False,
            "top_composite_score": 0.3,
        },
    ],
    "honorable_mentions": [
        {
            "rank": 1,
            "gene_symbol": "RET",
            "fisher_p_value": 0.00041966557652448966,
            "min_fdr_adjusted_q_value": 0.11976161828924137,
            "n_events_analyzed": 194,
            "in_frame_percent": 75.25773195876289,
            "domain_retention_percent": 92.26804123711341,
            "note": (
                "Did not survive genome-wide multiple-testing correction (FDR-adjusted "
                "q-value at or above the significance threshold), but ranks highly by raw "
                "p-value among the non-FDR-significant genes and may warrant targeted "
                "follow-up. This is NOT a claim of statistical significance."
            ),
        }
    ],
}

PAYLOAD_WITH_NO_SIGNIFICANT_OR_HONORABLE_MENTIONS = {
    "study_id": "some_other_study",
    "min_distinct_patients": 5,
    "total_genes_before_gating": 10,
    "genes_after_gating": 1,
    "curated_gene_count": 0,
    "auto_config_gene_count": 1,
    "unresolved_gene_count": 0,
    "significant_genes": [],
    "significance_level": 0.05,
    "generated_at": None,
    "genes": [
        {
            "gene_symbol": "FAKE1",
            "config_source": "auto",
            "status": "ok",
            "n_events_analyzed": 3,
            "in_frame_percent": 66.6,
            "domain_retention_percent": None,
            "fisher_p_value": None,
            "min_fdr_adjusted_q_value": None,
            "fdr_significant": False,
            "top_composite_score": None,
        }
    ],
    "honorable_mentions": [],
}

# A run that never recorded FDR-significance/honorable-mentions data at all
# (e.g. an older/different-schema run artifact) -- distinct from the payload
# above, which explicitly recorded a confirmed-empty result.
PAYLOAD_MISSING_SIGNIFICANCE_FIELDS = {
    key: value
    for key, value in PAYLOAD_WITH_NO_SIGNIFICANT_OR_HONORABLE_MENTIONS.items()
    if key not in ("significant_genes", "honorable_mentions")
}

# A run where ``max_genes`` capped how many gated genes were actually
# attempted: 10 genes passed the recurrence gate, but only 3 were attempted.
PAYLOAD_CAPPED_BY_MAX_GENES = {
    "study_id": "capped_study",
    "min_distinct_patients": 5,
    "total_genes_before_gating": 50,
    "genes_after_gating": 10,
    "curated_gene_count": 0,
    "auto_config_gene_count": 3,
    "unresolved_gene_count": 0,
    "significant_genes": [],
    "significance_level": 0.05,
    "generated_at": None,
    "algorithms_run": ["domain_retention"],
    "genes": [
        {
            "gene_symbol": f"GENE{i}",
            "config_source": "auto",
            "status": "ok",
            "n_events_analyzed": 5,
            "in_frame_percent": 50.0,
            "domain_retention_percent": 50.0,
            "fisher_p_value": 0.5,
            "min_fdr_adjusted_q_value": 0.9,
            "fdr_significant": False,
            "top_composite_score": 0.1,
        }
        for i in range(3)
    ],
    "honorable_mentions": [],
}

# A gene that has a real Fisher p-value from ``domain_retention`` (so it
# entered the BH correction and shows up in the plottable Manhattan count)
# but no composite score at all (e.g. ``composite_score`` didn't run for
# it) -- used to prove the "q-value-bearing genes" count (Methods/
# Discussion) is never smaller than the "plottable" count (Results caption),
# which additionally requires a composite score.
_GENE_WITH_Q_BUT_NO_COMPOSITE = {
    "gene_symbol": "NOCOMPOSITE",
    "config_source": "auto",
    "status": "ok",
    "n_events_analyzed": 8,
    "in_frame_percent": 40.0,
    "domain_retention_percent": 40.0,
    "fisher_p_value": 0.2,
    "min_fdr_adjusted_q_value": 0.6,
    "fdr_significant": False,
    "top_composite_score": None,
}


def test_render_manuscript_title_uses_study_id():
    assert (
        render_manuscript_title(PAYLOAD_WITH_SIGNIFICANT_AND_HONORABLE_MENTIONS)
        == "Genome-wide fusion-hotspot analysis of msk_impact_50k_2026"
    )


def test_manuscript_abstract_with_significant_and_honorable_mentions_matches_exact_text():
    abstract = render_manuscript_abstract(PAYLOAD_WITH_SIGNIFICANT_AND_HONORABLE_MENTIONS)

    assert abstract == (
        "This manuscript synthesizes a genome-wide fusion-hotspot cohort scan of "
        "msk_impact_50k_2026: 3919 genes carried at least one structural-variant record in "
        "the cohort, of which 2 passed the >= 5-distinct-patient recurrence gate. All 2 gated "
        "genes were attempted: 1 using hand-curated gene configs and 1 auto-configured, with "
        "0 gated-in gene(s) left unresolvable; 2 of the 2 attempted genes were successfully "
        "analyzed with the full registered algorithm suite. 1 gene reached genome-wide "
        "Benjamini-Hochberg FDR significance (q < 0.05): ETV6 (90 events, 71.1% in-frame, "
        "75.6% domain-retained, Fisher p=5.12326e-06, q=0.00433427). 1 additional gene forms "
        "a highly ranked non-FDR-significant tier flagged for targeted follow-up (see "
        "Honorable mentions, below)."
    )


def test_manuscript_abstract_never_claims_a_gated_gene_was_analyzed_when_it_was_not():
    """Regression test: the abstract used to say ALL of ``genes_after_gating``
    "were analyzed with the full registered algorithm suite", conflating the
    recurrence-gate count with the actually-successfully-analyzed count (a
    gene can be gated in but left unresolvable, or resolved but still fail
    during analysis). The real committed msk_impact_50k_2026 run has exactly
    this shape: 544 genes passed the gate, but only 524 were successfully
    analyzed -- the other 20 had no resolvable canonical transcript in
    Genome Nexus (``config_source == "unresolved"``). The exact count (and
    which specific genes/failure mode account for the gap) shifts whenever
    the run is regenerated against live cBioPortal/Genome Nexus data; what
    must always hold is that the abstract states whatever the real analyzed
    count is, never the full gated count.
    """
    payload = json.loads(REAL_COHORT_SCAN_SUMMARY_JSON.read_text())
    abstract = render_manuscript_abstract(payload)
    methods = render_manuscript_methods(payload)

    n_scanned = len(payload["genes"])
    n_analyzed = sum(1 for row in payload["genes"] if row.get("status") == "ok")
    assert n_scanned == payload["genes_after_gating"] == 544
    assert n_analyzed == 524
    assert n_analyzed < n_scanned  # gap: 20 genes unresolved in Genome Nexus.

    for text in (abstract, methods):
        assert "524 of the 544 attempted genes were successfully analyzed" in text
        # Never claims the full 544 (or the 524 curated+auto) were analyzed.
        assert "544 gated genes were successfully analyzed" not in text
        assert "544 were analyzed" not in text


def test_manuscript_abstract_flags_a_max_genes_capped_run():
    """When a run was capped with ``max_genes``, ``genes_after_gating`` (the
    full recurrence-gate count) legitimately exceeds ``len(payload["genes"])``
    (what was actually attempted) -- the abstract must say so explicitly,
    never silently present the smaller attempted count as if it were the
    full gated count."""
    abstract = render_manuscript_abstract(PAYLOAD_CAPPED_BY_MAX_GENES)

    assert "3 of those gated genes were attempted in this run" in abstract
    assert "capped by this run's configured gene limit" in abstract
    # Never claims all 10 gated genes were attempted.
    assert "All 10 gated genes" not in abstract


def test_manuscript_abstract_states_no_significant_genes_and_no_mentions_explicitly():
    abstract = render_manuscript_abstract(PAYLOAD_WITH_NO_SIGNIFICANT_OR_HONORABLE_MENTIONS)

    assert (
        "No gene reached genome-wide Benjamini-Hochberg FDR significance (q < 0.05) in this "
        "scan." in abstract
    )
    assert "No honorable-mentions tier was produced for this scan." in abstract
    # Never invents a claim of significance or a mention count that isn't there.
    assert "0 additional" not in abstract


def test_manuscript_abstract_distinguishes_missing_significance_data_from_confirmed_zero():
    """Regression test: a payload where ``significant_genes``/
    ``honorable_mentions`` are entirely ABSENT (e.g. an older/different-schema
    run artifact) must not be rendered identically to a payload that
    explicitly recorded zero of each -- "we don't know" and "we know it's
    zero" are different facts."""
    abstract_missing = render_manuscript_abstract(PAYLOAD_MISSING_SIGNIFICANCE_FIELDS)
    abstract_confirmed_empty = render_manuscript_abstract(
        PAYLOAD_WITH_NO_SIGNIFICANT_OR_HONORABLE_MENTIONS
    )

    assert "FDR-significance data is unavailable for this run" in abstract_missing
    assert "Honorable-mentions tier data is unavailable for this run" in abstract_missing
    assert "No gene reached genome-wide Benjamini-Hochberg" not in abstract_missing
    assert "No honorable-mentions tier was produced" not in abstract_missing

    # The confirmed-empty payload keeps its existing, different wording.
    assert "No gene reached genome-wide Benjamini-Hochberg" in abstract_confirmed_empty
    assert "No honorable-mentions tier was produced" in abstract_confirmed_empty
    assert "unavailable for this run" not in abstract_confirmed_empty


def test_manuscript_abstract_never_contains_overclaiming_honorable_mention_language():
    abstract = render_manuscript_abstract(PAYLOAD_WITH_SIGNIFICANT_AND_HONORABLE_MENTIONS)
    assert "near-significant" not in abstract.lower()
    assert "near significant" not in abstract.lower()


def test_manuscript_methods_lists_algorithms_recorded_in_this_runs_own_data():
    """Regression test: Methods must list the algorithm suite recorded in
    THIS run's own ``payload["algorithms_run"]`` -- never a live import of
    ``cfh.algorithms.registry`` -- so it stays correct for a programmatic run
    that used a restricted algorithm subset, and re-rendering the SAME
    summary.json later (after the registry has since changed) still
    describes what this run actually did. Using a deliberately fake
    algorithm name here (not a real registered algorithm) proves the text
    is sourced from the payload, not the registry."""
    payload = {
        **PAYLOAD_WITH_SIGNIFICANT_AND_HONORABLE_MENTIONS,
        "algorithms_run": ["totally_fake_algorithm", "another_fake_one"],
    }
    methods = render_manuscript_methods(payload)

    assert "totally_fake_algorithm" in methods
    assert "another_fake_one" in methods
    assert "Benjamini-Hochberg" in methods
    assert "Fisher's exact test" in methods
    assert "permutation test" in methods
    assert "q < 0.05" in methods


def test_manuscript_methods_states_algorithms_not_recorded_when_field_absent():
    payload = {
        key: value
        for key, value in PAYLOAD_WITH_SIGNIFICANT_AND_HONORABLE_MENTIONS.items()
        if key != "algorithms_run"
    }
    methods = render_manuscript_methods(payload)
    assert "Which algorithms ran for this scan is not recorded" in methods


def test_manuscript_methods_uses_real_gating_and_config_counts():
    methods = render_manuscript_methods(PAYLOAD_WITH_SIGNIFICANT_AND_HONORABLE_MENTIONS)
    assert "at least 5 distinct patient(s)" in methods
    assert "2 of 3919 genes passed the gate" in methods
    assert "1 using hand-curated gene configs" in methods
    assert "1 auto-configured" in methods
    assert "2 of the 2 attempted genes were successfully analyzed" in methods


def test_manuscript_methods_p_value_correction_count_matches_q_bearing_genes_not_all_scanned():
    """Regression test for the real internal contradiction found in the
    committed msk_impact_50k_2026 report: Methods said all 544 scanned
    genes' p-values were corrected while the Results caption said only 359
    genes had a q-value -- Methods must instead cite the actual number of
    genes that produced a computable p-value (which is what the correction
    was really applied across), not the raw scanned-gene count."""
    payload = json.loads(REAL_COHORT_SCAN_SUMMARY_JSON.read_text())
    methods = render_manuscript_methods(payload)
    caption = render_manhattan_caption(payload)

    n_with_q = sum(1 for row in payload["genes"] if row.get("min_fdr_adjusted_q_value") is not None)
    assert n_with_q == 359
    assert f"across the {n_with_q} genes that produced at least one computable p-value" in methods
    # Never claims correction was applied across all 544 scanned genes.
    assert "across all 544" not in methods
    assert "across the 544" not in methods

    plottable_match = re.search(r"Genome-wide summary of (\d+) scanned gene", caption)
    assert plottable_match is not None
    plottable_count = int(plottable_match.group(1))
    # The Results caption's plottable count (q AND composite score both
    # present) can never exceed the Methods/Discussion's p-value-bearing
    # count (q present) -- it is a strict subset by construction.
    assert plottable_count <= n_with_q


def test_manhattan_caption_plottable_count_never_exceeds_q_bearing_gene_count():
    """Same self-consistency invariant as above, exercised on a synthetic
    payload where a gene has a real q-value but no composite score, so the
    two counts actually differ (not just coincidentally equal, as they are
    in the real committed run)."""
    payload = {
        **PAYLOAD_WITH_SIGNIFICANT_AND_HONORABLE_MENTIONS,
        "genes": [
            *PAYLOAD_WITH_SIGNIFICANT_AND_HONORABLE_MENTIONS["genes"],
            _GENE_WITH_Q_BUT_NO_COMPOSITE,
        ],
    }
    methods = render_manuscript_methods(payload)
    caption = render_manhattan_caption(payload)

    n_with_q = sum(1 for row in payload["genes"] if row.get("min_fdr_adjusted_q_value") is not None)
    assert n_with_q == 3  # ETV6, RET, NOCOMPOSITE

    plottable_match = re.search(r"Genome-wide summary of (\d+) scanned gene", caption)
    plottable_count = int(plottable_match.group(1))
    assert plottable_count == 2  # NOCOMPOSITE excluded: no composite score.

    assert f"across the {n_with_q} genes that produced at least one computable p-value" in methods
    assert plottable_count < n_with_q


def test_manhattan_caption_counts_only_plottable_points_above_threshold():
    caption = render_manhattan_caption(PAYLOAD_WITH_SIGNIFICANT_AND_HONORABLE_MENTIONS)

    # Both ETV6 and RET have a q-value and a composite score -> plottable.
    assert "2 scanned genes" in caption
    assert "q=0.05" in caption
    # Only ETV6 is FDR-significant among the plottable points.
    assert "1 gene above it" in caption


def test_manhattan_caption_handles_no_plottable_points():
    caption = render_manhattan_caption(PAYLOAD_WITH_NO_SIGNIFICANT_OR_HONORABLE_MENTIONS)
    assert "0 scanned genes" in caption
    assert "0 genes above it" in caption


def test_gene_highlight_significant_gene_matches_exact_text():
    row = PAYLOAD_WITH_SIGNIFICANT_AND_HONORABLE_MENTIONS["genes"][0]
    paragraph = render_gene_highlight(row, n_fdr_corrected_genes=2)

    assert paragraph == (
        "ETV6 was analyzed across 90 fusion events, 71.1% in-frame and 75.6% domain-retained. "
        "Domain-retention Fisher's exact test on ETV6's own fusion events alone gives "
        "p=5.12326e-06 (statistically significant at alpha=0.05). The Sterile alpha motif "
        "(SAM) domain appears to be required for retention. Among 2 genes scanned "
        "genome-wide in this cohort run, ETV6's Benjamini-Hochberg-adjusted q=0.00433427 "
        "(reaches genome-wide FDR significance)."
    )


def test_gene_highlight_honorable_mention_appends_note_verbatim():
    row = PAYLOAD_WITH_SIGNIFICANT_AND_HONORABLE_MENTIONS["genes"][1]
    mention = PAYLOAD_WITH_SIGNIFICANT_AND_HONORABLE_MENTIONS["honorable_mentions"][0]

    paragraph = render_gene_highlight(
        row, n_fdr_corrected_genes=2, honorable_mention_note=mention["note"]
    )

    assert paragraph.endswith(mention["note"])
    # RET's raw Fisher p-value (0.00042) is significant at alpha=0.05 even
    # though its FDR-adjusted q-value (0.12) is not -- exact text below
    # pins down that each verdict is stated in its own sentence, explicitly
    # labeled with the statistic (within-gene Fisher p vs. cohort-scale
    # BH-adjusted q) it describes, so neither verdict can be misread as
    # describing the other statistic or as a retraction of it.
    assert (
        "Domain-retention Fisher's exact test on RET's own fusion events alone gives "
        "p=0.000419666 (statistically significant at alpha=0.05). The Protein kinase domain "
        "appears to be required for retention. Among 2 genes scanned genome-wide in this "
        "cohort run, RET's Benjamini-Hochberg-adjusted q=0.119762 (does not reach genome-wide "
        "FDR significance)."
    ) in paragraph


def test_gene_highlight_never_labels_an_fdr_nonsignificant_q_value_as_significant():
    """Regression test: previously the raw p-value's significance clause was
    appended directly after the q-value with no statistic label, so a gene
    with a significant raw p-value but a non-significant FDR-adjusted
    q-value (RET: p=0.00042 < 0.05 <= q=0.12) rendered as
    '...q=0.119762 (statistically significant at alpha=0.05)' -- misreadable
    as claiming the *q-value* was significant, which is false and
    contradicts genome-wide FDR significance semantics.
    """
    row = PAYLOAD_WITH_SIGNIFICANT_AND_HONORABLE_MENTIONS["genes"][1]
    paragraph = render_gene_highlight(row, n_fdr_corrected_genes=2)

    assert "q=0.119762 (statistically significant" not in paragraph
    assert "q=0.119762 (does not reach genome-wide FDR significance)" in paragraph
    assert "p=0.000419666 (statistically significant at alpha=0.05)" in paragraph


def test_gene_highlight_p_and_q_sentences_read_as_answering_different_questions():
    """The Fisher p-value sentence answers a purely within-gene question --
    "is this domain-retention finding real for this gene" -- computed from
    RET's own fusion events alone, and carries no cohort-scale/cross-gene
    language. The FDR q-value sentence answers the separate, cohort-scale
    question "is this gene notable among every gene scanned genome-wide in
    this run", and says so explicitly (naming the actual gene count) rather
    than sharing an ambiguous pronoun/adjacency with the p-value sentence
    that a reader could misread as one sentence qualifying the other.
    """
    row = PAYLOAD_WITH_SIGNIFICANT_AND_HONORABLE_MENTIONS["genes"][1]
    paragraph = render_gene_highlight(row, n_fdr_corrected_genes=2)

    fisher_start = paragraph.index("Domain-retention Fisher's exact test")
    fisher_end = paragraph.index(". ", fisher_start) + 1
    fisher_sentence = paragraph[fisher_start:fisher_end]

    q_start = paragraph.index("Among 2 genes")
    q_sentence = paragraph[q_start:]

    # Cohort-scale scope language ("among the N genes scanned genome-wide")
    # appears only in the q-value sentence, never in the p-value sentence.
    assert "genes scanned genome-wide" in q_sentence
    assert "genes scanned genome-wide" not in fisher_sentence
    assert "among" not in fisher_sentence.lower()

    # Each sentence names its own subject explicitly instead of leaning on a
    # shared pronoun ("it"/"this") that could let the two verdicts blur.
    assert fisher_sentence.startswith("Domain-retention Fisher's exact test on RET's own")
    assert q_sentence.startswith("Among 2 genes scanned genome-wide in this cohort run, RET's")


def test_gene_highlight_states_missing_data_explicitly_not_invented():
    row = PAYLOAD_WITH_NO_SIGNIFICANT_OR_HONORABLE_MENTIONS["genes"][0]
    paragraph = render_gene_highlight(row)

    assert (
        "No domain-retention statistical test could be computed for FAKE1 in this scan."
        in paragraph
    )
    # Never invents a p-value/q-value that isn't there.
    assert "p=unavailable" not in paragraph


def test_gene_highlight_omits_fdr_verdict_when_fdr_significant_flag_is_absent():
    """Never invents an FDR-significance verdict for q from the q-value's
    own magnitude when the precomputed ``fdr_significant`` flag isn't in the
    row -- omission over invention, same discipline as everywhere else."""
    row = {
        **PAYLOAD_WITH_SIGNIFICANT_AND_HONORABLE_MENTIONS["genes"][1],
        "fdr_significant": None,
    }
    paragraph = render_gene_highlight(row, n_fdr_corrected_genes=2)

    assert "Among 2 genes scanned genome-wide in this cohort run, RET's " in paragraph
    assert "Benjamini-Hochberg-adjusted q=0.119762." in paragraph
    assert "reach genome-wide FDR significance" not in paragraph


def test_gene_highlight_omits_gene_count_when_not_provided():
    """When the caller doesn't have the full cohort payload's scanned-gene
    count handy, the FDR sentence still carries explicit cohort-scope
    language, but never invents a specific count -- it states the count is
    unknown rather than guessing (same omission-over-invention discipline
    as :func:`_count_display` elsewhere in this module)."""
    row = PAYLOAD_WITH_SIGNIFICANT_AND_HONORABLE_MENTIONS["genes"][1]
    paragraph = render_gene_highlight(row)

    assert "Among an unknown number of genes scanned genome-wide in this cohort run" in paragraph


def test_discussion_bullets_are_gated_on_real_fields():
    bullets = render_discussion_bullets(PAYLOAD_WITH_SIGNIFICANT_AND_HONORABLE_MENTIONS)

    joined = " ".join(bullets)
    assert "Benjamini-Hochberg FDR correction" in joined
    assert "conservative correction" in joined
    assert "2 genes that produced at least one computable p-value" in joined
    assert "1 gene(s) used hand-curated" in joined
    assert "1 gene(s) were auto-configured" in joined
    assert "kinase/catalytic-keyword heuristic" in joined
    assert "2026-09-04T14:36:53.295005+00:00" in joined
    assert "not validated clinical calls" in joined


def test_discussion_bullets_p_value_count_matches_methods_not_raw_scanned_count():
    """Same self-consistency requirement as Methods: the Discussion's FDR
    correction bullet must cite the same p-value-bearing gene count as
    Methods, not the raw number of scanned-gene rows."""
    payload = json.loads(REAL_COHORT_SCAN_SUMMARY_JSON.read_text())
    bullets = render_discussion_bullets(payload)
    methods = render_manuscript_methods(payload)

    joined = " ".join(bullets)
    assert "across the 359 genes that produced at least one computable p-value" in joined
    assert "across the 359 genes that produced at least one computable p-value" in methods


def test_discussion_bullets_omit_auto_config_caveat_when_no_auto_configured_genes():
    payload = {**PAYLOAD_WITH_SIGNIFICANT_AND_HONORABLE_MENTIONS, "auto_config_gene_count": 0}
    bullets = render_discussion_bullets(payload)
    assert not any("auto-configured" in bullet for bullet in bullets)


def test_discussion_bullets_omit_generated_at_caveat_when_absent():
    bullets = render_discussion_bullets(PAYLOAD_WITH_NO_SIGNIFICANT_OR_HONORABLE_MENTIONS)
    assert not any("retrieved live" in bullet for bullet in bullets)
    # The standing clinical-use caveat is always present.
    assert any("not validated clinical calls" in bullet for bullet in bullets)
