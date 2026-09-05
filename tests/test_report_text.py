"""Unit tests for the deterministic results.json -> PDF-report text templates.

Every assertion here checks literal, byte-for-byte sentence text produced
from a synthetic ``results.json``-shaped fixture, per the requirement that
the PDF report's abstract/results-summary text is templated (never
LLM-generated, never a fixed generic sentence) directly from the numeric
fields already present in a run's ``results.json``.
"""

from __future__ import annotations

from cfh.reporting.text import (
    format_percent,
    format_stat,
    render_abstract,
    render_results_summary,
    significance_clause,
)

BRAF_WITH_REFERENCE = {
    "gene_symbol": "BRAF",
    "study_id": "msk_impact_2017",
    "summary": {
        "total_fusions": 4,
        "in_frame_count": 3,
        "in_frame_percent": 75.0,
        "domain_accession": "PF07714",
        "kinase_retained_count": 2,
        "kinase_retained_percent": 50.0,
        "in_frame_kinase_retained_count": 2,
    },
    "reference": {
        "citation": "PMC5461196",
        "fusion_count": 33,
        "in_frame_percent": 100.0,
        "domain_retained_percent": 100.0,
    },
    "algorithm_results": [
        {
            "Algorithm": "frequency",
            "Summary": {"analyzed_event_count": 4, "unique_partner_gene_count": 3},
            "Tables": {
                "Partner_gene_counts": [
                    {"Partner_gene": "AGK", "Event_count": 2},
                    {"Partner_gene": "SND1", "Event_count": 1},
                    {"Partner_gene": "KIAA1549", "Event_count": 1},
                ]
            },
            "Warnings": [],
        },
        {
            "Algorithm": "domain_retention",
            "Summary": {
                "fisher_odds_ratio": 3.5,
                "fisher_p_value": 0.02,
                "permutation_empirical_p_value": 0.15,
                "observed_in_frame_retention_rate": 0.6667,
            },
            "Tables": {"frame_domain_contingency_table": [[2, 0], [1, 1]]},
            "Warnings": [],
        },
    ],
}

RET_WITHOUT_REFERENCE = {
    "gene_symbol": "RET",
    "study_id": "msk_impact_50k_2026",
    "summary": {
        "total_fusions": 10,
        "in_frame_count": 8,
        "in_frame_percent": 80.0,
        "domain_accession": "PF07714",
        "kinase_retained_count": 6,
        "kinase_retained_percent": 60.0,
        "in_frame_kinase_retained_count": 6,
    },
    "reference": None,
    "gene_track": {
        "protein_length": 1114,
        "domains": [],
        "exon_boundaries_aa": [
            {"exon_rank": 14, "start_aa": 950, "end_aa": 1000},
            {"exon_rank": 15, "start_aa": 1000, "end_aa": 1114},
        ],
    },
    "algorithm_results": [
        {
            "Algorithm": "frequency",
            "Summary": {"analyzed_event_count": 10, "unique_partner_gene_count": 2},
            "Tables": {
                "Partner_gene_counts": [
                    {"Partner_gene": "CCDC6", "Event_count": 7},
                    {"Partner_gene": "KIF5B", "Event_count": 3},
                ]
            },
            "Warnings": [],
        },
        {
            "Algorithm": "domain_retention",
            "Summary": {
                "fisher_odds_ratio": 4.0,
                "fisher_p_value": 0.004,
                "permutation_empirical_p_value": 0.01,
                "observed_in_frame_retention_rate": 0.75,
            },
            "Tables": {"frame_domain_contingency_table": [[6, 0], [2, 2]]},
            "Warnings": [],
        },
        {
            # RET has no disruption_required_domains configured: an
            # inapplicable/excluded algorithm, present in the run but with
            # no statistics computed.
            "Algorithm": "domain_disruption",
            "Summary": {
                "fisher_odds_ratio": None,
                "fisher_p_value": None,
                "permutation_empirical_p_value": None,
                "observed_in_frame_disruption_rate": None,
            },
            "Tables": {},
            "Warnings": [
                "RET has no disruption_required_domains configured; "
                "domain-disruption analysis was skipped."
            ],
        },
        {
            "Algorithm": "cutpoint_detection",
            "Summary": {
                "determinable": True,
                "reason": None,
                "n_events_analyzed": 10,
                "inferred_cutpoint_aa": 1063,
                "observed_statistic_neg_log10_p": 5.0,
                "observed_p_value": 1e-5,
                "observed_odds_ratio": None,
                "corrected_p_value": 0.0009,
                "known_domain_boundary_comparison": {
                    "nearest_known_domain_boundary_aa": 1000,
                    "distance_aa": 63,
                },
            },
            "Tables": {"cutpoint_scan": []},
            "Warnings": [],
        },
        {
            # confidence_stats was run with no group_field configured, so it
            # failed and carries no usable statistics: another
            # inapplicable/excluded case, this time an outright plugin
            # failure rather than a graceful gene-config skip.
            "Algorithm": "confidence_stats",
            "Summary": {"Runtime_seconds": 0.001},
            "Warnings": ["Algorithm failed: ValueError: params['group_field'] is required"],
        },
    ],
}


def test_format_helpers_render_unavailable_for_missing_or_non_finite_values():
    assert format_stat(None) == "unavailable"
    assert format_stat(float("nan")) == "unavailable"
    assert format_stat(float("inf")) == "unavailable"
    assert format_stat(0.010084033613445377) == "0.010084"
    assert format_percent(None) == "unavailable"
    assert format_percent(94.28571428571429) == "94.3%"
    assert significance_clause(None) is None
    assert significance_clause(0.049) == "statistically significant at alpha=0.05"
    assert significance_clause(0.05) == "not statistically significant at alpha=0.05"


def test_abstract_with_configured_benchmark_reference_matches_exact_text():
    abstract = render_abstract(BRAF_WITH_REFERENCE)

    assert abstract == (
        "This report analyzes BRAF gene fusions from the msk_impact_2017 study, covering "
        "4 fusion events. 3/4 fusions (75.0%) were in-frame. 2/4 fusions (50.0%) retained "
        "the PF07714 domain. Domain retention was tested with Fisher's exact test "
        "(p=0.02, statistically significant at alpha=0.05). Compared to the literature "
        "benchmark (PMC5461196: 100.0% in-frame, 100.0% domain-retained), this run "
        "observed 75.0% in-frame and 50.0% domain retention."
    )


def test_abstract_without_benchmark_reference_states_none_configured_explicitly():
    abstract = render_abstract(RET_WITHOUT_REFERENCE)

    assert abstract == (
        "This report analyzes RET gene fusions from the msk_impact_50k_2026 study, "
        "covering 10 fusion events. 8/10 fusions (80.0%) were in-frame. 6/10 fusions "
        "(60.0%) retained the PF07714 domain. Domain retention was tested with Fisher's "
        "exact test (p=0.004, statistically significant at alpha=0.05). No literature "
        "benchmark is configured for RET."
    )
    assert "No literature benchmark is configured for RET." in abstract


def test_results_summary_frequency_paragraph_braf():
    sections = {s["algorithm"]: s["paragraph"] for s in render_results_summary(BRAF_WITH_REFERENCE)}

    assert sections["frequency"] == (
        "Fusion-partner frequency was tabulated across 4 analyzed fusion events, "
        "identifying 3 distinct fusion partner genes. The most frequent partner was AGK, "
        "observed in 2 events."
    )


def test_results_summary_domain_retention_paragraph_braf():
    sections = {s["algorithm"]: s["paragraph"] for s in render_results_summary(BRAF_WITH_REFERENCE)}

    assert sections["domain_retention"] == (
        "PF07714 domain retention was tested with Fisher's exact test comparing "
        "in-frame fusions against all others; 2/3 (66.7%) of in-frame fusions retained "
        "the domain, p=0.02 (statistically significant at alpha=0.05). A "
        "breakpoint-position permutation test produced a corroborating empirical "
        "p-value of 0.15."
    )


def test_results_summary_ret_orders_canonically_and_covers_excluded_algorithms():
    sections = render_results_summary(RET_WITHOUT_REFERENCE)
    order = [s["algorithm"] for s in sections]

    assert order == [
        "frequency",
        "domain_retention",
        "domain_disruption",
        "cutpoint_detection",
        "confidence_stats",
    ]

    by_name = {s["algorithm"]: s["paragraph"] for s in sections}

    assert by_name["frequency"] == (
        "Fusion-partner frequency was tabulated across 10 analyzed fusion events, "
        "identifying 2 distinct fusion partner genes. The most frequent partner was "
        "CCDC6, observed in 7 events."
    )
    assert by_name["domain_retention"] == (
        "PF07714 domain retention was tested with Fisher's exact test comparing "
        "in-frame fusions against all others; 6/8 (75.0%) of in-frame fusions retained "
        "the domain, p=0.004 (statistically significant at alpha=0.05). A "
        "breakpoint-position permutation test produced a corroborating empirical "
        "p-value of 0.01."
    )
    # Inapplicable/excluded: no disruption_required_domains configured for RET.
    assert by_name["domain_disruption"] == (
        "Domain-disruption analysis was skipped: RET has no disruption_required_domains "
        "configured; domain-disruption analysis was skipped."
    )
    assert by_name["cutpoint_detection"] == (
        "Cutpoint detection scanned 10 mapped breakpoints for the protein position that "
        "best separates domain-retained from lost/disrupted fusions; the inferred "
        "cutpoint was position 1063 aa (exon 15; permutation-corrected p=0.0009, "
        "statistically significant at alpha=0.05). This is 63 aa from the nearest "
        "configured domain boundary at 1000 aa."
    )
    # Inapplicable/excluded: confidence_stats failed with no group_field configured.
    assert by_name["confidence_stats"] == (
        "Corroborating confidence statistics were not computed for this run: Algorithm "
        "failed: ValueError: params['group_field'] is required"
    )


def test_cutpoint_detection_paragraph_states_exon_boundary_when_position_is_ambiguous():
    """Regression for BRAF: aa 327 is exactly the shared coordinate between
    exon 7's end and exon 8's start (a clamped intronic-breakpoint estimate
    artifact), so the paragraph must name both flanking exons rather than
    fabricating a single one.
    """
    payload = {
        "gene_track": {
            "protein_length": 766,
            "domains": [],
            "exon_boundaries_aa": [
                {"exon_rank": 7, "start_aa": 287, "end_aa": 327},
                {"exon_rank": 8, "start_aa": 327, "end_aa": 380},
            ],
        },
        "algorithm_results": [
            {
                "Algorithm": "cutpoint_detection",
                "Summary": {
                    "determinable": True,
                    "reason": None,
                    "n_events_analyzed": 179,
                    "inferred_cutpoint_aa": 327,
                    "observed_statistic_neg_log10_p": 1.4,
                    "observed_p_value": 0.04,
                    "observed_odds_ratio": 2.0,
                    "corrected_p_value": 0.035964,
                    "known_domain_boundary_comparison": None,
                },
                "Tables": {"cutpoint_scan": []},
                "Warnings": [],
            }
        ],
    }

    paragraph = render_results_summary(payload)[0]["paragraph"]

    assert paragraph == (
        "Cutpoint detection scanned 179 mapped breakpoints for the protein position that "
        "best separates domain-retained from lost/disrupted fusions; the inferred "
        "cutpoint was position 327 aa (exon 7/8 boundary; permutation-corrected "
        "p=0.035964, statistically significant at alpha=0.05)."
    )


def test_cutpoint_detection_paragraph_states_exon_unavailable_without_gene_track():
    payload = {
        "algorithm_results": [
            {
                "Algorithm": "cutpoint_detection",
                "Summary": {
                    "determinable": True,
                    "reason": None,
                    "n_events_analyzed": 5,
                    "inferred_cutpoint_aa": 200,
                    "corrected_p_value": 0.2,
                    "known_domain_boundary_comparison": None,
                },
                "Tables": {"cutpoint_scan": []},
                "Warnings": [],
            }
        ],
    }

    paragraph = render_results_summary(payload)[0]["paragraph"]

    assert (
        "(exon position unavailable; permutation-corrected p=0.2, "
        "not statistically significant at alpha=0.05)." in paragraph
    )


def test_cutpoint_detection_not_determinable_states_reason_explicitly():
    payload = {
        "algorithm_results": [
            {
                "Algorithm": "cutpoint_detection",
                "Summary": {
                    "determinable": False,
                    "reason": "fewer than 2 distinct domain-retention statuses were observed",
                },
            }
        ]
    }

    paragraph = render_results_summary(payload)[0]["paragraph"]

    assert paragraph == (
        "Cutpoint detection could not determine a breakpoint boundary for this run: "
        "fewer than 2 distinct domain-retention statuses were observed"
    )


def test_domain_retention_unavailable_states_reason_explicitly_not_invented():
    payload = {
        "summary": {},
        "algorithm_results": [
            {
                "Algorithm": "domain_retention",
                "Summary": {"fisher_p_value": None},
                "Warnings": [
                    "Domain-retention statistics are unavailable because no mapped "
                    "in-frame protein-fusion record has a known domain state."
                ],
            }
        ],
    }

    paragraph = render_results_summary(payload)[0]["paragraph"]

    assert paragraph == (
        "Domain-retention statistics were not computed for this run: Domain-retention "
        "statistics are unavailable because no mapped in-frame protein-fusion record "
        "has a known domain state."
    )


def test_frequency_partner_tie_break_is_deterministic_alphabetical():
    payload = {
        "algorithm_results": [
            {
                "Algorithm": "frequency",
                "Summary": {"analyzed_event_count": 4, "unique_partner_gene_count": 2},
                "Tables": {
                    "Partner_gene_counts": [
                        {"Partner_gene": "ZZZ", "Event_count": 2},
                        {"Partner_gene": "AAA", "Event_count": 2},
                    ]
                },
                "Warnings": [],
            }
        ]
    }

    paragraph = render_results_summary(payload)[0]["paragraph"]

    assert "The most frequent partner was AAA, observed in 2 events." in paragraph


def test_composite_score_forward_compatibility_hook_uses_only_present_fields():
    """No composite_score algorithm exists yet on this branch, but if a run's
    results.json ever carries one, the report must still say something
    accurate about it -- derived only from that result's own Tables/Warnings,
    never inventing a schema that doesn't exist.
    """
    payload = {
        "algorithm_results": [
            {
                "Algorithm": "composite_score",
                "Summary": {"n_ranked": 3},
                "Tables": {
                    "ranked_events": [
                        {"event_id": "E1", "score": 0.9},
                        {"event_id": "E2", "score": 0.7},
                        {"event_id": "E3", "score": 0.5},
                    ]
                },
                "Warnings": [],
            }
        ]
    }

    section = render_results_summary(payload)[0]

    assert section["algorithm"] == "composite_score"
    assert section["heading"] == "Composite evidence score"
    assert section["paragraph"] == (
        "Composite evidence score produced results for this run (3 row(s) in "
        "ranked_events); see the accompanying table."
    )


def test_unknown_future_algorithm_with_warning_reports_the_warning_verbatim():
    payload = {
        "algorithm_results": [
            {
                "Algorithm": "joint_partner",
                "Summary": {},
                "Tables": {},
                "Warnings": [
                    "Algorithm failed: ValueError: joint_partner requires a GeneConfig "
                    "with gene_pair"
                ],
            }
        ]
    }

    section = render_results_summary(payload)[0]

    assert section["heading"] == "Joint partner"
    assert section["paragraph"] == (
        "Joint partner reported: Algorithm failed: ValueError: joint_partner requires a "
        "GeneConfig with gene_pair"
    )
