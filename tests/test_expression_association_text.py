from cfh.reporting.expression_association_text import expression_association_sentence


def test_none_when_summary_is_empty():
    assert expression_association_sentence("BRAF", {}) is None
    assert expression_association_sentence("BRAF", None) is None


def test_none_when_neither_comparison_produced_a_p_value():
    assert expression_association_sentence("BRAF", {"expression_field": "mRNA expression"}) is None


def test_fusion_positive_vs_negative_sentence_states_direction_and_significance():
    summary = {
        "fusion_positive_vs_negative": {
            "n_fusion_positive": 15,
            "n_fusion_negative": 483,
            "test": "mann_whitney_u",
            "statistic": 12.0,
            "p_value": 0.001,
            "mean_a": 2.5,
            "mean_b": 0.1,
        }
    }
    sentence = expression_association_sentence("BRAF", summary)
    assert (
        "BRAF mRNA expression is significantly higher in fusion-positive samples (n=15)" in sentence
    )
    assert "fusion-negative samples (n=483)" in sentence
    assert "mann_whitney_u, p=0.001" in sentence


def test_not_significant_result_is_stated_as_such_not_omitted():
    summary = {
        "fusion_positive_vs_negative": {
            "n_fusion_positive": 3,
            "n_fusion_negative": 20,
            "test": "welch_t_test",
            "statistic": 0.4,
            "p_value": 0.7,
            "mean_a": 1.0,
            "mean_b": 0.9,
        }
    }
    sentence = expression_association_sentence("RET", summary)
    assert "not significantly higher" in sentence


def test_domain_retention_split_sentence():
    summary = {
        "domain_retention_split": {
            "group_a_label": "retained",
            "group_b_label": "not_retained",
            "n_a": 9,
            "n_b": 6,
            "test": "mann_whitney_u",
            "statistic": 5.0,
            "p_value": 0.02,
            "mean_a": 3.0,
            "mean_b": -0.5,
        }
    }
    sentence = expression_association_sentence("BRAF", summary)
    assert "Among fusion-positive samples, BRAF expression" in sentence
    assert "kinase-domain-retained fusions (n=9)" in sentence
    assert "not-retained fusions (n=6)" in sentence


def test_both_comparisons_render_as_two_sentences():
    summary = {
        "fusion_positive_vs_negative": {
            "n_fusion_positive": 15,
            "n_fusion_negative": 483,
            "test": "mann_whitney_u",
            "statistic": 1.0,
            "p_value": 0.01,
            "mean_a": 2.0,
            "mean_b": 0.0,
        },
        "domain_retention_split": {
            "group_a_label": "retained",
            "group_b_label": "not_retained",
            "n_a": 9,
            "n_b": 6,
            "test": "mann_whitney_u",
            "statistic": 1.0,
            "p_value": 0.5,
            "mean_a": 2.1,
            "mean_b": 1.9,
        },
    }
    sentence = expression_association_sentence("BRAF", summary)
    assert "fusion-positive samples (n=15)" in sentence
    assert "kinase-domain-retained fusions (n=9)" in sentence
    # Both comparisons rendered as separate sentences, joined by ". ".
    assert sentence.count(". ") == 1
