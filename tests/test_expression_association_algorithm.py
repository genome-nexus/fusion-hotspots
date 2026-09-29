from cfh.algorithms import registry
from cfh.algorithms.expression_association import ExpressionAssociationAlgorithm
from cfh.genes.registry import load_gene_config
from cfh.model.algorithm_result import AlgorithmResult
from cfh.model.fusion_event import FusionEvent
from cfh.model.fusion_feature import FusionFeature


def test_algorithm_registered():
    assert registry.get("expression_association") is ExpressionAssociationAlgorithm
    assert "expression_association" in registry.list_algorithms()


def test_no_expression_data_is_a_clean_noop_not_a_raise():
    events = [FusionEvent(Event_id="e0", Cohort="c", Sample_id="S1")]

    result = ExpressionAssociationAlgorithm().run(events, [], load_gene_config("RET"), {})

    assert isinstance(result, AlgorithmResult)
    assert result.Algorithm == "expression_association"
    assert result.Summary == {}
    assert result.Warnings == [
        "No mRNA expression data was available for RET in this cohort; "
        "expression-association analysis was skipped."
    ]


def test_empty_expression_by_sample_dict_is_also_a_clean_noop():
    events = [FusionEvent(Event_id="e0", Cohort="c", Sample_id="S1")]
    result = ExpressionAssociationAlgorithm().run(
        events, [], load_gene_config("RET"), {"expression_by_sample": {}}
    )
    assert result.Summary == {}
    assert result.Warnings


def test_fusion_positive_vs_negative_uses_welch_when_both_groups_are_large_and_normal():
    events = [FusionEvent(Event_id=f"pos{i}", Cohort="c", Sample_id=f"POS{i}") for i in range(6)]
    # Fusion-positive samples clearly overexpress relative to negatives, and
    # both groups are large/roughly-normal-shaped so Welch's t-test applies.
    positive_values = [4.0, 4.2, 3.8, 4.1, 3.9, 4.3]
    negative_values = [0.1, -0.2, 0.0, 0.3, -0.1, 0.2, 0.05, -0.05]
    expression_by_sample = {
        **{f"POS{i}": value for i, value in enumerate(positive_values)},
        **{f"NEG{i}": value for i, value in enumerate(negative_values)},
    }
    cohort_sample_ids = list(expression_by_sample.keys())

    result = ExpressionAssociationAlgorithm().run(
        events,
        [],
        None,
        {
            "expression_by_sample": expression_by_sample,
            "cohort_sample_ids": cohort_sample_ids,
        },
    )

    block = result.Summary["fusion_positive_vs_negative"]
    assert block["n_fusion_positive"] == 6
    assert block["n_fusion_negative"] == 8
    assert block["test"] == "welch_t_test"
    assert block["p_value"] < 0.01
    assert block["mean_a"] > block["mean_b"]


def test_fusion_positive_vs_negative_falls_back_to_mann_whitney_below_shapiro_minimum():
    """Below n=3 per group, Shapiro-Wilk cannot assess normality, so the
    algorithm uses Mann-Whitney U unconditionally rather than Welch's."""
    events = [
        FusionEvent(Event_id="e0", Cohort="c", Sample_id="POS0"),
        FusionEvent(Event_id="e1", Cohort="c", Sample_id="POS1"),
    ]
    expression_by_sample = {"POS0": 5.0, "POS1": 5.2, "NEG0": 0.1, "NEG1": 0.2}
    result = ExpressionAssociationAlgorithm().run(
        events,
        [],
        None,
        {
            "expression_by_sample": expression_by_sample,
            "cohort_sample_ids": list(expression_by_sample.keys()),
        },
    )
    block = result.Summary["fusion_positive_vs_negative"]
    assert block["n_fusion_positive"] == 2
    assert block["n_fusion_negative"] == 2
    assert block["test"] == "mann_whitney_u"


def test_single_fusion_positive_sample_skips_the_comparison_not_raises():
    events = [FusionEvent(Event_id="e0", Cohort="c", Sample_id="POS0")]
    expression_by_sample = {"POS0": 5.0, "NEG0": 0.1}
    result = ExpressionAssociationAlgorithm().run(
        events,
        [],
        None,
        {
            "expression_by_sample": expression_by_sample,
            "cohort_sample_ids": list(expression_by_sample.keys()),
        },
    )
    # Only 1 fusion-positive sample -- below the minimum of 2, so this
    # comparison is skipped rather than raising.
    assert "fusion_positive_vs_negative" not in result.Summary
    assert any("fusion-positive-vs-negative" in warning for warning in result.Warnings)


def test_fusion_positive_vs_negative_skipped_without_cohort_sample_ids():
    events = [FusionEvent(Event_id="e0", Cohort="c", Sample_id="POS0")]
    result = ExpressionAssociationAlgorithm().run(
        events,
        [],
        None,
        {"expression_by_sample": {"POS0": 1.0, "NEG0": 0.0}},
    )
    assert "fusion_positive_vs_negative" not in result.Summary
    assert any("no cohort_sample_ids" in warning for warning in result.Warnings)


def test_domain_retention_split_reuses_confidence_stats_default_grouping():
    config = load_gene_config("RET")
    statuses = ("retained", "retained", "retained", "lost", "disrupted", "unknown")
    events = [
        FusionEvent(Event_id=f"e{i}", Cohort="c", Sample_id=f"S{i}") for i in range(len(statuses))
    ]
    features = [
        FusionFeature(
            Event_id=event.Event_id, Gene="RET", Domain_retention_flags={"kinase": status}
        )
        for event, status in zip(events, statuses, strict=True)
    ]
    expression_by_sample = {
        "S0": 3.0,
        "S1": 3.2,
        "S2": 2.9,  # retained group: high expression
        "S3": 0.1,
        "S4": -0.2,  # lost/disrupted -> not_retained group: low expression
        "S5": 9.9,  # unknown -- excluded from both groups
    }

    result = ExpressionAssociationAlgorithm().run(
        events, features, config, {"expression_by_sample": expression_by_sample}
    )

    split = result.Summary["domain_retention_split"]
    assert split["group_a_label"] == "retained"
    assert split["group_b_label"] == "not_retained"
    assert split["n_a"] == 3
    assert split["n_b"] == 2
    assert split["mean_a"] > split["mean_b"]


def test_domain_split_counts_each_sample_once_and_excludes_conflicting_statuses():
    config = load_gene_config("RET")
    sample_statuses = [
        ("S0", "retained"),
        ("S0", "retained"),
        ("S1", "retained"),
        ("S2", "retained"),
        ("S3", "lost"),
        ("S4", "disrupted"),
        ("S5", "retained"),
        ("S5", "lost"),
    ]
    events = [
        FusionEvent(Event_id=f"e{i}", Cohort="c", Sample_id=sample_id)
        for i, (sample_id, _) in enumerate(sample_statuses)
    ]
    features = [
        FusionFeature(
            Event_id=event.Event_id, Gene="RET", Domain_retention_flags={"kinase": status}
        )
        for event, (_, status) in zip(events, sample_statuses, strict=True)
    ]
    result = ExpressionAssociationAlgorithm().run(
        events,
        features,
        config,
        {"expression_by_sample": {f"S{i}": float(i) for i in range(6)}},
    )

    split = result.Summary["domain_retention_split"]
    assert split["analysis_unit"] == "sample"
    assert split["n_a"] == 3
    assert split["n_b"] == 2
    assert split["n_ambiguous_samples_excluded"] == 1
    assert any(
        "1 sample with both retained and not-retained" in warning for warning in result.Warnings
    )


def test_expression_warns_when_known_patient_has_repeated_profiled_biopsies():
    events = [
        FusionEvent(Event_id="e0", Cohort="c", Patient_id="P0", Sample_id="POS0"),
        FusionEvent(Event_id="e1", Cohort="c", Patient_id="P0", Sample_id="POS1"),
    ]
    expression_by_sample = {"POS0": 5.0, "POS1": 5.2, "NEG0": 0.1, "NEG1": 0.2}
    result = ExpressionAssociationAlgorithm().run(
        events,
        [],
        None,
        {
            "expression_by_sample": expression_by_sample,
            "cohort_sample_ids": list(expression_by_sample),
        },
    )

    assert result.Summary["fusion_positive_vs_negative"]["analysis_unit"] == "sample"
    assert result.Summary["fusion_positive_vs_negative"]["n_fusion_positive"] == 2
    assert any(
        "1 known patient has multiple profiled samples" in warning for warning in result.Warnings
    )


def test_domain_retention_split_skipped_when_gene_has_no_key_domains():
    from cfh.genes.registry import GeneConfig

    config = GeneConfig(
        gene_symbol="NOKINASE",
        canonical_transcript_id="NM_000000",
        protein_id="P00000",
    )
    events = [FusionEvent(Event_id="e0", Cohort="c", Sample_id="S0")]
    result = ExpressionAssociationAlgorithm().run(
        events, [], config, {"expression_by_sample": {"S0": 1.0}}
    )
    assert "domain_retention_split" not in result.Summary
    assert any("no key_domains configured" in warning for warning in result.Warnings)


def test_result_schema_matches_canonical_algorithm_result_fields():
    events = [FusionEvent(Event_id="e0", Cohort="c", Sample_id="S0")]
    result = ExpressionAssociationAlgorithm().run(
        events, [], None, {"expression_by_sample": {"S0": 1.0}}
    )
    assert set(type(result).model_fields) == set(AlgorithmResult.model_fields)


def test_expression_by_sample_is_never_echoed_into_parameters():
    """The per-sample expression map can be large (a whole cohort) and must
    never be dumped verbatim into the result's Parameters -- unlike
    algorithms whose params are small/inspectable by design."""
    events = [FusionEvent(Event_id="e0", Cohort="c", Sample_id="S0")]
    big_map = {f"S{i}": float(i) for i in range(500)}
    big_map["S0"] = 1.0
    result = ExpressionAssociationAlgorithm().run(
        events, [], None, {"expression_by_sample": big_map, "cohort_sample_ids": list(big_map)}
    )
    assert "expression_by_sample" not in result.Parameters
    assert "cohort_sample_ids" not in result.Parameters
