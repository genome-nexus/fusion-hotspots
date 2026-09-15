from cfh.algorithms import registry
from cfh.algorithms.mutation_cooccurrence import (
    MutationCooccurrenceAlgorithm,
    comparator_target_key,
)
from cfh.genes.registry import GeneConfig, MutualExclusivityTarget, load_gene_config
from cfh.model.algorithm_result import AlgorithmResult
from cfh.model.fusion_event import FusionEvent

_TARGET_GENE_WITH_COMPARATOR = GeneConfig(
    gene_symbol="FAKE1",
    canonical_transcript_id="NM_000001",
    protein_id="P00001",
    mutual_exclusivity_targets=[
        MutualExclusivityTarget(
            gene="FAKE1",
            alteration_type="point_mutation",
            entrez_gene_id=99999,
            protein_change="X100Y",
        )
    ],
)

_TARGET_GENE_NO_COMPARATOR = GeneConfig(
    gene_symbol="FAKE2",
    canonical_transcript_id="NM_000002",
    protein_id="P00002",
)


def _event(event_id: str, sample_id: str) -> FusionEvent:
    return FusionEvent(Event_id=event_id, Cohort="synthetic", Sample_id=sample_id)


def _comparator_row(sample_id: str, protein_change: str = "X100Y") -> dict:
    return {
        "Sample_id": sample_id,
        "Gene": "FAKE1",
        "Alteration_type": "point_mutation",
        "Protein_change": protein_change,
    }


def test_algorithm_registered():
    assert registry.get("mutation_cooccurrence") is MutationCooccurrenceAlgorithm
    assert "mutation_cooccurrence" in registry.list_algorithms()


def test_no_op_when_gene_leaves_mutual_exclusivity_targets_unset():
    """Same graceful-skip pattern already proven for domain_disruption/
    exon_retention/joint_partner: a gene that never opts in produces a
    documented no-op, not an error."""
    result = MutationCooccurrenceAlgorithm().run([], [], _TARGET_GENE_NO_COMPARATOR, {})

    assert result.Algorithm == "mutation_cooccurrence"
    assert result.Summary == {"targets": []}
    assert result.Tables == {}
    assert result.Warnings == [
        "FAKE2 has no mutual_exclusivity_targets configured; co-occurrence/"
        "mutual-exclusivity analysis was skipped."
    ]
    assert set(type(result).model_fields) == set(AlgorithmResult.model_fields)


def test_no_op_when_configured_but_no_cohort_data_supplied():
    """A gene CAN opt in without a caller supplying the live-fetched
    cohort/comparator data (e.g. the default orchestrator sweep) -- this
    must still no-op gracefully with a distinct warning, never crash."""
    result = MutationCooccurrenceAlgorithm().run([], [], _TARGET_GENE_WITH_COMPARATOR, {})

    assert result.Summary == {"targets": []}
    assert result.Warnings
    assert "no cohort_sample_ids were supplied" in result.Warnings[0]


def test_detects_clean_mutual_exclusivity_on_synthetic_data():
    cohort_sample_ids = [f"S{i}" for i in range(1, 21)]
    fusion_events = [_event(f"evt-{i}", f"S{i}") for i in range(1, 9)]  # S1-S8 fusion+
    comparator_alterations = [_comparator_row(f"S{i}") for i in range(9, 17)]  # S9-S16 mutated

    result = MutationCooccurrenceAlgorithm().run(
        fusion_events,
        [],
        _TARGET_GENE_WITH_COMPARATOR,
        {
            "cohort_sample_ids": cohort_sample_ids,
            "comparator_alterations": comparator_alterations,
        },
    )

    assert result.Algorithm == "mutation_cooccurrence"
    row = result.Summary["targets"][0]
    assert row["contingency_table"] == [[0, 8], [8, 4]]
    assert row["direction"] == "mutually_exclusive"
    assert row["odds_ratio"] == 0.0
    assert row["p_value"] < 0.05
    assert row["comparator_gene"] == "FAKE1"
    assert row["protein_change"] == "X100Y"
    assert result.Tables["cooccurrence_results"] == result.Summary["targets"]
    assert set(type(result).model_fields) == set(AlgorithmResult.model_fields)


def test_comparator_alterations_filtered_by_gene_alteration_type_and_protein_change():
    cohort_sample_ids = [f"S{i}" for i in range(1, 11)]
    fusion_events = [_event("evt-1", "S1")]
    comparator_alterations = [
        _comparator_row("S2"),  # matches
        {**_comparator_row("S3"), "Protein_change": "OTHER"},  # wrong protein change
        {**_comparator_row("S4"), "Gene": "OTHERGENE"},  # wrong gene
        {**_comparator_row("S5"), "Alteration_type": "cna_amp"},  # wrong alteration type
    ]

    result = MutationCooccurrenceAlgorithm().run(
        fusion_events,
        [],
        _TARGET_GENE_WITH_COMPARATOR,
        {
            "cohort_sample_ids": cohort_sample_ids,
            "comparator_alterations": comparator_alterations,
        },
    )

    row = result.Summary["targets"][0]
    assert row["comparator_altered_sample_count"] == 1


def test_no_comparator_samples_still_computes_and_warns():
    cohort_sample_ids = [f"S{i}" for i in range(1, 6)]
    fusion_events = [_event("evt-1", "S1")]

    result = MutationCooccurrenceAlgorithm().run(
        fusion_events,
        [],
        _TARGET_GENE_WITH_COMPARATOR,
        {"cohort_sample_ids": cohort_sample_ids, "comparator_alterations": []},
    )

    row = result.Summary["targets"][0]
    assert row["comparator_altered_sample_count"] == 0
    assert any("No cohort samples carried" in warning for warning in result.Warnings)


def test_successful_zero_data_target_with_availability_sidecar_still_computes():
    target = _TARGET_GENE_WITH_COMPARATOR.mutual_exclusivity_targets[0]
    result = MutationCooccurrenceAlgorithm().run(
        [_event("evt-1", "S1")],
        [],
        _TARGET_GENE_WITH_COMPARATOR,
        {
            "cohort_sample_ids": ["S1", "S2"],
            "comparator_alterations": [],
            "comparator_availability": {comparator_target_key(target): True},
        },
    )

    assert result.Summary["targets"][0]["comparator_altered_sample_count"] == 0
    assert any("No cohort samples carried" in warning for warning in result.Warnings)


def test_none_comparator_data_skips_instead_of_reporting_confirmed_zero_calls():
    result = MutationCooccurrenceAlgorithm().run(
        [_event("evt-1", "S1")],
        [],
        _TARGET_GENE_WITH_COMPARATOR,
        {"cohort_sample_ids": ["S1", "S2"], "comparator_alterations": None},
    )

    assert result.Summary == {"targets": []}
    assert result.Tables == {}
    assert "comparator alteration data were unavailable" in result.Warnings[0]


def test_availability_sidecar_skips_only_failed_comparator_targets():
    second_target = MutualExclusivityTarget(
        gene="FAKE2",
        alteration_type="point_mutation",
        entrez_gene_id=88888,
    )
    config = _TARGET_GENE_WITH_COMPARATOR.model_copy(
        update={
            "mutual_exclusivity_targets": [
                _TARGET_GENE_WITH_COMPARATOR.mutual_exclusivity_targets[0],
                second_target,
            ]
        }
    )

    result = MutationCooccurrenceAlgorithm().run(
        [_event("evt-1", "S1")],
        [],
        config,
        {
            "cohort_sample_ids": ["S1", "S2", "S3"],
            "comparator_alterations": [_comparator_row("S2")],
            "comparator_availability": {
                comparator_target_key(config.mutual_exclusivity_targets[0]): True,
                comparator_target_key(second_target): False,
            },
        },
    )

    assert [row["comparator_gene"] for row in result.Summary["targets"]] == ["FAKE1"]
    assert any("FAKE2 point_mutation were unavailable" in warning for warning in result.Warnings)


def test_sidecar_missing_a_target_fails_closed():
    result = MutationCooccurrenceAlgorithm().run(
        [_event("evt-1", "S1")],
        [],
        _TARGET_GENE_WITH_COMPARATOR,
        {
            "cohort_sample_ids": ["S1", "S2"],
            "comparator_alterations": [],
            "comparator_availability": {},
        },
    )

    assert result.Summary == {"targets": []}
    assert any("were unavailable" in warning for warning in result.Warnings)


def test_real_braf_config_opts_in_and_ret_does_not():
    """The real curated configs: BRAF's own V600E hotspot is configured as
    a comparator (the biology this evidence layer was built to test);
    RET's config never opted in, so it stays a clean no-op -- proving this
    feature is genuinely gene-agnostic/opt-in rather than hardcoded."""
    braf_config = load_gene_config("braf")
    ret_config = load_gene_config("ret")

    assert braf_config.mutual_exclusivity_targets == [
        MutualExclusivityTarget(
            gene="BRAF",
            alteration_type="point_mutation",
            entrez_gene_id=673,
            protein_change="V600E",
        )
    ]
    assert ret_config.mutual_exclusivity_targets == []

    ret_result = MutationCooccurrenceAlgorithm().run([], [], ret_config, {})
    assert ret_result.Summary == {"targets": []}
