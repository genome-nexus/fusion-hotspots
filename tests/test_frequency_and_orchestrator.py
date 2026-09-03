import json
import time

from cfh.algorithms import ExonRetentionAnalysis, FrequencyAnalysis
from cfh.algorithms.base import Algorithm
from cfh.algorithms.registry import register, unregister
from cfh.genes.registry import GeneConfig
from cfh.model.algorithm_result import AlgorithmResult
from cfh.model.fusion_event import FusionEvent
from cfh.model.fusion_feature import FusionFeature
from cfh.orchestrator.run import _params_for, results_to_json, run_algorithms


def _gene_config() -> GeneConfig:
    return GeneConfig(
        gene_symbol="BRAF",
        canonical_transcript_id="NM_004333",
        protein_id="P15056",
    )


def _events() -> list[FusionEvent]:
    return [
        FusionEvent(
            Event_id="event-1",
            Cohort="test",
            Patient_id="patient-1",
            Sample_id="sample-1",
            Site1_gene="BRAF",
            Site2_gene="AGK",
        ),
        FusionEvent(
            Event_id="event-2",
            Cohort="test",
            Patient_id="patient-1",
            Sample_id="sample-2",
            Site1_gene="BRAF",
            Site2_gene="AGK",
        ),
        FusionEvent(
            Event_id="event-3",
            Cohort="test",
            Patient_id="patient-2",
            Sample_id="sample-3",
            Site1_gene="BRAF",
            Site2_gene="KIAA1549",
        ),
        FusionEvent(
            Event_id="event-4",
            Cohort="test",
            Patient_id="patient-3",
            Sample_id="sample-4",
            Site1_gene="BRAF",
            Site2_gene="SND1",
        ),
    ]


def test_frequency_counts_all_events_and_patient_deduplication_changes_counts():
    events = _events()
    config = _gene_config()

    all_events = FrequencyAnalysis().run(events, [], config, {"dedup_by_patient": False})
    deduplicated = FrequencyAnalysis().run(events, [], config, {"dedup_by_patient": True})

    assert (
        sum(row["Event_count"] for row in all_events.Tables["Partner_gene_counts"])
        == len(events)
    )
    assert all_events.Tables["Partner_gene_counts"] != deduplicated.Tables["Partner_gene_counts"]
    assert sum(row["Event_count"] for row in deduplicated.Tables["Partner_gene_counts"]) == 3


def test_exon_retention_fraction_is_exact_for_hand_built_features():
    events = _events()
    features = [
        FusionFeature(Event_id=event.Event_id, Gene="BRAF", Retained_exons=[8])
        for event in events[:3]
    ] + [FusionFeature(Event_id=events[3].Event_id, Gene="BRAF", Retained_exons=[7])]

    result = ExonRetentionAnalysis().run(events, features, _gene_config(), {"target_exon": 8})

    assert result.Summary["retained_fraction"] == 0.75


def test_orchestrator_runs_plugins_concurrently_and_isolates_failures():
    @register("slow-success")
    class SlowSuccess(Algorithm):
        def run(self, events, features, gene_config, params) -> AlgorithmResult:
            time.sleep(0.2)
            return AlgorithmResult(Algorithm="slow-success", Algorithm_version="test")

    @register("slow-failure")
    class SlowFailure(Algorithm):
        def run(self, events, features, gene_config, params) -> AlgorithmResult:
            time.sleep(0.2)
            raise RuntimeError("intentional failure")

    try:
        started_at = time.perf_counter()
        results = run_algorithms(["slow-success", "slow-failure"], _events(), [], _gene_config())
        elapsed = time.perf_counter() - started_at
    finally:
        unregister("slow-success")
        unregister("slow-failure")

    by_name = {result.Algorithm: result for result in results}
    assert elapsed < 0.35
    assert by_name["slow-success"].Warnings == []
    assert "intentional failure" in by_name["slow-failure"].Warnings[0]
    assert len({result.Input_fingerprint for result in results}) == 1


def test_orchestrator_results_are_json_serializable():
    results = run_algorithms(
        ["frequency", "exon_retention"],
        _events(),
        [FusionFeature(Event_id="event-1", Gene="BRAF", Retained_exons=[8])],
        _gene_config(),
        params={"exon_retention": {"target_exon": 8}},
    )

    emitted = results_to_json(results)
    assert [result["Algorithm"] for result in json.loads(emitted)] == [
        "frequency",
        "exon_retention",
    ]


def test_params_are_isolated_without_breaking_flat_algorithm_params():
    assert _params_for("frequency", {"exon_retention": {"target_exon": 8}}) == {}
    assert _params_for("frequency", {"dedup_by_patient": True}) == {"dedup_by_patient": True}
