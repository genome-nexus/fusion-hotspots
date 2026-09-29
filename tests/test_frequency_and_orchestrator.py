import json
import time

import pytest

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

    assert sum(row["Event_count"] for row in all_events.Tables["Partner_gene_counts"]) == len(
        events
    )
    assert all_events.Tables["Partner_gene_counts"] != deduplicated.Tables["Partner_gene_counts"]
    assert sum(row["Event_count"] for row in deduplicated.Tables["Partner_gene_counts"]) == 3


def test_frequency_patient_dedup_keeps_distinct_partners_and_reports_identity_units():
    events = _events()
    events[1] = events[1].model_copy(update={"Site2_gene": "NTRK1"})
    events.extend(
        [
            FusionEvent(
                Event_id="event-5",
                Cohort="test",
                Sample_id="sample-5",
                Site1_gene="BRAF",
                Site2_gene="AGK",
            ),
            FusionEvent(
                Event_id="event-6",
                Cohort="test",
                Sample_id="sample-5",
                Site1_gene="BRAF",
                Site2_gene="AGK",
            ),
            FusionEvent(
                Event_id="event-7",
                Cohort="test",
                Site1_gene="BRAF",
                Site2_gene="AGK",
            ),
        ]
    )

    result = FrequencyAnalysis().run(events, [], _gene_config(), {"dedup_by_patient": True})
    rows = {row["Partner_gene"]: row for row in result.Tables["Partner_gene_counts"]}

    assert rows["NTRK1"]["Event_count"] == 1
    assert rows["AGK"]["Event_count"] == 3  # patient, sample fallback, event fallback
    assert rows["AGK"]["Raw_event_count"] == 4
    assert rows["AGK"]["Sample_count"] == 2
    assert rows["AGK"]["Patient_count"] == 1
    assert result.Summary["input_event_count"] == 7
    assert result.Summary["input_sample_count"] == 5
    assert result.Summary["input_patient_count"] == 3
    assert result.Summary["analyzed_event_count"] == 6
    assert result.Summary["count_unit"] == "patient_partner"
    assert "Sample_id for 2 and Event_id for 1" in result.Warnings[0]


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


def test_orchestrator_defers_a_dependent_algorithm_and_injects_dependency_results():
    """A DEPENDS_ON-declaring algorithm must run in a later wave than its
    also-requested dependency, with that dependency's real AlgorithmResult
    automatically merged into its own params -- the real orchestrator wiring
    composite_score relies on, exercised here against a minimal synthetic
    plugin instead of composite_score's own more complex logic.
    """

    @register("consumer")
    class Consumer(Algorithm):
        DEPENDS_ON = ("frequency",)

        def run(self, events, features, gene_config, params) -> AlgorithmResult:
            supplied = params.get("algorithm_results") or []
            frequency_result = next((r for r in supplied if r.Algorithm == "frequency"), None)
            total = None
            if frequency_result is not None:
                total = sum(
                    row["Event_count"] for row in frequency_result.Tables["Partner_gene_counts"]
                )
            return AlgorithmResult(Algorithm="consumer", Summary={"total_from_frequency": total})

    try:
        results = run_algorithms(["consumer", "frequency"], _events(), [], _gene_config())
    finally:
        unregister("consumer")

    by_name = {result.Algorithm: result for result in results}
    assert by_name["consumer"].Summary["total_from_frequency"] == len(_events())
    # Output stays in request order even though "consumer" was scheduled
    # into a later dependency-respecting wave than "frequency".
    assert [result.Algorithm for result in results] == ["consumer", "frequency"]


def test_orchestrator_leaves_a_dependency_unavailable_when_not_requested():
    """A declared dependency that was not itself requested must not block
    scheduling or crash the dependent algorithm -- it is simply unavailable,
    exactly like composite_score's own graceful-exclusion behavior.
    """

    @register("consumer")
    class Consumer(Algorithm):
        DEPENDS_ON = ("frequency",)

        def run(self, events, features, gene_config, params) -> AlgorithmResult:
            supplied = params.get("algorithm_results") or []
            return AlgorithmResult(Algorithm="consumer", Summary={"n_supplied": len(supplied)})

    try:
        results = run_algorithms(["consumer"], _events(), [], _gene_config())
    finally:
        unregister("consumer")

    assert results[0].Summary["n_supplied"] == 0
    assert results[0].Warnings == []


def test_orchestrator_extra_results_seed_dependency_injection_without_rerunning():
    """``extra_results`` lets a caller that pre-computed one algorithm
    separately (as real_benchmark.py does for domain_retention) make that
    result available for dependency injection without it being re-run or
    re-returned by this call.
    """

    @register("consumer")
    class Consumer(Algorithm):
        DEPENDS_ON = ("frequency",)

        def run(self, events, features, gene_config, params) -> AlgorithmResult:
            supplied = params.get("algorithm_results") or []
            frequency_result = next((r for r in supplied if r.Algorithm == "frequency"), None)
            return AlgorithmResult(
                Algorithm="consumer",
                Summary={"saw_frequency": frequency_result is not None},
            )

    precomputed_frequency = AlgorithmResult(
        Algorithm="frequency",
        Tables={"Partner_gene_counts": [{"Partner_gene": "X", "Event_count": 7}]},
    )
    try:
        results = run_algorithms(
            ["consumer"],
            _events(),
            [],
            _gene_config(),
            extra_results=[precomputed_frequency],
        )
    finally:
        unregister("consumer")

    assert results == [
        result for result in results if result.Algorithm != "frequency"
    ]  # frequency was not re-run
    assert results[0].Summary["saw_frequency"] is True


def test_orchestrator_respects_caller_supplied_algorithm_results_over_auto_injection():
    """If the caller already set ``params[name]["algorithm_results"]``
    explicitly, the orchestrator must not overwrite it with auto-injected
    dependency results.
    """

    @register("consumer")
    class Consumer(Algorithm):
        DEPENDS_ON = ("frequency",)

        def run(self, events, features, gene_config, params) -> AlgorithmResult:
            supplied = params.get("algorithm_results") or []
            return AlgorithmResult(
                Algorithm="consumer",
                Summary={"algorithms_seen": sorted(r.Algorithm for r in supplied)},
            )

    explicit_override = AlgorithmResult(Algorithm="custom_source", Summary={})
    try:
        results = run_algorithms(
            ["consumer", "frequency"],
            _events(),
            [],
            _gene_config(),
            params={"consumer": {"algorithm_results": [explicit_override]}},
        )
    finally:
        unregister("consumer")

    by_name = {result.Algorithm: result for result in results}
    assert by_name["consumer"].Summary["algorithms_seen"] == ["custom_source"]


def test_orchestrator_raises_on_circular_depends_on_instead_of_silently_degrading():
    """A -> B -> A can never make scheduling progress. Silently running
    both in one wave (as an earlier implementation did) would inject EMPTY
    dependency results into each -- a silent DEPENDS_ON violation that
    looks like a normal run. This must instead raise a clear, explicit
    configuration error naming the algorithms involved.
    """

    @register("cycle-a")
    class CycleA(Algorithm):
        DEPENDS_ON = ("cycle-b",)

        def run(self, events, features, gene_config, params) -> AlgorithmResult:
            return AlgorithmResult(Algorithm="cycle-a")

    @register("cycle-b")
    class CycleB(Algorithm):
        DEPENDS_ON = ("cycle-a",)

        def run(self, events, features, gene_config, params) -> AlgorithmResult:
            return AlgorithmResult(Algorithm="cycle-b")

    try:
        with pytest.raises(ValueError, match="circular"):
            run_algorithms(["cycle-a", "cycle-b"], _events(), [], _gene_config())
    finally:
        unregister("cycle-a")
        unregister("cycle-b")
