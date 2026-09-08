from cfh.algorithms import registry
from cfh.algorithms.window_detection import WindowDetectionAlgorithm, _mapping_sensitivity_summary
from cfh.genes.registry import GeneConfig, KeyDomain
from cfh.model.algorithm_result import AlgorithmResult
from cfh.model.fusion_event import FusionEvent
from cfh.model.fusion_feature import FusionFeature
from cfh.orchestrator.run import run_algorithms

# A fabricated single-gene config, deliberately unrelated to any real gene
# in genes/configs/, to demonstrate the algorithm has no gene-specific
# hardcoding of its own.
_FAKE_GENE = GeneConfig(
    gene_symbol="FAKE1",
    canonical_transcript_id="NM_000001",
    protein_id="P00001",
    key_domains=[KeyDomain(name="Made-up domain", source="test", key="made_up")],
)


def _events_and_features(records: list[tuple[str, int, str]]) -> tuple[list, list]:
    events = []
    features = []
    for event_id, position, status in records:
        events.append(FusionEvent(Event_id=event_id, Cohort="synthetic"))
        features.append(
            FusionFeature(
                Event_id=event_id,
                Gene="FAKE1",
                Junction_position_aa=position,
                Domain_retention_flags={"made_up": status},
            )
        )
    return events, features


def _internal_region_records() -> list[tuple[str, int, str]]:
    """Retained breakpoints cluster inside [500, 550]; lost breakpoints sit
    well outside that band on both sides -- an internal region bounded on
    both sides, unlike a single terminal cutpoint (which cutpoint_detection
    already covers)."""
    retained_positions = [505, 510, 520, 530, 545]
    lost_positions = [100, 150, 700, 750, 800]
    retained = [(f"r{i}", position, "retained") for i, position in enumerate(retained_positions)]
    lost = [(f"l{i}", position, "lost") for i, position in enumerate(lost_positions)]
    return retained + lost


def _clamped_pile_records() -> list[tuple[str, int, str]]:
    """Most events clamp onto a single position (the historical ALK-style
    intronic-clamping artifact) -- many candidate windows around that pile
    contain the exact same set of events."""
    clamped = [(f"c{i}", 400, "retained") for i in range(8)]
    other = [
        ("o0", 50, "lost"),
        ("o1", 60, "lost"),
        ("o2", 900, "lost"),
        ("o3", 950, "lost"),
    ]
    return clamped + other


def test_algorithm_registered():
    assert registry.get("window_detection") is WindowDetectionAlgorithm
    assert "window_detection" in registry.list_algorithms()


def test_recovers_internal_window_from_synthetic_clean_separation():
    events, features = _events_and_features(_internal_region_records())

    result = WindowDetectionAlgorithm().run(
        events, features, _FAKE_GENE, {"seed": 42, "n_permutations": 300}
    )

    assert isinstance(result, AlgorithmResult)
    assert result.Algorithm == "window_detection"
    assert result.Summary["determinable"] is True
    best = result.Summary["best_window"]
    assert best["start_aa"] <= 505
    assert best["end_aa"] >= 545
    assert result.Summary["n_events_analyzed"] == 10
    assert 0.0 <= result.Summary["corrected_p_value"] <= 1.0
    assert result.Warnings == [
        "mapping-sensitivity information was not supplied; the fraction of "
        "events using a clamped/approximate breakpoint position vs. an "
        "exact-genomic-mapped one could not be computed"
    ]
    assert result.Tables["window_scan"]
    assert result.Tables["top_windows"]
    assert result.Summary["known_domain_boundary_comparison"] is None


def test_domain_boundary_comparison_reports_both_edges():
    events, features = _events_and_features(_internal_region_records())

    result = WindowDetectionAlgorithm().run(
        events,
        features,
        _FAKE_GENE,
        {"seed": 42, "n_permutations": 300, "domain_boundaries": [503, 547, 900]},
    )

    comparison = result.Summary["known_domain_boundary_comparison"]
    assert comparison["start_aa"]["nearest_known_domain_boundary_aa"] == 503
    assert comparison["end_aa"]["nearest_known_domain_boundary_aa"] == 547


def test_mapping_sensitivity_unavailable_by_default():
    events, features = _events_and_features(_internal_region_records())

    result = WindowDetectionAlgorithm().run(
        events, features, _FAKE_GENE, {"seed": 42, "n_permutations": 100}
    )

    sensitivity = result.Summary["mapping_sensitivity"]
    assert sensitivity["available"] is False
    assert sensitivity["fraction_clamped_or_approximate"] is None
    assert result.Summary["best_window_mapping_sensitivity"]["available"] is False


def test_mapping_sensitivity_reports_clamped_fraction_when_supplied():
    events, features = _events_and_features(_internal_region_records())
    # Every "retained" (inside-window) event was clamped/approximate;
    # every "lost" (outside) event was mapped exactly.
    mapping_sensitivity = {f"r{i}": True for i in range(5)} | {f"l{i}": False for i in range(5)}

    result = WindowDetectionAlgorithm().run(
        events,
        features,
        _FAKE_GENE,
        {"seed": 42, "n_permutations": 300, "mapping_sensitivity": mapping_sensitivity},
    )

    overall = result.Summary["mapping_sensitivity"]
    assert overall["available"] is True
    assert overall["n_exact"] == 5
    assert overall["n_clamped_or_approximate"] == 5
    assert overall["fraction_clamped_or_approximate"] == 0.5

    best_window_sensitivity = result.Summary["best_window_mapping_sensitivity"]
    assert best_window_sensitivity["available"] is True
    assert best_window_sensitivity["fraction_clamped_or_approximate"] == 1.0
    assert result.Warnings == []


def test_mapping_sensitivity_summary_treats_missing_events_as_unknown():
    summary = _mapping_sensitivity_summary(["a", "b", "c"], {"a": True})
    assert summary["available"] is True
    assert summary["n_clamped_or_approximate"] == 1
    assert summary["n_unknown"] == 2
    assert summary["fraction_clamped_or_approximate"] == 1.0


def test_top_windows_table_dedups_a_clamped_pile():
    events, features = _events_and_features(_clamped_pile_records())

    result = WindowDetectionAlgorithm().run(
        events, features, _FAKE_GENE, {"seed": 1, "n_permutations": 200}
    )

    assert result.Summary["determinable"] is True
    top_windows = result.Tables["top_windows"]
    masks = [frozenset(row["event_ids_inside"]) for row in top_windows]
    assert len(masks) == len(set(masks))
    assert (
        result.Summary["n_distinct_candidate_windows_by_event_mask"]
        < result.Summary["n_candidate_windows"]
    )


def test_degenerate_input_is_reported_as_not_determinable_without_raising():
    events, features = _events_and_features(
        [("a", 100, "retained"), ("b", 110, "retained"), ("c", 120, "retained")]
    )

    result = WindowDetectionAlgorithm().run(
        events, features, _FAKE_GENE, {"seed": 42, "n_permutations": 50}
    )

    assert result.Summary["determinable"] is False
    assert result.Summary["best_window"] is None
    assert result.Warnings
    assert result.Summary["known_domain_boundary_comparison"] is None
    assert result.Summary["best_window_mapping_sensitivity"] is None


def test_gene_with_no_key_domain_is_a_graceful_no_op_through_the_orchestrator():
    """A gene with no configured key domain has nothing for this algorithm to
    scan (there is no domain-retention status to separate on). Mirrors
    ``cutpoint_detection``'s existing behavior for the same input shape:
    the per-gene lookup this relies on raises, and it is the orchestrator's
    job (:func:`cfh.orchestrator.run.run_algorithms`) to turn that into a
    non-crashing, warning-carrying result rather than aborting the run --
    exactly the same contract every other algorithm here already has."""
    no_domain_gene = GeneConfig(
        gene_symbol="NODOM1", canonical_transcript_id="NM_999999", protein_id="P99999"
    )
    events = [FusionEvent(Event_id="a", Cohort="synthetic")]
    features = [
        FusionFeature(
            Event_id="a",
            Gene="NODOM1",
            Junction_position_aa=100,
            Domain_retention_flags={"made_up": "retained"},
        )
    ]

    (result,) = run_algorithms(["window_detection"], events, features, no_domain_gene, {})

    assert isinstance(result, AlgorithmResult)
    assert result.Algorithm == "window_detection"
    assert result.Warnings
    assert "no configured domains" in result.Warnings[0].lower()


def test_events_missing_mapped_positions_are_excluded_not_crashed_on():
    events = [FusionEvent(Event_id="unmapped", Cohort="synthetic")]
    features = [
        FusionFeature(
            Event_id="unmapped",
            Gene="FAKE1",
            Junction_position_aa=None,
            Domain_retention_flags={"made_up": "retained"},
        )
    ]

    result = WindowDetectionAlgorithm().run(
        events, features, _FAKE_GENE, {"seed": 42, "n_permutations": 50}
    )

    assert result.Summary["n_events_analyzed"] == 0
    assert result.Summary["determinable"] is False


def test_result_schema_matches_canonical_algorithm_result_fields():
    events, features = _events_and_features(_internal_region_records())
    result = WindowDetectionAlgorithm().run(
        events, features, _FAKE_GENE, {"seed": 42, "n_permutations": 50}
    )
    assert set(type(result).model_fields) == set(AlgorithmResult.model_fields)


def test_custom_widths_and_min_events_are_honored():
    events, features = _events_and_features(_internal_region_records())

    result = WindowDetectionAlgorithm().run(
        events,
        features,
        _FAKE_GENE,
        {"seed": 42, "n_permutations": 100, "widths": [10, 20], "min_events_per_window": 5},
    )

    assert result.Summary["widths_tested_aa"] == [10, 20]
    assert result.Summary["min_events_per_window"] == 5
    for row in result.Tables["window_scan"]:
        assert row["width_aa"] in (10, 20)
        assert row["n_events_inside"] >= 5
        assert row["n_events_outside"] >= 5
