from cfh.algorithms import registry
from cfh.algorithms.genomic_position_recurrence import GenomicPositionRecurrenceAlgorithm
from cfh.genes.registry import GeneConfig
from cfh.model.algorithm_result import AlgorithmResult
from cfh.model.fusion_event import FusionEvent
from cfh.model.fusion_feature import FusionFeature

# A fabricated single-gene config, deliberately unrelated to any real gene in
# genes/configs/, to demonstrate the algorithm has no gene-specific
# hardcoding of its own -- same convention used by the other gene-agnostic
# algorithm test modules (see test_window_detection_algorithm.py).
_FAKE_GENE = GeneConfig(
    gene_symbol="FAKE1",
    canonical_transcript_id="NM_000001",
    protein_id="P00001",
)


def _event_and_feature(event_id: str, aa: int) -> tuple[FusionEvent, FusionFeature]:
    event = FusionEvent(Event_id=event_id, Cohort="synthetic")
    feature = FusionFeature(Event_id=event_id, Gene="FAKE1", Junction_position_aa=aa)
    return event, feature


def _breakpoint(chromosome: str, position: int, build: str = "GRCh37") -> dict:
    return {"chromosome": chromosome, "position": position, "build": build}


def test_algorithm_registered():
    assert registry.get("genomic_position_recurrence") is GenomicPositionRecurrenceAlgorithm
    assert "genomic_position_recurrence" in registry.list_algorithms()


def test_depends_on_declares_cutpoint_detection():
    assert GenomicPositionRecurrenceAlgorithm.DEPENDS_ON == ("cutpoint_detection",)


def test_no_genomic_breakpoints_supplied_is_a_graceful_no_op():
    event, feature = _event_and_feature("a", 439)

    result = GenomicPositionRecurrenceAlgorithm().run([event], [feature], _FAKE_GENE, {})

    assert isinstance(result, AlgorithmResult)
    assert result.Algorithm == "genomic_position_recurrence"
    assert result.Summary["determinable"] is False
    assert result.Summary["n_events_analyzed"] == 0
    assert result.Warnings
    assert result.Tables["genomic_bin_recurrence"] == []


def test_no_target_gene_event_has_a_known_position_is_a_graceful_no_op():
    event, feature = _event_and_feature("a", 439)
    # Genomic breakpoints supplied, but for an event id that never shows up
    # among this gene's mapped features.
    params = {"genomic_breakpoints": {"unrelated-event": _breakpoint("7", 100)}}

    result = GenomicPositionRecurrenceAlgorithm().run([event], [feature], _FAKE_GENE, params)

    assert result.Summary["determinable"] is False
    assert result.Summary["n_events_analyzed"] == 0


def test_events_missing_chromosome_or_position_are_excluded_not_crashed_on():
    events, features = zip(
        *[_event_and_feature("a", 439), _event_and_feature("b", 439), _event_and_feature("c", 439)]
    )
    params = {
        "genomic_breakpoints": {
            "a": _breakpoint("7", 100_000),
            "b": _breakpoint("7", 100_000),
            "c": {"chromosome": None, "position": None, "build": "GRCh37"},
        }
    }

    result = GenomicPositionRecurrenceAlgorithm().run(
        list(events), list(features), _FAKE_GENE, params
    )

    assert result.Summary["determinable"] is True
    assert result.Summary["n_events_analyzed"] == 2
    assert result.Summary["n_events_excluded_missing_position"] == 1


def test_shared_protein_position_with_identical_genomic_breakpoint_is_a_real_hotspot():
    """Three events all mapped to the same protein position AND the same
    exact genomic coordinate -- the "real DNA-level hotspot" case."""
    records = [
        _event_and_feature("h0", 439),
        _event_and_feature("h1", 439),
        _event_and_feature("h2", 439),
        _event_and_feature("s0", 200),
    ]
    events, features = zip(*records)
    params = {
        "genomic_breakpoints": {
            "h0": _breakpoint("7", 100_000),
            "h1": _breakpoint("7", 100_000),
            "h2": _breakpoint("7", 100_000),
            "s0": _breakpoint("7", 500),
        }
    }

    result = GenomicPositionRecurrenceAlgorithm().run(
        list(events), list(features), _FAKE_GENE, params
    )

    assert result.Summary["determinable"] is True
    assert result.Summary["n_events_analyzed"] == 4
    assert result.Summary["reference_build"] == "GRCh37"
    assert result.Summary["n_distinct_genomic_positions"] == 2

    exact_table = result.Tables["exact_position_recurrence"]
    assert exact_table == [
        {
            "chromosome": "7",
            "position": 100_000,
            "n_events": 3,
            "event_ids": ["h0", "h1", "h2"],
        }
    ]

    spread_table = result.Tables["protein_position_genomic_spread"]
    assert spread_table == [
        {
            "protein_position_aa": 439,
            "n_events": 3,
            "n_distinct_genomic_positions": 1,
            "genomic_span_bp": 0,
            "spans_multiple_chromosomes": False,
        }
    ]
    note = result.Summary["genomic_vs_protein_clustering_note"]
    assert "share the exact same genomic breakpoint" in note
    assert "real DNA-level recurrent breakpoint" in note

    top_bin = result.Summary["top_genomic_bin"]
    assert top_bin["n_events"] == 3
    assert top_bin["chromosome"] == "7"


def test_shared_protein_position_with_scattered_genomic_breakpoints_is_not_a_hotspot():
    """Three events share one protein position but land at three distinct,
    widely-spread genomic coordinates -- the mapping/clamping-artifact case
    the aa-439-vs-458 investigation was concerned about."""
    records = [
        _event_and_feature("s0", 439),
        _event_and_feature("s1", 439),
        _event_and_feature("s2", 439),
    ]
    events, features = zip(*records)
    params = {
        "genomic_breakpoints": {
            "s0": _breakpoint("7", 100_000),
            "s1": _breakpoint("7", 105_000),
            "s2": _breakpoint("7", 120_000),
        }
    }

    result = GenomicPositionRecurrenceAlgorithm().run(
        list(events), list(features), _FAKE_GENE, params
    )

    spread_table = result.Tables["protein_position_genomic_spread"]
    assert spread_table == [
        {
            "protein_position_aa": 439,
            "n_events": 3,
            "n_distinct_genomic_positions": 3,
            "genomic_span_bp": 20_000,
            "spans_multiple_chromosomes": False,
        }
    ]
    note = result.Summary["genomic_vs_protein_clustering_note"]
    assert "3 distinct genomic positions spanning 20000 bp" in note
    assert "not one shared DNA lesion" in note or "rather than one shared DNA lesion" in note


def test_clustering_note_anchors_on_cutpoint_detection_result_when_available():
    """Without a cutpoint_detection cross-reference, the note describes the
    *most-recurrent* protein position (aa 300, a real 4-event hotspot).
    Supplying a cutpoint_detection result inferring aa 439 instead switches
    the note to describe aa 439's cluster (scattered, 2 events) -- proving
    the DEPENDS_ON cross-reference actually changes the reported note."""
    records = [
        _event_and_feature("hot0", 300),
        _event_and_feature("hot1", 300),
        _event_and_feature("hot2", 300),
        _event_and_feature("hot3", 300),
        _event_and_feature("sc0", 439),
        _event_and_feature("sc1", 439),
    ]
    events, features = zip(*records)
    genomic_breakpoints = {
        "hot0": _breakpoint("7", 50_000),
        "hot1": _breakpoint("7", 50_000),
        "hot2": _breakpoint("7", 50_000),
        "hot3": _breakpoint("7", 50_000),
        "sc0": _breakpoint("7", 100_000),
        "sc1": _breakpoint("7", 100_100),
    }

    without_cutpoint = GenomicPositionRecurrenceAlgorithm().run(
        list(events), list(features), _FAKE_GENE, {"genomic_breakpoints": genomic_breakpoints}
    )
    default_note = without_cutpoint.Summary["genomic_vs_protein_clustering_note"]
    assert "share the exact same genomic breakpoint" in default_note

    cutpoint_result = AlgorithmResult(
        Algorithm="cutpoint_detection",
        Summary={"determinable": True, "inferred_cutpoint_aa": 439},
    )
    with_cutpoint = GenomicPositionRecurrenceAlgorithm().run(
        list(events),
        list(features),
        _FAKE_GENE,
        {
            "genomic_breakpoints": genomic_breakpoints,
            "algorithm_results": [cutpoint_result],
        },
    )
    anchored_note = with_cutpoint.Summary["genomic_vs_protein_clustering_note"]
    assert "2 distinct genomic positions spanning 100 bp" in anchored_note


def test_reference_build_mismatch_excludes_minority_build_and_warns():
    records = [
        _event_and_feature("a", 439),
        _event_and_feature("b", 439),
        _event_and_feature("c", 439),
    ]
    events, features = zip(*records)
    params = {
        "genomic_breakpoints": {
            "a": _breakpoint("7", 100_000, build="GRCh37"),
            "b": _breakpoint("7", 100_000, build="GRCh37"),
            "c": _breakpoint("7", 999_000, build="GRCh38"),
        }
    }

    result = GenomicPositionRecurrenceAlgorithm().run(
        list(events), list(features), _FAKE_GENE, params
    )

    assert result.Summary["determinable"] is True
    assert result.Summary["reference_build"] == "GRCh37"
    assert result.Summary["n_events_analyzed"] == 2
    assert result.Summary["n_events_excluded_other_build"] == 1
    assert result.Summary["reference_build_counts"] == {"GRCh37": 2, "GRCh38": 1}
    assert any("more than one reference genome build" in warning for warning in result.Warnings)
    # The excluded GRCh38 record must never be pooled into a table alongside
    # GRCh37 positions.
    for row in result.Tables["genomic_bin_recurrence"]:
        assert "c" not in row["event_ids"]


def test_custom_bin_size_is_honored():
    records = [_event_and_feature("a", 439), _event_and_feature("b", 439)]
    events, features = zip(*records)
    params = {
        "genomic_breakpoints": {
            "a": _breakpoint("7", 100_000),
            "b": _breakpoint("7", 100_500),
        },
        "bin_size_bp": 2_000,
    }

    result = GenomicPositionRecurrenceAlgorithm().run(
        list(events), list(features), _FAKE_GENE, params
    )

    assert result.Summary["bin_size_bp"] == 2_000
    bins = result.Tables["genomic_bin_recurrence"]
    assert len(bins) == 1
    assert bins[0]["n_events"] == 2
    assert bins[0]["bin_start"] == 100_000
    assert bins[0]["bin_end"] == 102_000


def test_result_schema_matches_canonical_algorithm_result_fields():
    event, feature = _event_and_feature("a", 439)
    result = GenomicPositionRecurrenceAlgorithm().run(
        [event],
        [feature],
        _FAKE_GENE,
        {"genomic_breakpoints": {"a": _breakpoint("7", 100)}},
    )
    assert set(type(result).model_fields) == set(AlgorithmResult.model_fields)
