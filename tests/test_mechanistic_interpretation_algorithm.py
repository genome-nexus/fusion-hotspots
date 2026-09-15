from __future__ import annotations

import pytest

from cfh.algorithms import registry
from cfh.algorithms.mechanistic_interpretation import MechanisticInterpretationAlgorithm
from cfh.genes.registry import GeneConfig, KeyDomain
from cfh.model.algorithm_result import AlgorithmResult
from cfh.model.fusion_event import FusionEvent
from cfh.model.fusion_feature import FusionFeature
from cfh.orchestrator.run import run_algorithms

_RETENTION_DOMAIN = KeyDomain(name="Kinase domain", source="test", key="kinase")
_DISRUPTION_DOMAIN = KeyDomain(name="Regulatory domain", source="test", key="reg")

_GENE_BOTH_DOMAINS = GeneConfig(
    gene_symbol="FAKE1",
    canonical_transcript_id="NM_000001",
    protein_id="P00001",
    key_domains=[_RETENTION_DOMAIN],
    disruption_required_domains=[_DISRUPTION_DOMAIN],
)

_GENE_NO_DOMAINS = GeneConfig(
    gene_symbol="FAKE2",
    canonical_transcript_id="NM_000002",
    protein_id="P00002",
)


def _event(event_id: str, *, frame_status: str, partner: str) -> FusionEvent:
    return FusionEvent(
        Event_id=event_id,
        Cohort="synthetic",
        Sample_id=f"S-{event_id}",
        Frame_status=frame_status,
        Is_protein_fusion=True,
        Site1_gene="FAKE1",
        Site2_gene=partner,
    )


def _feature(event_id: str, *, domain_key: str, status: str, position: int) -> FusionFeature:
    return FusionFeature(
        Event_id=event_id,
        Gene="FAKE1",
        Junction_position_aa=position,
        Domain_retention_flags={domain_key: status},
    )


def _retention_significant_with_subcluster() -> tuple[list[FusionEvent], list[FusionFeature]]:
    """8 in-frame+retained, 3 in-frame+lost/disrupted (2 sharing partner
    SHARED, 1 with a distinct partner), 5 out-of-frame+lost/disrupted --
    a clean-separation table (OR=inf, p tiny) like
    ``test_domain_disruption_algorithm``'s own fixture, with a
    recurrent-partner counter-intuitive subcluster baked in."""
    events: list[FusionEvent] = []
    features: list[FusionFeature] = []
    for index in range(8):
        event_id = f"r{index}"
        events.append(_event(event_id, frame_status="in-frame", partner=f"P{index}"))
        features.append(_feature(event_id, domain_key="kinase", status="retained", position=100 + index))
    for index, partner in enumerate(["SHARED", "SHARED", "LONER"]):
        event_id = f"c{index}"
        events.append(_event(event_id, frame_status="in-frame", partner=partner))
        features.append(
            _feature(
                event_id,
                domain_key="kinase",
                status="lost" if index != 1 else "disrupted",
                position=200 + index,
            )
        )
    for index in range(5):
        event_id = f"o{index}"
        events.append(_event(event_id, frame_status="out-of-frame", partner=f"Q{index}"))
        features.append(_feature(event_id, domain_key="kinase", status="lost", position=300 + index))
    return events, features


def _retention_significant_no_counter_events() -> tuple[list[FusionEvent], list[FusionFeature]]:
    events: list[FusionEvent] = []
    features: list[FusionFeature] = []
    for index in range(8):
        event_id = f"r{index}"
        events.append(_event(event_id, frame_status="in-frame", partner=f"P{index}"))
        features.append(_feature(event_id, domain_key="kinase", status="retained", position=100 + index))
    for index in range(5):
        event_id = f"o{index}"
        events.append(_event(event_id, frame_status="out-of-frame", partner=f"Q{index}"))
        features.append(_feature(event_id, domain_key="kinase", status="lost", position=300 + index))
    return events, features


def _retention_significant_diffuse_counter_events() -> tuple[list[FusionEvent], list[FusionFeature]]:
    """Same clean-separation shape, but the 2 counter events have distinct
    partners -- no partner recurs, so this must land in
    ``"insufficient_recurrence"`` rather than ``"possible_subcluster"``."""
    events: list[FusionEvent] = []
    features: list[FusionFeature] = []
    for index in range(8):
        event_id = f"r{index}"
        events.append(_event(event_id, frame_status="in-frame", partner=f"P{index}"))
        features.append(_feature(event_id, domain_key="kinase", status="retained", position=100 + index))
    for index, partner in enumerate(["ALONE1", "ALONE2"]):
        event_id = f"c{index}"
        events.append(_event(event_id, frame_status="in-frame", partner=partner))
        features.append(_feature(event_id, domain_key="kinase", status="lost", position=200 + index))
    for index in range(5):
        event_id = f"o{index}"
        events.append(_event(event_id, frame_status="out-of-frame", partner=f"Q{index}"))
        features.append(_feature(event_id, domain_key="kinase", status="lost", position=300 + index))
    return events, features


def _retention_not_significant() -> tuple[list[FusionEvent], list[FusionFeature]]:
    """Retained/lost split evenly between in-frame and out-of-frame --
    no enrichment in either direction."""
    events: list[FusionEvent] = []
    features: list[FusionFeature] = []
    statuses = ["retained", "lost", "retained", "lost"]
    frames = ["in-frame", "in-frame", "out-of-frame", "out-of-frame"]
    for index, (frame_status, status) in enumerate(zip(frames, statuses, strict=True)):
        event_id = f"e{index}"
        events.append(_event(event_id, frame_status=frame_status, partner=f"P{index}"))
        features.append(_feature(event_id, domain_key="kinase", status=status, position=100 + index))
    return events, features


def _run_domain_results(
    events: list[FusionEvent], features: list[FusionFeature], gene_config: GeneConfig
) -> list[AlgorithmResult]:
    """Route through ``run_algorithms`` (not a direct ``.run()`` call) so a
    domain with no matching feature data in this fixture (e.g. the
    disruption domain in a retention-only fixture) degrades to a
    warning-carrying result, the same as in production, rather than
    raising and taking down the whole test."""
    return run_algorithms(
        ["domain_retention", "domain_disruption"],
        events,
        features,
        gene_config,
        {
            "domain_retention": {"seed": 42, "n_permutations": 200},
            "domain_disruption": {"seed": 42, "n_permutations": 200},
        },
    )


def test_algorithm_registered():
    assert registry.get("mechanistic_interpretation") is MechanisticInterpretationAlgorithm
    assert "mechanistic_interpretation" in registry.list_algorithms()


def test_result_schema_matches_canonical_algorithm_result_fields():
    events, features = _retention_significant_no_counter_events()
    upstream = _run_domain_results(events, features, _GENE_BOTH_DOMAINS)
    result = MechanisticInterpretationAlgorithm().run(
        events, features, _GENE_BOTH_DOMAINS, {"algorithm_results": upstream}
    )
    assert set(type(result).model_fields) == set(AlgorithmResult.model_fields)


def test_not_statistically_supported_is_too_weak_to_flag_anything():
    """The 'too weak to associate causation with' branch: when the
    underlying domain_retention finding itself doesn't clear p<0.05/OR>1,
    no counter-intuitive analysis is attempted at all."""
    events, features = _retention_not_significant()
    upstream = _run_domain_results(events, features, _GENE_BOTH_DOMAINS)
    result = MechanisticInterpretationAlgorithm().run(
        events, features, _GENE_BOTH_DOMAINS, {"algorithm_results": upstream}
    )
    assert result.Summary["retention_statistically_supported"] is False
    assert result.Summary["retention_counter_intuitive_confidence"] == "not_applicable"
    assert result.Summary["retention_counter_intuitive_count"] == 0
    assert "retention_counter_intuitive_events" not in (result.Tables or {})


def test_supported_with_no_counter_events_is_confidence_none():
    events, features = _retention_significant_no_counter_events()
    upstream = _run_domain_results(events, features, _GENE_BOTH_DOMAINS)
    result = MechanisticInterpretationAlgorithm().run(
        events, features, _GENE_BOTH_DOMAINS, {"algorithm_results": upstream}
    )
    assert result.Summary["retention_statistically_supported"] is True
    assert result.Summary["retention_fisher_p_value"] < 0.05
    assert result.Summary["retention_counter_intuitive_confidence"] == "none"
    assert result.Summary["retention_counter_intuitive_count"] == 0


def test_diffuse_counter_events_are_insufficient_recurrence():
    events, features = _retention_significant_diffuse_counter_events()
    upstream = _run_domain_results(events, features, _GENE_BOTH_DOMAINS)
    result = MechanisticInterpretationAlgorithm().run(
        events, features, _GENE_BOTH_DOMAINS, {"algorithm_results": upstream}
    )
    assert result.Summary["retention_statistically_supported"] is True
    assert result.Summary["retention_counter_intuitive_count"] == 2
    assert result.Summary["retention_counter_intuitive_confidence"] == "insufficient_recurrence"
    assert "retention_counter_intuitive_recurrent_partners" not in (result.Tables or {})
    events_table = result.Tables["retention_counter_intuitive_events"]
    assert {row["partner_gene"] for row in events_table} == {"ALONE1", "ALONE2"}


def test_recurrent_partner_among_counter_events_is_possible_subcluster():
    events, features = _retention_significant_with_subcluster()
    upstream = _run_domain_results(events, features, _GENE_BOTH_DOMAINS)
    result = MechanisticInterpretationAlgorithm().run(
        events, features, _GENE_BOTH_DOMAINS, {"algorithm_results": upstream}
    )
    assert result.Summary["retention_statistically_supported"] is True
    assert result.Summary["retention_counter_intuitive_count"] == 3
    # Only the 11 in-frame events (8 retained + 3 counter) count toward the
    # denominator -- the 5 out-of-frame events were never part of the
    # in-frame retention claim's population.
    assert result.Summary["retention_counter_intuitive_total"] == 11
    assert result.Summary["retention_counter_intuitive_percent"] == pytest.approx(100.0 * 3 / 11)
    assert result.Summary["retention_counter_intuitive_confidence"] == "possible_subcluster"
    recurrent = result.Tables["retention_counter_intuitive_recurrent_partners"]
    assert recurrent == [{"partner_gene": "SHARED", "count": 2}]
    events_table = result.Tables["retention_counter_intuitive_events"]
    assert {row["event_id"] for row in events_table} == {"c0", "c1", "c2"}
    assert all(row["domain_status"] in {"lost", "disrupted"} for row in events_table)


def test_disruption_effect_mirrors_retention_with_opposite_contradicting_status():
    """Same shape, but for a disruption_required_domain: the 'counter'
    status is 'retained' rather than 'lost'/'disrupted'."""
    events: list[FusionEvent] = []
    features: list[FusionFeature] = []
    for index in range(8):
        event_id = f"d{index}"
        events.append(_event(event_id, frame_status="in-frame", partner=f"P{index}"))
        features.append(_feature(event_id, domain_key="reg", status="lost", position=100 + index))
    for index, partner in enumerate(["SHARED", "SHARED"]):
        event_id = f"c{index}"
        events.append(_event(event_id, frame_status="in-frame", partner=partner))
        features.append(_feature(event_id, domain_key="reg", status="retained", position=200 + index))
    for index in range(5):
        event_id = f"o{index}"
        events.append(_event(event_id, frame_status="out-of-frame", partner=f"Q{index}"))
        features.append(_feature(event_id, domain_key="reg", status="retained", position=300 + index))

    upstream = _run_domain_results(events, features, _GENE_BOTH_DOMAINS)
    result = MechanisticInterpretationAlgorithm().run(
        events, features, _GENE_BOTH_DOMAINS, {"algorithm_results": upstream}
    )
    assert result.Summary["disruption_statistically_supported"] is True
    assert result.Summary["disruption_counter_intuitive_confidence"] == "possible_subcluster"
    recurrent = result.Tables["disruption_counter_intuitive_recurrent_partners"]
    assert recurrent == [{"partner_gene": "SHARED", "count": 2}]
    events_table = result.Tables["disruption_counter_intuitive_events"]
    assert all(row["domain_status"] == "retained" for row in events_table)


def test_no_domains_configured_is_not_applicable_for_both_effects():
    result = MechanisticInterpretationAlgorithm().run([], [], _GENE_NO_DOMAINS, {})
    assert result.Summary["retention_domains_configured"] is False
    assert result.Summary["retention_counter_intuitive_confidence"] == "not_applicable"
    assert result.Summary["disruption_domains_configured"] is False
    assert result.Summary["disruption_counter_intuitive_confidence"] == "not_applicable"
    assert result.Tables == {}


def test_no_op_for_gene_pair_configs():
    config = GeneConfig(gene_pair=("TMPRSS2", "ERG"), analysis_modes=["promoter_swap"])
    result = MechanisticInterpretationAlgorithm().run([], [], config, {})
    assert result.Summary == {}
    assert result.Tables == {}
    assert result.Warnings
    assert "gene_pair" in result.Warnings[0]


def test_missing_upstream_results_is_treated_as_not_supported_not_a_crash():
    """A caller that requests mechanistic_interpretation without also
    requesting domain_retention/domain_disruption gets a graceful
    not_applicable verdict, never a KeyError."""
    events, features = _retention_significant_with_subcluster()
    result = MechanisticInterpretationAlgorithm().run(events, features, _GENE_BOTH_DOMAINS, {})
    assert result.Summary["retention_counter_intuitive_confidence"] == "not_applicable"
    assert result.Summary["disruption_counter_intuitive_confidence"] == "not_applicable"


def test_wired_through_orchestrator_dependency_injection():
    """End-to-end through run_algorithms: DEPENDS_ON auto-injects the
    already-computed domain_retention/domain_disruption results, the same
    mechanism composite_score relies on."""
    events, features = _retention_significant_with_subcluster()
    results = run_algorithms(
        ["domain_retention", "domain_disruption", "mechanistic_interpretation"],
        events,
        features,
        _GENE_BOTH_DOMAINS,
        {
            "domain_retention": {"seed": 42, "n_permutations": 200},
            "domain_disruption": {"seed": 42, "n_permutations": 200},
        },
    )
    by_name = {result.Algorithm: result for result in results}
    mechanistic = by_name["mechanistic_interpretation"]
    assert not mechanistic.Warnings
    assert mechanistic.Summary["retention_counter_intuitive_confidence"] == "possible_subcluster"
