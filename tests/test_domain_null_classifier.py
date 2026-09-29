"""The domain-status null classifies simulated breakpoints from coordinates."""

from __future__ import annotations

from cfh.genes.registry import load_gene_config
from cfh.model.fusion_event import FusionEvent
from cfh.model.fusion_feature import DomainRetentionDetail, FusionFeature
from cfh.stats.breakpoint_tests import permutation_null_classifier, permutation_null_test

CONFIG = load_gene_config("BRAF")
KEY = CONFIG.key_domains[0].key or CONFIG.key_domains[0].name


def _record(event_id, *, frame, role, status, position=300, details=True):
    event = FusionEvent(Event_id=event_id, Cohort="c", Is_protein_fusion=True, Frame_status=frame)
    feature = FusionFeature(
        Event_id=event_id,
        Gene="BRAF",
        Role=role,
        Junction_position_aa=position,
        Domain_retention_flags={KEY: status},
        Domain_retention_details=(
            {KEY: DomainRetentionDetail(Domain_start_aa=400, Domain_end_aa=600)}
            if details
            else None
        ),
    )
    return event, feature


def _split(records):
    return [event for event, _ in records], [feature for _, feature in records]


def test_null_uses_coordinates_and_observed_roles_not_nearest_label():
    # A 5' partner breaking at aa 300 loses a 400-600 domain; a 3' partner
    # breaking at the same residue retains it. The first observation at
    # aa 300 is the 5' "lost" one, so the old nearest-label rule would call
    # every simulated in-frame breakpoint "lost".
    records = [
        _record("OUT", frame="out-of-frame", role="five_prime", status="lost"),
        _record("IN1", frame="in-frame", role="three_prime", status="retained"),
        _record("IN2", frame="in-frame", role="three_prime", status="retained"),
    ]
    events, features = _split(records)
    assert permutation_null_classifier(features, CONFIG) == "domain_coordinates"
    _, observed, null_rates = permutation_null_test(
        events, features, CONFIG, seed=1, n_permutations=20
    )
    assert observed == 1.0
    assert set(null_rates) == {1.0}


def test_missing_coordinates_fall_back_to_nearest_label():
    records = [
        _record("IN1", frame="in-frame", role="three_prime", status="retained", details=False),
        _record("IN2", frame="in-frame", role=None, status="lost", position=700),
    ]
    _, features = _split(records)
    assert permutation_null_classifier(features, CONFIG) == "nearest_observed_label"


def test_retention_result_reports_its_null_classifier():
    from cfh.algorithms.domain_retention import DomainRetentionAlgorithm

    records = [
        _record("IN1", frame="in-frame", role="three_prime", status="retained"),
        _record("OUT1", frame="out-of-frame", role="three_prime", status="lost", position=500),
    ]
    events, features = _split(records)
    result = DomainRetentionAlgorithm().run(events, features, CONFIG, {"n_permutations": 5})
    assert result.Summary["permutation_null_classifier"] == "domain_coordinates"
