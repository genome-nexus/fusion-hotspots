"""Junction concentration: breakpoints redrawn along the gene as the null."""

from __future__ import annotations

from unittest.mock import MagicMock

import numpy as np
import pytest

from cfh.algorithms.junction_concentration import JunctionConcentrationAlgorithm
from cfh.gene_comparison import collect_p_values_from_algorithm_results
from cfh.genes.registry import GeneConfig
from cfh.mapping.genome_nexus_source import ExonRecord
from cfh.model.fusion_event import FusionEvent
from cfh.model.fusion_feature import FusionFeature
from cfh.stats.junction_concentration import (
    PANEL_CAPTURE_CAVEAT,
    junction_concentration_test,
    junction_probabilities,
)

# Four 100 bp exons separated by 900, 2,900 and 900 bp introns: 5,100 bp.
_EXONS = [
    (1, 1_000, 1_099),
    (2, 2_000, 2_099),
    (3, 5_000, 5_099),
    (4, 6_000, 6_099),
]


def _exon_records(strand: int = 1) -> list[ExonRecord]:
    return [ExonRecord(f"E{rank}", start, end, rank, strand) for rank, start, end in _EXONS]


def _payload() -> dict:
    return {
        "transcriptId": "ENST_TEST",
        "exons": [
            {"exonId": f"E{rank}", "exonStart": start, "exonEnd": end, "rank": rank, "strand": 1}
            for rank, start, end in _EXONS
        ],
        "utrs": [],
    }


def test_three_prime_junction_is_the_first_exon_after_a_random_break():
    probabilities = junction_probabilities(_exon_records(), "three_prime")

    # Exon 1 itself, then each intron with the exon it leads into.
    assert probabilities == pytest.approx(
        {1: 100 / 5100, 2: 1000 / 5100, 3: 3000 / 5100, 4: 1000 / 5100}
    )


def test_five_prime_junction_is_the_last_exon_before_a_random_break():
    probabilities = junction_probabilities(_exon_records(), "five_prime")

    assert probabilities == pytest.approx(
        {1: 1000 / 5100, 2: 3000 / 5100, 3: 1000 / 5100, 4: 100 / 5100}
    )


def test_a_pile_up_at_a_short_intron_exon_is_detected():
    probabilities = {"three_prime": junction_probabilities(_exon_records(), "three_prime")}
    # 30 of 40 fusions at exon 4, which random breakage reaches 20% of the time.
    junctions = [("three_prime", 4)] * 30 + [("three_prime", 3)] * 10

    result = junction_concentration_test(junctions, probabilities, n_permutations=999)

    assert result["top_junction"]["junction_exon"] == 4
    assert result["top_junction"]["observed_share"] == pytest.approx(0.75)
    assert result["top_junction"]["random_breakage_share"] == pytest.approx(1000 / 5100)
    assert result["permutation_empirical_p_value"] == pytest.approx(1 / 1000)


def test_the_largest_exon_is_not_significant_just_for_being_the_most_common():
    probabilities = {"three_prime": junction_probabilities(_exon_records(), "three_prime")}
    # Exon 3 is the modal junction, but at its random-breakage rate.
    junctions = [("three_prime", 3)] * 24 + [("three_prime", 2)] * 8 + [("three_prime", 4)] * 8

    result = junction_concentration_test(junctions, probabilities, n_permutations=999)

    assert result["permutation_empirical_p_value"] > 0.2


def test_random_breakpoints_reject_at_about_the_nominal_rate():
    probabilities = junction_probabilities(_exon_records(), "three_prime")
    rng = np.random.default_rng(3)
    exons, weights = list(probabilities), list(probabilities.values())
    rejections = 0
    for replicate in range(200):
        draw = rng.choice(exons, size=30, p=weights)
        result = junction_concentration_test(
            [("three_prime", int(exon)) for exon in draw],
            {"three_prime": probabilities},
            seed=replicate,
            n_permutations=199,
        )
        rejections += result["permutation_empirical_p_value"] <= 0.05
    assert rejections / 200 <= 0.09


def _event(event_id: str) -> FusionEvent:
    return FusionEvent(
        Event_id=event_id,
        Cohort="c",
        Patient_id=event_id,
        Site1_gene="PARTNER",
        Site2_gene="GENE",
        Frame_status="in-frame",
    )


def _feature(event_id: str, exon: int) -> FusionFeature:
    return FusionFeature(Event_id=event_id, Gene="GENE", Role="three_prime", Breakpoint_exon=exon)


def _config() -> GeneConfig:
    return GeneConfig(gene_symbol="GENE", canonical_transcript_id="ENST_TEST", protein_id="P1")


def _client() -> MagicMock:
    client = MagicMock()
    client.fetch_canonical_transcript.return_value = _payload()
    return client


def test_algorithm_tests_inside_the_transcript_and_reports_the_caveat():
    exons = [4] * 12 + [3] * 3 + [1] * 5
    events = [_event(f"E{i}") for i in range(len(exons))]
    features = [_feature(f"E{i}", exon) for i, exon in enumerate(exons)]
    # The five exon-1 calls sit upstream of the transcript, clamped by the mapper.
    breakpoints = {
        f"E{i}": {"position": 900 if exon == 1 else 5_500} for i, exon in enumerate(exons)
    }

    result = JunctionConcentrationAlgorithm().run(
        events,
        features,
        _config(),
        {
            "genome_nexus_client": _client(),
            "genomic_breakpoints": breakpoints,
            "n_permutations": 999,
        },
    )

    summary = result.Summary
    assert summary["determinable"] is True
    assert summary["n_events_analyzed"] == 15
    assert summary["n_outside_transcript"] == 5
    assert summary["top_junction"]["junction_exon"] == 4
    assert summary["permutation_empirical_p_value"] < 0.01
    assert summary["panel_capture_caveat"] == PANEL_CAPTURE_CAVEAT
    # Corrected by the cohort scan as its own family, not the main one.
    assert collect_p_values_from_algorithm_results("GENE", "study", [result], source="t") == []


def test_algorithm_is_not_determinable_without_a_client_or_enough_events():
    events = [_event(f"E{i}") for i in range(3)]
    features = [_feature(f"E{i}", 4) for i in range(3)]

    without_client = JunctionConcentrationAlgorithm().run(events, features, _config(), {})
    too_few = JunctionConcentrationAlgorithm().run(
        events, features, _config(), {"genome_nexus_client": _client()}
    )

    assert without_client.Summary["determinable"] is False
    assert "Genome Nexus" in without_client.Summary["reason"]
    assert too_few.Summary["determinable"] is False
    assert too_few.Summary["permutation_empirical_p_value"] is None
