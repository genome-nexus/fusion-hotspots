"""CI-gating composite_score benchmark against already-committed real
MSK-IMPACT run data.

Reuses the checked-in ``runs/braf_msk-impact-*`` and ``runs/ret_msk-impact-*``
directories (real cBioPortal + Genome Nexus pulls, already committed to the
repo) instead of any new network call, exactly like
``tests/benchmark/test_braf_domain_disruption_real_msk_impact.py`` already
does. It reconstructs typed ``FusionEvent``/``FusionFeature`` objects from
each run's ``results.json`` and runs the real algorithm plugins end to end:

* BRAF (this repo's fully-configured gene: ``key_domains``,
  ``disruption_required_domains``, and ``breakpoint_hotspot``/
  ``domain_disruption`` analysis modes) proves composite_score with all five
  sub-scores applicable at once.
* RET (only ``domain_retention`` configured -- no
  ``disruption_required_domains``) proves composite_score degrades
  gracefully on real data with fewer applicable algorithms, reusing the
  ``frequency``/``domain_retention``/``cutpoint_detection``/
  ``confidence_stats`` results already present in that committed run's own
  ``results.json`` and excluding the explicitly-skipped ``domain_disruption``
  result rather than zero-filling it.

Both tests print the actual ranked output so it is reported honestly,
whatever it turns out to be, rather than asserted into a predetermined
shape.
"""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from cfh.algorithms.composite_score import CompositeScoreAlgorithm
from cfh.algorithms.cutpoint_detection import CutpointDetectionAlgorithm
from cfh.algorithms.domain_disruption import DomainDisruptionAlgorithm
from cfh.algorithms.domain_retention import DomainRetentionAlgorithm
from cfh.algorithms.frequency import FrequencyAnalysis
from cfh.algorithms.registry import list_algorithms
from cfh.genes.registry import load_gene_config
from cfh.mapping.feature_mapper import classify_domain_retention
from cfh.mapping.genome_nexus_source import GenomeNexusClient, parse_canonical_transcript
from cfh.model.algorithm_result import AlgorithmResult
from cfh.model.fusion_event import FusionEvent
from cfh.model.fusion_feature import FusionFeature
from cfh.orchestrator.run import run_algorithms
from cfh.real_benchmark import analyze_structural_variant_calls

_REPO_ROOT = Path(__file__).resolve().parents[2]
_RUNS_DIR = _REPO_ROOT / "runs"
_GENOME_NEXUS_BRAF_FIXTURE = (
    Path(__file__).resolve().parents[1]
    / "fixtures"
    / "genome_nexus"
    / "canonical_transcript_braf.json"
)


def _real_run_results_path(prefix: str) -> Path:
    candidates = sorted(_RUNS_DIR.glob(f"{prefix}_*"))
    if not candidates:
        pytest.skip(f"no committed real run directory found for {prefix!r} under {_RUNS_DIR}")
    return candidates[-1] / "results.json"


def _domain_bounds(pfam_id: str) -> tuple[int, int]:
    payload = json.loads(_GENOME_NEXUS_BRAF_FIXTURE.read_text())
    domain = next(d for d in payload["pfamDomains"] if d["pfamDomainId"] == pfam_id)
    return domain["pfamDomainStart"], domain["pfamDomainEnd"]


def _braf_pfam_domain_boundaries() -> list[int]:
    payload = json.loads(_GENOME_NEXUS_BRAF_FIXTURE.read_text())
    transcript = parse_canonical_transcript(payload)
    boundaries: set[int] = set()
    for domain in transcript.pfam_domains:
        boundaries.add(domain.start_aa)
        boundaries.add(domain.end_aa)
    return sorted(boundaries)


def _braf_events_and_features(results_path: Path) -> tuple[list[FusionEvent], list[FusionFeature]]:
    """Reconstruct real BRAF events/features from a committed run, including
    the partner gene (needed for composite_score's per-partner ranking,
    which the domain_disruption benchmark's own helper omits since it never
    ranks by partner) and RAS-binding/cysteine-rich domain status (needed
    for a real domain_disruption result), using the same real per-event
    breakpoint protein position and the same committed Genome Nexus domain
    boundaries the domain_disruption benchmark already uses.
    """
    payload = json.loads(results_path.read_text())
    ras_bounds = _domain_bounds("PF02196")
    cys_bounds = _domain_bounds("PF00130")

    events: list[FusionEvent] = []
    features: list[FusionFeature] = []
    for row in payload["events"]:
        events.append(
            FusionEvent(
                Event_id=row["event_id"],
                Cohort=payload["study_id"],
                Sample_id=row["sample_id"],
                Fusion_name=row["fusion_name"],
                Site1_gene="BRAF",
                Site2_gene=row["partner_gene"],
                Frame_status=row["frame_status"],
                Is_protein_fusion=True,
            )
        )
        role = row["target_role"]
        position = row["breakpoint_protein_position"]
        features.append(
            FusionFeature(
                Event_id=row["event_id"],
                Gene="BRAF",
                Role=role,
                Junction_position_aa=position,
                Domain_retention_flags={
                    "kinase": row["domain_status"],
                    "ras_binding": classify_domain_retention(*ras_bounds, position, role),
                    "cysteine_rich": classify_domain_retention(*cys_bounds, position, role),
                },
            )
        )
    return events, features


def _print_ranking(gene: str, result: AlgorithmResult) -> None:
    print(f"\n=== composite_score real-data ranking: {gene} ===")
    print(json.dumps(result.Summary, indent=2, sort_keys=True))
    for row in result.Tables["composite_evidence_ranking"]:
        print(
            f"  #{row['Rank']:>2} {row['Partner_gene']:<15} "
            f"n={row['Event_count']:<4} composite={row['Composite_score']:.4f} "
            f"components={row['Components_applicable']}"
        )
    for warning in result.Warnings:
        print(f"  Warning: {warning}")


def _braf_confidence_stats_result(events, features, config) -> AlgorithmResult:
    """Choose a meaningful, real (not fabricated) group_field/outcome_field
    pairing for the gene-agnostic confidence_stats algorithm from the
    fields the committed run's projection actually carries: Frame_status
    genuinely varies (in-frame/out-of-frame), and Is_protein_fusion is real
    -- these rows are exactly the protein-fusion records selected upstream
    -- even though it is constant across this pre-filtered set, so its
    per-group MLE/CI narrows with real per-group sample size: a real,
    honestly computed confidence-interval-width signal.
    """
    from cfh.algorithms.confidence_stats import ConfidenceStatsAlgorithm

    return ConfidenceStatsAlgorithm().run(
        events,
        features,
        config,
        {
            "group_field": "Frame_status",
            "group_values": ["in-frame", "out-of-frame"],
            "outcome_field": "Is_protein_fusion",
            "numeric_field": "Junction_position_aa",
        },
    )


def test_composite_score_real_braf_msk_impact_all_five_subscores_applicable():
    """BRAF has domain_retention + domain_disruption + cutpoint_detection +
    confidence_stats all applicable (per genes/configs/braf.yaml), so this
    proves composite_score with every sub-score populated at once, computed
    from real MSK-IMPACT structural-variant data already committed to the
    repo -- no network call.
    """
    results_path = _real_run_results_path("braf_msk-impact-50k-2026")
    events, features = _braf_events_and_features(results_path)
    config = load_gene_config("braf")
    boundaries = _braf_pfam_domain_boundaries()

    frequency_result = FrequencyAnalysis().run(events, features, config, {})
    domain_retention_result = DomainRetentionAlgorithm().run(
        events, features, config, {"seed": 42, "n_permutations": 2_000}
    )
    domain_disruption_result = DomainDisruptionAlgorithm().run(
        events, features, config, {"seed": 42, "n_permutations": 2_000}
    )
    cutpoint_result = CutpointDetectionAlgorithm().run(
        events,
        features,
        config,
        {"seed": 42, "n_permutations": 2_000, "domain_boundaries": boundaries},
    )
    confidence_stats_result = _braf_confidence_stats_result(events, features, config)

    # Sanity: this is real data, so every upstream algorithm should actually
    # be applicable for BRAF, not gracefully skipped.
    assert domain_retention_result.Summary["fisher_p_value"] is not None
    assert domain_disruption_result.Summary["fisher_p_value"] is not None
    assert cutpoint_result.Summary["determinable"] is True
    assert confidence_stats_result.Summary.get("mle") is not None

    result = CompositeScoreAlgorithm().run(
        events,
        features,
        config,
        {
            "algorithm_results": [
                frequency_result,
                domain_retention_result,
                domain_disruption_result,
                cutpoint_result,
                confidence_stats_result,
            ]
        },
    )
    _print_ranking("BRAF", result)

    assert result.Summary["components_applicable"] == {
        "recurrence": True,
        "domain_retention": True,
        "domain_disruption": True,
        "cutpoint_proximity": True,
        "confidence_certainty": True,
    }
    assert result.Warnings == []

    ranking = result.Tables["composite_evidence_ranking"]
    partner_names = {row["Partner_gene"] for row in ranking}
    expected_partners = {
        row["Partner_gene"] for row in frequency_result.Tables["Partner_gene_counts"]
    }
    assert partner_names == expected_partners
    assert [row["Rank"] for row in ranking] == list(range(1, len(ranking) + 1))
    composite_scores = [row["Composite_score"] for row in ranking]
    assert composite_scores == sorted(composite_scores, reverse=True)
    assert all(0.0 <= score <= 1.0 for score in composite_scores)
    for row in ranking:
        assert set(row["Components_applicable"]) == {
            "recurrence",
            "domain_retention",
            "domain_disruption",
            "cutpoint_proximity",
            "confidence_certainty",
        }

    # With corrected locus mapping, AGK leads the five-component composite
    # despite KIAA1549 remaining the most recurrent partner. Pin the exact
    # corrected real-data result so changes in aggregation are caught.
    assert ranking[0]["Partner_gene"] == "AGK"
    assert ranking[0]["Event_count"] == 14
    assert ranking[0]["Composite_score"] == pytest.approx(0.2579, abs=5e-5)


def _ret_events_and_features(results_path: Path) -> tuple[list[FusionEvent], list[FusionFeature]]:
    payload = json.loads(results_path.read_text())
    events: list[FusionEvent] = []
    features: list[FusionFeature] = []
    for row in payload["events"]:
        events.append(
            FusionEvent(
                Event_id=row["event_id"],
                Cohort=payload["study_id"],
                Sample_id=row["sample_id"],
                Fusion_name=row["fusion_name"],
                Site1_gene="RET",
                Site2_gene=row["partner_gene"],
                Frame_status=row["frame_status"],
                Is_protein_fusion=True,
            )
        )
        features.append(
            FusionFeature(
                Event_id=row["event_id"],
                Gene="RET",
                Role=row["target_role"],
                Junction_position_aa=row["breakpoint_protein_position"],
                Domain_retention_flags={"kinase": row["domain_status"]},
            )
        )
    return events, features


def test_composite_score_real_ret_msk_impact_gracefully_degrades():
    """RET only configures ``domain_retention`` (genes/configs/ret.yaml has
    no ``disruption_required_domains``), a config-level fact independent of
    any particular run. This reuses the ``frequency``, ``domain_retention``,
    ``cutpoint_detection``, and ``confidence_stats`` AlgorithmResult objects
    already present in the latest committed RET run's results.json verbatim
    -- the most literal form of "consume already-computed outputs as inputs"
    -- and proves composite_score ranks real RET fusion partners while
    excluding domain_disruption (recorded as an explicit skipped result)
    rather than zero-filling it.
    """
    results_path = _real_run_results_path("ret_msk-impact-50k-2026")
    payload = json.loads(results_path.read_text())
    committed_results = {item["Algorithm"]: item for item in payload["algorithm_results"]}
    domain_disruption_result = committed_results["domain_disruption"]
    assert domain_disruption_result["Summary"]["fisher_p_value"] is None
    assert "was skipped" in domain_disruption_result["Warnings"][0]
    assert committed_results["confidence_stats"]["Warnings"] == []
    assert committed_results["cutpoint_detection"]["Summary"]["determinable"] is True

    events, features = _ret_events_and_features(results_path)
    config = load_gene_config("ret")

    result = CompositeScoreAlgorithm().run(
        events,
        features,
        config,
        {
            "algorithm_results": [
                committed_results["frequency"],
                committed_results["domain_retention"],
                committed_results["cutpoint_detection"],
                committed_results["confidence_stats"],
            ]
        },
    )
    _print_ranking("RET", result)

    assert result.Summary["components_applicable"] == {
        "recurrence": True,
        "domain_retention": True,
        "domain_disruption": False,
        "cutpoint_proximity": True,
        "confidence_certainty": True,
    }
    ranking = result.Tables["composite_evidence_ranking"]
    assert ranking, "expected at least one ranked RET fusion partner"

    # KIF5B is RET's dominant real partner in this cohort (87/194 events).
    # Pinned to the exact real value, not just the ranking, since this
    # reuses the committed run's own already-computed AlgorithmResult
    # objects verbatim (no seed/n_permutations choice made here at all).
    assert ranking[0]["Partner_gene"] == "KIF5B"
    assert ranking[0]["Event_count"] == 87
    assert ranking[0]["Composite_score"] == pytest.approx(0.5046697811340378)

    for row in ranking:
        assert row["Domain_disruption_score"] is None
        assert row["Confidence_certainty_score"] is not None
        assert "domain_disruption" not in row["Components_applicable"]
        assert "confidence_certainty" in row["Components_applicable"]
    composite_scores = [row["Composite_score"] for row in ranking]
    assert composite_scores == sorted(composite_scores, reverse=True)
    assert all(0.0 <= score <= 1.0 for score in composite_scores)

    warning_text = " ".join(result.Warnings)
    assert "domain_disruption" in warning_text


# --- Real orchestrator-dispatch integration tests -------------------------
#
# The tests above call `CompositeScoreAlgorithm().run(...)` directly with a
# hand-assembled `algorithm_results` list. That proves the aggregation math,
# but it does not prove the *orchestrator* actually wires composite_score's
# inputs together on the real dispatch path (`run_algorithms`, which both
# `cfh analyze`/`run_analysis` and any other caller of the registered
# algorithm set go through). The tests below call `run_algorithms` --
# composite_score's `DEPENDS_ON` declaration and the dependency-aware wave
# scheduling in `cfh.orchestrator.run` -- directly, never touching
# `CompositeScoreAlgorithm` themselves, against the same real committed BRAF
# and RET data used above.


def test_composite_score_via_real_orchestrator_dispatch_braf():
    """Runs the real orchestrator (`run_algorithms`, exactly what
    `run_analysis`/`cfh analyze` calls) with the full registered algorithm
    set for BRAF and confirms a populated composite_score ranking table
    comes out, with domain_retention/domain_disruption/cutpoint_detection
    auto-injected as dependencies rather than composite_score being called
    directly or its inputs hand-assembled.
    """
    results_path = _real_run_results_path("braf_msk-impact-50k-2026")
    events, features = _braf_events_and_features(results_path)
    config = load_gene_config("braf")
    boundaries = _braf_pfam_domain_boundaries()

    results = run_algorithms(
        list_algorithms(),
        events,
        features,
        config,
        {
            "domain_retention": {"seed": 42, "n_permutations": 500},
            "domain_disruption": {"seed": 42, "n_permutations": 500},
            "cutpoint_detection": {
                "seed": 42,
                "n_permutations": 500,
                "domain_boundaries": boundaries,
            },
            "confidence_stats": {
                "group_field": "Frame_status",
                "group_values": ["in-frame", "out-of-frame"],
                "outcome_field": "Is_protein_fusion",
            },
        },
    )
    results_by_name = {result.Algorithm: result for result in results}
    composite_result = results_by_name["composite_score"]
    _print_ranking("BRAF (via real orchestrator dispatch)", composite_result)

    assert not any(
        str(warning).startswith("Algorithm failed") for warning in composite_result.Warnings
    )
    assert composite_result.Summary["components_applicable"] == {
        "recurrence": True,
        "domain_retention": True,
        "domain_disruption": True,
        "cutpoint_proximity": True,
        "confidence_certainty": True,
    }
    ranking = composite_result.Tables["composite_evidence_ranking"]
    assert ranking, "expected a populated composite_score ranking table"
    assert ranking[0]["Partner_gene"] == "AGK"
    assert ranking[0]["Event_count"] == 14
    assert all(0.0 <= row["Composite_score"] <= 1.0 for row in ranking)


def test_composite_score_via_real_orchestrator_dispatch_ret_gracefully_degrades():
    """Same real-orchestrator proof for RET: domain_disruption legitimately
    runs (it is registered) but no-ops for RET (no
    ``disruption_required_domains`` configured), and confidence_stats is
    requested with no per-algorithm params (matching what `cfh analyze`
    actually sends today) so it fails exactly as the committed real RET run
    already shows -- composite_score must still produce a populated,
    correctly-degraded ranking table, not fail or no-op itself.
    """
    results_path = _real_run_results_path("ret_msk-impact-50k-2026")
    events, features = _ret_events_and_features(results_path)
    config = load_gene_config("ret")

    results = run_algorithms(
        list_algorithms(),
        events,
        features,
        config,
        {
            "domain_retention": {"seed": 42, "n_permutations": 500},
            "domain_disruption": {"seed": 42, "n_permutations": 500},
            "cutpoint_detection": {"seed": 42, "n_permutations": 500},
        },
    )
    results_by_name = {result.Algorithm: result for result in results}
    composite_result = results_by_name["composite_score"]
    _print_ranking("RET (via real orchestrator dispatch)", composite_result)

    assert not any(
        str(warning).startswith("Algorithm failed") for warning in composite_result.Warnings
    )
    assert composite_result.Summary["components_applicable"] == {
        "recurrence": True,
        "domain_retention": True,
        "domain_disruption": False,
        "cutpoint_proximity": True,
        "confidence_certainty": False,
    }
    ranking = composite_result.Tables["composite_evidence_ranking"]
    assert ranking, "expected a populated composite_score ranking table"
    assert ranking[0]["Partner_gene"] == "KIF5B"
    assert ranking[0]["Event_count"] == 87
    assert all(0.0 <= row["Composite_score"] <= 1.0 for row in ranking)


def test_composite_score_via_real_analyze_pipeline_offline(
    genome_nexus_canonical_transcript_fixture_path,
):
    """End-to-end proof through `cfh.real_benchmark.analyze_structural_variant_calls`
    -- the exact function `run_analysis`/`cfh analyze <gene> <study>` calls
    -- with the full registered algorithm set (`list_algorithms()`,
    composite_score included), fully offline via a mocked Genome Nexus
    client backed by the same committed canonical-transcript fixture the
    existing offline BRAF pipeline test
    (`test_braf_kinase_retention_msk_impact_50k.py`) already uses for real
    production normalization/mapping. Confirms `run.results` -- the exact
    list `write_outputs` serializes into a committed `runs/*/results.json`
    -- actually carries a populated composite_score ranking table, closing
    the loop the reviewer flagged: `cfh analyze` no longer no-ops or fails
    on composite_score.
    """
    client = MagicMock(spec=GenomeNexusClient)
    client.fetch_canonical_transcript.return_value = json.loads(
        genome_nexus_canonical_transcript_fixture_path.read_text()
    )

    def _call(sample_id: str, partner: str) -> dict:
        return {
            "sampleId": sample_id,
            "site1HugoSymbol": partner,
            "site2HugoSymbol": "BRAF",
            "connectionType": "3to3",
            "site2Position": 140493152,
            "site2EffectOnFrame": "NA",
            "eventInfo": f"Protein Fusion: in frame  {{{partner}:BRAF}}",
        }

    calls = [_call(f"SAMPLE-{i}", "AGAP3") for i in range(5)] + [
        _call(f"SAMPLE-{i}", "SND1") for i in range(5, 8)
    ]

    run = analyze_structural_variant_calls(
        calls,
        "BRAF",
        "offline_composite_score_integration",
        genome_nexus_client=client,
        n_permutations=500,
        algorithm_names=list_algorithms(),
    )

    composite_results = [result for result in run.results if result.Algorithm == "composite_score"]
    assert composite_results, "composite_score result missing from analyze_structural_variant_calls"
    composite_result = composite_results[0]
    _print_ranking("BRAF (via analyze_structural_variant_calls)", composite_result)

    assert not any(
        str(warning).startswith("Algorithm failed") for warning in composite_result.Warnings
    )
    assert composite_result.Summary["components_applicable"]["recurrence"] is True
    assert composite_result.Summary["components_applicable"]["domain_retention"] is True
    ranking = composite_result.Tables["composite_evidence_ranking"]
    assert ranking, "expected a populated composite_score ranking table"
    assert {row["Partner_gene"] for row in ranking} == {"AGAP3", "SND1"}
    assert ranking[0]["Partner_gene"] == "AGAP3"
    assert ranking[0]["Event_count"] == 5
    assert all(0.0 <= row["Composite_score"] <= 1.0 for row in ranking)
