import csv
import json
from unittest.mock import MagicMock

import pytest

from cfh.algorithms.domain_retention import DomainRetentionAlgorithm
from cfh.genes.registry import load_gene_config
from cfh.mapping.domain_source import ProteinDomain
from cfh.mapping.feature_mapper import calculate_domain_retention, map_event
from cfh.model.fusion_event import FusionEvent
from cfh.model.fusion_feature import FusionFeature
from conftest import RUNS_DIR, latest_run_dir


@pytest.mark.parametrize(
    ("role", "breakpoint", "interval", "fraction", "truncated"),
    [
        ("five_prime", 457, (None, None), 0.0, False),
        ("five_prime", 458, (458, 458), 1 / 255, True),
        ("five_prime", 600, (458, 600), 143 / 255, True),
        ("five_prime", 712, (458, 712), 1.0, False),
        ("five_prime", 800, (458, 712), 1.0, False),
        ("three_prime", 400, (458, 712), 1.0, False),
        ("three_prime", 458, (458, 712), 1.0, False),
        ("three_prime", 600, (600, 712), 113 / 255, True),
        ("three_prime", 712, (712, 712), 1 / 255, True),
        ("three_prime", 713, (None, None), 0.0, False),
    ],
)
def test_braf_pf07714_retained_interval_fraction_and_truncation(
    role, breakpoint, interval, fraction, truncated
):
    detail = calculate_domain_retention(458, 712, breakpoint, role)

    assert (detail.Retained_start_aa, detail.Retained_end_aa) == interval
    assert detail.Retained_fraction == pytest.approx(fraction)
    assert detail.Is_truncated is truncated


def test_mapping_adds_quantitative_detail_without_changing_binary_calls():
    config = load_gene_config("braf")
    source = MagicMock()
    source.fetch.return_value = [
        ProteinDomain(
            name="PF07714",
            accession="PF07714",
            start_aa=458,
            end_aa=712,
            source="genome_nexus",
        )
    ]
    event = FusionEvent(Event_id="braf-truncated", Cohort="test")

    feature = map_event(
        event,
        config,
        role="five_prime",
        junction_position_aa=600,
        domain_source=source,
    )

    assert feature.Domain_retention_flags["kinase"] == "disrupted"
    assert feature.Retained_domains == []
    assert feature.Lost_domains == []
    assert feature.Disrupted_domains == ["Protein kinase domain"]
    detail = feature.Domain_retention_details["kinase"]
    assert (detail.Retained_start_aa, detail.Retained_end_aa) == (458, 600)
    assert detail.Retained_fraction == pytest.approx(143 / 255)
    assert detail.Is_truncated is True


@pytest.mark.parametrize(
    ("run_prefix", "gene", "bounds", "expected"),
    [
        (
            "braf_msk-impact-50k-2026",
            "BRAF",
            (458, 712),
            (178, 179, 163, 91.06145251396649, [[142, 21], [9, 6]], 0.013367557978153668),
        ),
        (
            "ret_msk-impact-50k-2026",
            "RET",
            (724, 1005),
            (194, 194, 179, 92.26804123711341, [[141, 38], [5, 10]], 0.00041966557652448966),
        ),
    ],
)
def test_committed_benchmark_binary_results_are_unchanged_with_quantitative_details(
    run_prefix, gene, bounds, expected
):
    """Sanity-check reconstructed features against the committed live artifact rows."""
    payload = json.loads((latest_run_dir(run_prefix) / "results.json").read_text())
    events = []
    features = []
    for row in payload["events"]:
        events.append(
            FusionEvent(
                Event_id=row["event_id"],
                Cohort=payload["study_id"],
                Frame_status=row["frame_status"],
                Is_protein_fusion=True,
            )
        )
        detail = calculate_domain_retention(
            *bounds, row["breakpoint_protein_position"], row["target_role"]
        )
        features.append(
            FusionFeature(
                Event_id=row["event_id"],
                Gene=gene,
                Role=row["target_role"],
                Junction_position_aa=row["breakpoint_protein_position"],
                Domain_retention_flags={"kinase": row["domain_status"]},
                Domain_retention_details={"kinase": detail},
            )
        )

    (
        expected_mapped,
        expected_total,
        expected_retained,
        expected_percent,
        expected_table,
        expected_p,
    ) = expected
    result = DomainRetentionAlgorithm().run(
        events,
        features,
        load_gene_config(gene.lower()),
        {"seed": 42, "n_permutations": 10},
    )

    assert len(features) == expected_mapped
    retained = sum(f.Domain_retention_flags["kinase"] == "retained" for f in features)
    assert retained == expected_retained
    assert retained / expected_total * 100 == pytest.approx(expected_percent)
    assert result.Tables["frame_domain_contingency_table"] == expected_table
    assert result.Summary["fisher_p_value"] == pytest.approx(expected_p)


@pytest.mark.parametrize("gene", ["braf", "ret"])
def test_verified_live_benchmark_conclusions_are_unchanged(gene):
    """Replay real mapped inputs through domain mapping and the statistical pipeline.

    The committed runs preserve mapped event rows, not the raw API responses or
    canonical genomic transcripts. This regression therefore starts at the real
    protein breakpoints; it does not replay ingestion or genomic mapping.
    """
    candidates = sorted(RUNS_DIR.glob(f"{gene}_msk-impact-50k-2026_*"))
    assert candidates, f"Missing committed {gene.upper()} benchmark run"
    run_dir = candidates[-1]
    payload = json.loads((run_dir / "results.json").read_text())
    with (run_dir / "results.tsv").open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert payload["gene_symbol"] == gene.upper()
    assert payload["study_id"] == "msk_impact_50k_2026"
    assert rows
    assert len(rows) == len(payload["events"])
    config = load_gene_config(gene)
    source = MagicMock()
    source.fetch.return_value = [
        ProteinDomain(**domain, source="genome_nexus")
        for domain in payload["gene_track"]["domains"]
    ]
    events = []
    features = []
    for row, json_row in zip(rows, payload["events"], strict=True):
        event = FusionEvent(
            Event_id=row["event_id"],
            Cohort=payload["study_id"],
            Sample_id=row["sample_id"],
            Fusion_name=row["fusion_name"],
            Frame_status=row["frame_status"],
            Is_protein_fusion=True,
        )
        feature = map_event(
            event,
            config,
            role=row["target_role"],
            junction_position_aa=int(row["breakpoint_protein_position"]),
            domain_source=source,
        )
        # Expected domain calls are outputs only: never feed them into map_event.
        assert event.Event_id == json_row["event_id"]
        assert event.Frame_status == json_row["frame_status"]
        assert feature.Role == json_row["target_role"]
        assert feature.Junction_position_aa == json_row["breakpoint_protein_position"]
        assert (
            feature.Domain_retention_flags["kinase"]
            == (row["domain_status"])
            == json_row["domain_status"]
        )
        detail = feature.Domain_retention_details["kinase"]
        assert (
            detail.Retained_fraction
            == (float(row["domain_retained_fraction"]))
            == json_row["domain_retained_fraction"]
        )
        assert (
            detail.Is_truncated
            == (row["domain_is_truncated"] == "True")
            == json_row["domain_is_truncated"]
        )
        events.append(event)
        features.append(feature)

    expected = next(
        result
        for result in payload["algorithm_results"]
        if result["Algorithm"] == "domain_retention"
    )
    result = DomainRetentionAlgorithm().run(events, features, config, expected["Parameters"])
    summary = payload["summary"]
    assert len(features) == summary["mapped_fusions"]
    assert sum(event.Frame_status == "in-frame" for event in events) == summary["in_frame_count"]
    retained = sum(feature.Domain_retention_flags["kinase"] == "retained" for feature in features)
    assert retained == summary["kinase_retained_count"]
    # Unmapped records (one in BRAF) are absent from both event artifacts, so
    # percentages are denominated by the mapped population, not total_fusions
    # (see #30).
    assert 100.0 * retained / summary["mapped_fusions"] == summary["kinase_retained_percent"]
    assert (
        result.Tables["frame_domain_contingency_table"]
        == (summary["frame_domain_contingency_table"])
        == expected["Tables"]["frame_domain_contingency_table"]
    )
    assert (
        result.Tables["domain_retention_descriptives"]
        == (expected["Tables"]["domain_retention_descriptives"])
    )
    for key in ("fisher_p_value", "fisher_odds_ratio"):
        assert result.Summary[key] == summary[key] == expected["Summary"][key]
    assert (
        result.Summary["observed_in_frame_retention_rate"]
        == (expected["Summary"]["observed_in_frame_retention_rate"])
    )
    # Do not compare permutation p-values: live runs sample genomic breakpoints
    # using transcripts not persisted here; offline replay resamples observed
    # protein positions. These are different null models, not rounding error
    # that a pytest.approx tolerance could legitimately accommodate.


def test_corrected_ret_live_artifact_summary():
    """Lock RET conclusions directly to the post-locus-validation live artifact."""
    payload = json.loads((latest_run_dir("ret_msk-impact-50k-2026") / "results.json").read_text())
    summary = payload["summary"]

    assert summary["total_fusions"] == 194
    assert summary["mapped_fusions"] == 194
    assert summary["in_frame_count"] == 146
    assert summary["kinase_retained_count"] == 179
    assert summary["kinase_retained_percent"] == pytest.approx(92.26804123711341)
    assert summary["frame_domain_contingency_table"] == [[141, 38], [5, 10]]
    assert summary["fisher_p_value"] == pytest.approx(0.00041966557652448966)


def test_retention_descriptives_separate_truncated_from_fully_lost():
    config = load_gene_config("braf")
    events = []
    features = []
    for index, (breakpoint, status) in enumerate(
        [(712, "retained"), (600, "disrupted"), (457, "lost")]
    ):
        event_id = f"event-{index}"
        events.append(
            FusionEvent(
                Event_id=event_id,
                Cohort="test",
                Frame_status="in-frame",
                Is_protein_fusion=True,
            )
        )
        features.append(
            FusionFeature(
                Event_id=event_id,
                Gene="BRAF",
                Junction_position_aa=breakpoint,
                Domain_retention_flags={"kinase": status},
                Domain_retention_details={
                    "kinase": calculate_domain_retention(458, 712, breakpoint, "five_prime")
                },
            )
        )

    result = DomainRetentionAlgorithm().run(
        events, features, config, {"seed": 42, "n_permutations": 10}
    )
    row = result.Tables["domain_retention_descriptives"][0]

    assert row["Fully_retained_count"] == 1
    assert row["Truncated_count"] == 1
    assert row["Fully_lost_count"] == 1
    assert row["Mean_retained_fraction_among_non_retained_calls"] == pytest.approx((143 / 255) / 2)
