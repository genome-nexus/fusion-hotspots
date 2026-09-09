"""Offline artifact/CLI contract and regression on exactly three committed runs."""

import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from cfh.algorithms.cross_cohort_concordance import compare_cohort_runs
from cfh.cli import main

ROOT = Path(__file__).resolve().parents[1]
BRAF_RUNS = [
    ROOT / "runs/braf_msk-impact-50k-2026_20260909T181926Z",
    ROOT / "runs/braf_msk-impact-2017_20260905T012645Z",
    ROOT / "runs/braf_thca-tcga-pan-can-atlas-2018_20260909T034447Z",
]


def test_real_braf_artifacts():
    report = compare_cohort_runs(BRAF_RUNS)
    assert [s["frame_domain_contingency_table"] for s in report["strata"]] == [
        [[142, 21], [9, 6]],
        [[31, 2], [2, 5]],
        [[9, 0], [6, 0]],
    ]
    assert report["cmh"]["statistic"] == pytest.approx(21.36486008160812)
    assert report["cmh"]["p_value"] == pytest.approx(3.796664804980965e-06, rel=1e-10)
    assert report["cmh"]["informative_strata"] == 2
    assert report["cmh"]["total_strata"] == 3
    assert report["cmh"]["common_odds_ratio"] == pytest.approx(7.45527079303675)
    assert report["nominal_pooled_significant"] is True
    assert report["all_cohorts_positive_and_informative"] is False
    assert report["sample_overlaps"] == [
        {
            "study_ids": ["msk_impact_50k_2026", "msk_impact_2017"],
            "shared_sample_id_count": 34,
        }
    ]
    assert any("violate" in w for w in report["warnings"])


def test_cli_offline_json(tmp_path, monkeypatch):
    def forbidden(*args, **kwargs):
        pytest.fail("Offline comparison must not ingest data")

    monkeypatch.setattr("cfh.cli.run_analysis", forbidden)
    monkeypatch.setattr("cfh.cli.run_real_benchmark", forbidden)
    monkeypatch.setattr("requests.sessions.Session.request", forbidden)
    output = tmp_path / "report.json"
    result = CliRunner().invoke(
        main, ["compare-cohorts", *map(str, BRAF_RUNS), "--output", str(output)]
    )
    assert result.exit_code == 0, result.output
    assert json.loads(output.read_text()) == json.loads(result.output)


def write_run(tmp_path, study, *, gene="OTHER", table=None):
    path = tmp_path / f"{study}.json"
    path.write_text(
        json.dumps(
            {
                "gene_symbol": gene,
                "study_id": study,
                "summary": {"frame_domain_contingency_table": table or [[10, 2], [2, 10]]},
            }
        )
    )
    return path


def test_gene_agnostic_and_infinity_json(tmp_path):
    paths = [write_run(tmp_path, s, table=[[10, 0], [0, 10]]) for s in ("one", "two")]
    result = CliRunner().invoke(main, ["compare-cohorts", *map(str, paths)])
    assert result.exit_code == 0
    report = json.loads(result.output)
    assert report["gene_symbol"] == "OTHER"
    assert report["cmh"]["common_odds_ratio"] == "infinity"
    assert report["all_cohorts_positive_and_informative"] is True


def test_reject_duplicate_cohort_and_mixed_genes(tmp_path):
    first = write_run(tmp_path, "one")
    second = write_run(tmp_path, "two", gene="DIFFERENT")
    with pytest.raises(ValueError, match="Duplicate cohort"):
        compare_cohort_runs([first, first])
    with pytest.raises(ValueError, match="same gene"):
        compare_cohort_runs([first, second])
    with pytest.raises(ValueError, match="at least two"):
        compare_cohort_runs([first])


@pytest.mark.parametrize("payload", ["{}", "null", "not json", '{"gene_symbol": "X"}'])
def test_bad_artifacts_are_cli_errors(tmp_path, payload):
    bad = tmp_path / "bad.json"
    bad.write_text(payload)
    good = write_run(tmp_path, "good")
    result = CliRunner().invoke(main, ["compare-cohorts", str(bad), str(good)])
    assert result.exit_code == 1
    assert "Invalid run artifact" in result.output
