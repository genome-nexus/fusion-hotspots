import csv
import json

import pytest
from click.testing import CliRunner

from cfh.cli import main
from cfh.gene_comparison import collect_p_values, compare_gene_runs


def _write_artifact(path, gene, study, algorithm_results):
    path.mkdir()
    (path / "results.json").write_text(
        json.dumps(
            {
                "gene_symbol": gene,
                "study_id": study,
                "algorithm_results": algorithm_results,
            }
        )
    )


def test_compare_gene_runs_collects_applicable_final_p_values_and_adjusts_together(tmp_path):
    first = tmp_path / "first"
    second = tmp_path / "second"
    _write_artifact(
        first,
        "GENE_A",
        "study_one",
        [
            {
                "Algorithm": "domain_retention",
                "Summary": {
                    "fisher_p_value": 0.01,
                    "permutation_empirical_p_value": 0.04,
                },
            },
            {"Algorithm": "frequency", "Summary": {"event_count": 2}},
        ],
    )
    _write_artifact(
        second,
        "GENE_B",
        "study_two",
        [
            {
                "Algorithm": "cutpoint_detection",
                "Summary": {"observed_p_value": 0.0001, "corrected_p_value": 0.03},
            }
        ],
    )

    rows = compare_gene_runs([first, second])

    assert [(row["algorithm"], row["test"]) for row in rows] == [
        ("domain_retention", "fisher"),
        ("domain_retention", "permutation"),
        ("cutpoint_detection", "permutation_corrected"),
    ]
    assert [row["bh_adjusted_q"] for row in rows] == pytest.approx([0.03, 0.04, 0.04])
    assert [row["significant_at_q_lt_0_05"] for row in rows] == ["y", "y", "y"]


def test_collect_p_values_rejects_malformed_artifact(tmp_path):
    path = tmp_path / "results.json"
    path.write_text("{}")

    with pytest.raises(ValueError, match="gene_symbol and study_id"):
        collect_p_values([path])


def test_compare_genes_cli_writes_tsv(tmp_path):
    run = tmp_path / "run"
    _write_artifact(
        run,
        "GENE_A",
        "study_one",
        [
            {
                "Algorithm": "joint_partner",
                "Summary": {"p_value": 0.02},
            }
        ],
    )
    output = tmp_path / "comparison.tsv"

    result = CliRunner().invoke(main, ["compare-genes", str(run), "--output", str(output)])

    assert result.exit_code == 0, result.output
    assert "Adjusted 1 p-values" in result.output
    rows = list(csv.DictReader(output.open(), delimiter="\t"))
    assert rows[0]["raw_p"] == "0.02"
    assert rows[0]["bh_adjusted_q"] == "0.02"
    assert rows[0]["significant_at_q_lt_0_05"] == "y"


def test_collects_window_expression_and_each_mutation_comparator_once(tmp_path):
    run = tmp_path / "run"
    _write_artifact(
        run,
        "GENE_A",
        "study_one",
        [
            {"Algorithm": "window_detection", "Summary": {"corrected_p_value": 0.01}},
            {
                "Algorithm": "expression_association",
                "Summary": {
                    "fusion_positive_vs_negative": {"p_value": 0.02},
                    "domain_retention_split": {"p_value": 0.03},
                },
            },
            {
                "Algorithm": "mutation_cooccurrence",
                "Summary": {
                    "targets": [
                        {
                            "comparator_gene": "BRAF",
                            "alteration_type": "point_mutation",
                            "protein_change": "V600E",
                            "p_value": 0.04,
                        },
                        {
                            "comparator_gene": "BRAF",
                            "alteration_type": "cna_amp",
                            "protein_change": None,
                            "p_value": 0.05,
                        },
                    ]
                },
                "Tables": {
                    "cooccurrence_results": [
                        {"comparator_gene": "BRAF", "p_value": 0.04},
                        {"comparator_gene": "BRAF", "p_value": 0.05},
                    ]
                },
            },
        ],
    )

    rows = collect_p_values([run])

    assert [(row["algorithm"], row["test"], row["raw_p"]) for row in rows] == [
        ("window_detection", "permutation_corrected", 0.01),
        ("expression_association", "fusion_positive_vs_negative", 0.02),
        ("expression_association", "domain_retention_split", 0.03),
        ("mutation_cooccurrence", "cooccurrence:BRAF:point_mutation:V600E", 0.04),
        (
            "mutation_cooccurrence",
            "cooccurrence:BRAF:cna_amp:any_protein_change",
            0.05,
        ),
    ]


@pytest.mark.parametrize("p_value", [float("nan"), float("inf"), -float("inf"), 1.1, -0.1])
def test_collect_p_values_rejects_nonfinite_or_out_of_range_nested_values(tmp_path, p_value):
    run = tmp_path / "run"
    _write_artifact(
        run,
        "GENE_A",
        "study_one",
        [
            {
                "Algorithm": "expression_association",
                "Summary": {"fusion_positive_vs_negative": {"p_value": p_value}},
            }
        ],
    )

    with pytest.raises(ValueError, match="finite and between 0 and 1"):
        collect_p_values([run])


def test_missing_optional_comparators_are_skipped(tmp_path):
    run = tmp_path / "run"
    _write_artifact(
        run,
        "GENE_A",
        "study_one",
        [
            {
                "Algorithm": "expression_association",
                "Summary": {"fusion_positive_vs_negative": {"p_value": None}},
            },
            {"Algorithm": "mutation_cooccurrence", "Summary": {"targets": []}},
        ],
    )

    assert collect_p_values([run]) == []


def test_new_hypotheses_join_the_same_bh_family(tmp_path):
    first = tmp_path / "first"
    second = tmp_path / "second"
    _write_artifact(
        first,
        "GENE_A",
        "study_one",
        [
            {"Algorithm": "window_detection", "Summary": {"corrected_p_value": 0.01}},
            {
                "Algorithm": "mutation_cooccurrence",
                "Summary": {
                    "targets": [
                        {
                            "comparator_gene": "BRAF",
                            "alteration_type": "point_mutation",
                            "protein_change": "V600E",
                            "p_value": 0.04,
                        }
                    ]
                },
            },
        ],
    )
    _write_artifact(
        second,
        "GENE_B",
        "study_two",
        [
            {
                "Algorithm": "expression_association",
                "Summary": {"fusion_positive_vs_negative": {"p_value": 0.02}},
            }
        ],
    )

    rows = compare_gene_runs([first, second])

    assert [row["bh_adjusted_q"] for row in rows] == pytest.approx([0.03, 0.04, 0.03])
