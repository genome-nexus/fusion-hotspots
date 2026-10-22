"""Tests for scoring a cohort scan against textbook positive controls."""

from __future__ import annotations

import csv
import json

import pytest
from click.testing import CliRunner

from cfh.benchmarks.positive_controls import (
    FAIL,
    NOT_EVALUABLE,
    PASS,
    evaluate_positive_controls,
    load_positive_controls,
)
from cfh.cli import main

CONTROLS = {
    "breakpoint_controls": [
        {
            "gene": "ALK",
            "role": "three_prime",
            "junction_exons": [20],
            "partners": ["EML4"],
            "source": "Soda 2007",
        },
        {"gene": "NTRK3", "role": "three_prime", "partners": ["ETV6"], "source": "Knezevich"},
    ],
    "cooccurrence_controls": [
        {
            "fusion_gene": "ERG",
            "comparator_gene": "SPOP",
            "alteration_type": "mutation",
            "direction": "mutually_exclusive",
            "source": "Barbieri 2012",
        },
        {
            "fusion_gene": "ALK",
            "comparator_gene": "KRAS",
            "alteration_type": "mutation",
            "direction": "mutually_exclusive",
            "source": "Gainor 2013",
        },
    ],
}


def _write_tsv(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _event(exon, partner, role="three_prime"):
    return {"target_role": role, "breakpoint_exon": exon, "partner_gene": partner}


@pytest.fixture
def scan_dir(tmp_path):
    scan = tmp_path / "cohort-scan_x" / "cohort_scan"
    (scan / "gene_reports" / "alk").mkdir(parents=True)
    _write_tsv(
        scan / "fdr_family.tsv",
        [
            {
                "gene": "ALK",
                "algorithm": "domain_retention",
                "test": "permutation",
                "bh_adjusted_q": "0.001",
            },
            {
                "gene": "NTRK3",
                "algorithm": "domain_retention",
                "test": "permutation",
                "bh_adjusted_q": "0.3",
            },
        ],
    )
    _write_tsv(
        scan / "junction_concentration_family.tsv",
        [
            {
                "gene": "ALK",
                "algorithm": "junction_concentration",
                "test": "permutation",
                "bh_adjusted_q": "0.002",
            }
        ],
    )
    events = [_event(20, "EML4")] * 8 + [_event(19, "EML4"), _event(20, "KIF5B")]
    (scan / "gene_reports" / "alk" / "results.json").write_text(json.dumps({"events": events}))
    _write_tsv(
        scan / "cooccurrence_discovery.tsv",
        [
            {
                "fusion_gene": "ERG",
                "comparator_gene": "SPOP",
                "alteration_type": "mutation",
                "direction": "mutually_exclusive",
                "cmh_q_value": "1e-6",
                "both_samples": "2",
                "expected_both_samples": "40.5",
            },
        ],
    )
    return scan.parent


def _by_key(checks):
    return {(check["control"], check["check"]): check for check in checks}


def test_scores_each_check_against_the_scan(scan_dir):
    checks = _by_key(evaluate_positive_controls(scan_dir, controls=CONTROLS))

    assert checks[("ALK", "selection")]["status"] == PASS
    assert checks[("ALK", "junction")]["status"] == PASS
    assert "90% of 10" in checks[("ALK", "junction")]["detail"]
    assert checks[("ALK", "partner")]["status"] == PASS
    assert checks[("ALK", "concentration")]["status"] == PASS
    assert (
        "junction-concentration permutation q = 0.002"
        in (checks[("ALK", "concentration")]["detail"])
    )
    assert checks[("NTRK3", "selection")]["status"] == FAIL
    assert checks[("NTRK3", "concentration")]["status"] == NOT_EVALUABLE
    # No per-gene report: the partner check cannot be scored, and NTRK3 has no
    # junction expectation at all.
    assert checks[("NTRK3", "partner")]["status"] == NOT_EVALUABLE
    assert ("NTRK3", "junction") not in checks
    assert checks[("ERG fusion + SPOP mutation", "cooccurrence")]["status"] == PASS
    assert checks[("ALK fusion + KRAS mutation", "cooccurrence")]["status"] == FAIL


def test_junction_fails_when_the_modal_exon_is_wrong(scan_dir):
    controls = {
        "breakpoint_controls": [{**CONTROLS["breakpoint_controls"][0], "junction_exons": [19]}]
    }

    checks = _by_key(evaluate_positive_controls(scan_dir, controls=controls))

    assert checks[("ALK", "junction")]["status"] == FAIL


def test_cooccurrence_fails_in_the_wrong_direction(scan_dir):
    control = {**CONTROLS["cooccurrence_controls"][0], "direction": "co_occurring"}

    checks = evaluate_positive_controls(scan_dir, controls={"cooccurrence_controls": [control]})

    assert checks[0]["status"] == FAIL


def test_requires_the_fdr_family_table(tmp_path):
    with pytest.raises(FileNotFoundError, match="fdr_family.tsv"):
        evaluate_positive_controls(tmp_path, controls=CONTROLS)


def test_packaged_controls_are_well_formed():
    controls = load_positive_controls()

    for control in controls["breakpoint_controls"]:
        assert control["role"] in {"three_prime", "five_prime"}
        assert control["partners"] and control["source"]
    for control in controls["cooccurrence_controls"]:
        assert control["direction"] in {"co_occurring", "mutually_exclusive"}
        assert control["alteration_type"] in {"mutation", "amplification", "deep_deletion"}


def test_cli_writes_the_scored_table(scan_dir):
    result = CliRunner().invoke(main, ["eval-controls", str(scan_dir)])

    assert result.exit_code == 0, result.output
    written = scan_dir / "cohort_scan" / "positive_controls.tsv"
    assert written.exists()
    assert "passed" in result.output
