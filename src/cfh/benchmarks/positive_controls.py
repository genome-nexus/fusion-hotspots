"""Score a finished cohort scan against textbook fusion biology.

A scan that cannot recover canonical findings (EML4-ALK at ALK exon 20,
ETS fusions excluding SPOP mutations, ...) gives no reason to trust its
calls on genes nobody has curated. :data:`CONTROLS_PATH` lists those
findings with their literature sources; :func:`evaluate_positive_controls`
checks each against a scan's outputs and reports ``pass``, ``fail`` or
``not_evaluable`` with the evidence behind the verdict.
"""

from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path
from typing import Any

import yaml

CONTROLS_PATH = Path(__file__).parent / "positive_controls.yaml"
DEFAULT_MIN_JUNCTION_FRACTION = 0.5

PASS = "pass"
FAIL = "fail"
NOT_EVALUABLE = "not_evaluable"


def load_positive_controls(path: Path = CONTROLS_PATH) -> dict[str, list[dict]]:
    return yaml.safe_load(path.read_text())


def _scan_directory(run_directory: Path) -> Path:
    nested = run_directory / "cohort_scan"
    return nested if nested.is_dir() else run_directory


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _float(value: str | None) -> float | None:
    try:
        return float(value) if value not in (None, "") else None
    except ValueError:
        return None


def _check(control: str, check: str, status: str, detail: str, source: str) -> dict[str, str]:
    return {
        "control": control,
        "check": check,
        "status": status,
        "detail": detail,
        "source": source,
    }


def _fdr_q(fdr_rows: list[dict[str, str]], gene: str, algorithm: str, test: str) -> float | None:
    for row in fdr_rows:
        if row["gene"] == gene and row["algorithm"] == algorithm and row["test"] == test:
            return _float(row["bh_adjusted_q"])
    return None


# check name -> (algorithm, label) of the permutation test whose q-value it reads.
_FDR_CHECKS = {
    "selection": ("domain_retention", "retention"),
    "concentration": ("junction_concentration", "junction-concentration"),
}


def _fdr_check(control, fdr_rows, alpha, check) -> dict[str, str]:
    gene, source = control["gene"], control["source"]
    algorithm, label = _FDR_CHECKS[check]
    q_value = _fdr_q(fdr_rows, gene, algorithm, "permutation")
    if q_value is None:
        return _check(
            gene,
            check,
            NOT_EVALUABLE,
            f"no {label} permutation q-value in the scan",
            source,
        )
    status = PASS if q_value < alpha else FAIL
    return _check(gene, check, status, f"{label} permutation q = {q_value:.2g}", source)


def _junction_check(control, events, min_fraction) -> dict[str, str]:
    gene, source = control["gene"], control["source"]
    expected = set(control["junction_exons"])
    exons = Counter(
        event["breakpoint_exon"]
        for event in events
        if event.get("target_role") == control["role"] and event.get("breakpoint_exon") is not None
    )
    total = sum(exons.values())
    if not total:
        return _check(
            gene, "junction", NOT_EVALUABLE, f"no {control['role']} events with an exon", source
        )
    modal_exon, _ = exons.most_common(1)[0]
    fraction = sum(count for exon, count in exons.items() if exon in expected) / total
    status = PASS if modal_exon in expected and fraction >= min_fraction else FAIL
    expected_text = "/".join(str(exon) for exon in sorted(expected))
    detail = (
        f"modal exon {modal_exon}; {fraction:.0%} of {total} {control['role']} events "
        f"at exon {expected_text}"
    )
    return _check(gene, "junction", status, detail, source)


def _partner_check(control, events) -> dict[str, str]:
    gene, source = control["gene"], control["source"]
    partners = Counter(event["partner_gene"] for event in events if event.get("partner_gene"))
    if not partners:
        return _check(gene, "partner", NOT_EVALUABLE, "no events with a named partner", source)
    top, count = partners.most_common(1)[0]
    status = PASS if top in control["partners"] else FAIL
    detail = f"most frequent partner {top} ({count} of {sum(partners.values())} events)"
    return _check(gene, "partner", status, detail, source)


def _breakpoint_checks(control, scan_directory, fdr_rows, alpha, min_fraction):
    checks = [_fdr_check(control, fdr_rows, alpha, check) for check in _FDR_CHECKS]
    results_path = scan_directory / "gene_reports" / control["gene"].lower() / "results.json"
    needed = (["junction"] if "junction_exons" in control else []) + ["partner"]
    if not results_path.exists():
        for check in needed:
            checks.append(
                _check(
                    control["gene"],
                    check,
                    NOT_EVALUABLE,
                    "the scan wrote no per-gene report for this gene",
                    control["source"],
                )
            )
        return checks
    events = json.loads(results_path.read_text())["events"]
    if "junction_exons" in control:
        checks.append(_junction_check(control, events, min_fraction))
    checks.append(_partner_check(control, events))
    return checks


def _cooccurrence_check(control, discovery_rows, alpha) -> dict[str, str]:
    fusion, comparator = control["fusion_gene"], control["comparator_gene"]
    name = f"{fusion} fusion + {comparator} {control['alteration_type']}"
    source = control["source"]
    if discovery_rows is None:
        return _check(
            name, "cooccurrence", NOT_EVALUABLE, "scan ran without co-occurrence discovery", source
        )
    row = next(
        (
            row
            for row in discovery_rows
            if row["fusion_gene"] == control["fusion_gene"]
            and row["comparator_gene"] == control["comparator_gene"]
            and row["alteration_type"] == control["alteration_type"]
        ),
        None,
    )
    if row is None:
        return _check(
            name,
            "cooccurrence",
            FAIL,
            "not among the written discovery rows (q >= 0.25 and not a top pair, or untested)",
            source,
        )
    q_value = _float(row["cmh_q_value"])
    direction = row["direction"]
    passed = q_value is not None and q_value < alpha and direction == control["direction"]
    detail = (
        f"{direction}, {row['both_samples']} vs {_float(row['expected_both_samples']) or 0:.1f} "
        f"expected, q = {q_value:.2g}"
        if q_value is not None
        else f"{direction}, no q-value"
    )
    return _check(name, "cooccurrence", PASS if passed else FAIL, detail, source)


def evaluate_positive_controls(
    run_directory: Path,
    *,
    controls: dict[str, list[dict]] | None = None,
    alpha: float = 0.05,
    min_junction_fraction: float = DEFAULT_MIN_JUNCTION_FRACTION,
) -> list[dict[str, str]]:
    """Return one row per check: control, check, status, detail, source."""
    controls = controls or load_positive_controls()
    scan_directory = _scan_directory(Path(run_directory))
    fdr_path = scan_directory / "fdr_family.tsv"
    if not fdr_path.exists():
        raise FileNotFoundError(
            f"{fdr_path} not found; positive controls need a scan that writes fdr_family.tsv"
        )
    fdr_rows = _read_tsv(fdr_path)
    # Junction concentration is its own family; its q-values live beside the main one.
    concentration_path = scan_directory / "junction_concentration_family.tsv"
    if concentration_path.exists():
        fdr_rows += _read_tsv(concentration_path)
    discovery_path = scan_directory / "cooccurrence_discovery.tsv"
    discovery_rows = _read_tsv(discovery_path) if discovery_path.exists() else None

    checks: list[dict[str, str]] = []
    for control in controls.get("breakpoint_controls", []):
        checks.extend(
            _breakpoint_checks(control, scan_directory, fdr_rows, alpha, min_junction_fraction)
        )
    for control in controls.get("cooccurrence_controls", []):
        checks.append(_cooccurrence_check(control, discovery_rows, alpha))
    return checks


def summarize(checks: list[dict[str, Any]]) -> Counter:
    return Counter(check["status"] for check in checks)
