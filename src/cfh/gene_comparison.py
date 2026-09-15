"""Offline cross-gene multiple-testing reports from existing run artifacts."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any

from cfh.stats.multiple_testing import benjamini_hochberg

SIGNIFICANCE_LEVEL = 0.05

# These are the final inferential p-values exposed by registered algorithms.
# A label is included because some algorithms expose independent test families.
_SUMMARY_P_VALUES: dict[str, tuple[tuple[str, tuple[str, ...]], ...]] = {
    "domain_retention": (
        ("fisher", ("fisher_p_value",)),
        ("permutation", ("permutation_empirical_p_value",)),
    ),
    "domain_disruption": (
        ("fisher", ("fisher_p_value",)),
        ("permutation", ("permutation_empirical_p_value",)),
    ),
    # corrected_p_value already accounts for scanning candidate cutpoints.
    "cutpoint_detection": (("permutation_corrected", ("corrected_p_value",)),),
    # corrected_p_value already accounts for scanning candidate windows.
    "window_detection": (("permutation_corrected", ("corrected_p_value",)),),
    "joint_partner": (("enrichment", ("p_value",)),),
    "confidence_stats": (("welch_t_test", ("ttest", "p_value")),),
    "expression_association": (
        ("fusion_positive_vs_negative", ("fusion_positive_vs_negative", "p_value")),
        ("domain_retention_split", ("domain_retention_split", "p_value")),
    ),
}


def _nested_value(mapping: dict[str, Any], path: tuple[str, ...]) -> Any:
    value: Any = mapping
    for key in path:
        if not isinstance(value, dict):
            return None
        value = value.get(key)
    return value


def _results_path(run_artifact: Path) -> Path:
    return run_artifact / "results.json" if run_artifact.is_dir() else run_artifact


def _validate_p_value(
    raw_p: Any, *, algorithm: str, value_path: tuple[str, ...], source: str
) -> float | None:
    if raw_p is None:
        return None
    if isinstance(raw_p, bool) or not isinstance(raw_p, (int, float)):
        raise ValueError(
            f"P-value {'.'.join(value_path)} for {algorithm} in {source} must be numeric or null"
        )
    value = float(raw_p)
    if not math.isfinite(value) or not 0.0 <= value <= 1.0:
        raise ValueError(
            f"P-value {'.'.join(value_path)} for {algorithm} in {source} "
            f"must be finite and between 0 and 1; got {value!r}"
        )
    return value


def _mutation_target_test_label(target: dict[str, Any]) -> str:
    """Return a stable, unique label for one configured comparator target."""
    return "cooccurrence:{}:{}:{}".format(
        target.get("comparator_gene", ""),
        target.get("alteration_type", ""),
        target.get("protein_change") or "any_protein_change",
    )


def collect_p_values_from_algorithm_results(
    gene: str,
    study: str,
    algorithm_results: list[Any],
    *,
    source: str,
) -> list[dict[str, Any]]:
    """Collect final scalar p-values straight from already-computed results.

    Shared by :func:`collect_p_values` (which reads them back from an
    on-disk ``results.json``) and any in-memory caller -- e.g. a
    genome-wide scan that runs many genes through the same orchestrator
    without necessarily writing each one's own ``results.json`` first.
    ``algorithm_results`` entries may be plain dicts (as read from JSON) or
    objects exposing ``model_dump()`` (an :class:`~cfh.model.algorithm_result.AlgorithmResult`).
    """
    rows: list[dict[str, Any]] = []
    for result in algorithm_results:
        if not isinstance(result, dict):
            if hasattr(result, "model_dump"):
                result = result.model_dump(mode="json")
            else:
                continue
        algorithm = result.get("Algorithm")
        summary = result.get("Summary") or {}
        for test, value_path in _SUMMARY_P_VALUES.get(algorithm, ()):
            raw_p = _nested_value(summary, value_path)
            raw_p = _validate_p_value(
                raw_p, algorithm=algorithm, value_path=value_path, source=source
            )
            if raw_p is None:
                continue
            rows.append(
                {
                    "gene": gene,
                    "study": study,
                    "algorithm": algorithm,
                    "test": test,
                    "raw_p": raw_p,
                    "source": source,
                }
            )

        # ``Tables["cooccurrence_results"]`` repeats these rows verbatim. Read
        # the canonical Summary representation once so each comparator is one
        # hypothesis in the family.
        if algorithm == "mutation_cooccurrence":
            targets = summary.get("targets")
            if not isinstance(targets, list):
                continue
            for target in targets:
                if not isinstance(target, dict):
                    continue
                value_path = ("targets", "p_value")
                raw_p = _validate_p_value(
                    target.get("p_value"),
                    algorithm=algorithm,
                    value_path=value_path,
                    source=source,
                )
                if raw_p is None:
                    continue
                rows.append(
                    {
                        "gene": gene,
                        "study": study,
                        "algorithm": algorithm,
                        "test": _mutation_target_test_label(target),
                        "raw_p": raw_p,
                        "source": source,
                    }
                )
    return rows


def collect_p_values(run_artifacts: list[Path]) -> list[dict[str, Any]]:
    """Collect final scalar p-values from existing ``results.json`` artifacts."""
    rows: list[dict[str, Any]] = []
    for run_artifact in run_artifacts:
        path = _results_path(run_artifact)
        try:
            payload = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError) as exc:
            raise ValueError(f"Could not read run artifact {path}: {exc}") from exc

        gene = payload.get("gene_symbol")
        study = payload.get("study_id")
        algorithm_results = payload.get("algorithm_results")
        if not isinstance(gene, str) or not isinstance(study, str):
            raise ValueError(f"Run artifact {path} must contain gene_symbol and study_id")
        if not isinstance(algorithm_results, list):
            raise ValueError(f"Run artifact {path} must contain an algorithm_results list")

        rows.extend(
            collect_p_values_from_algorithm_results(
                gene, study, algorithm_results, source=str(path)
            )
        )
    return rows


def compare_gene_runs(run_artifacts: list[Path]) -> list[dict[str, Any]]:
    """Collect and BH-adjust all applicable p-values as one hypothesis family."""
    rows = collect_p_values(run_artifacts)
    hypotheses = [(row["gene"], f"{row['algorithm']}:{row['test']}", row["raw_p"]) for row in rows]
    adjusted = benjamini_hochberg(hypotheses)
    for row, (_, _, _, q_value) in zip(rows, adjusted, strict=True):
        row["bh_adjusted_q"] = q_value
        row["significant_at_q_lt_0_05"] = "y" if q_value < SIGNIFICANCE_LEVEL else "n"
    return rows


def write_comparison_tsv(rows: list[dict[str, Any]], output_path: Path) -> None:
    """Write a stable, machine-readable adjusted-p-value table."""
    fieldnames = [
        "gene",
        "study",
        "algorithm",
        "test",
        "raw_p",
        "bh_adjusted_q",
        "significant_at_q_lt_0_05",
        "source",
    ]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
