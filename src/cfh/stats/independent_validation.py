"""Offline, descriptive held-out ranking comparison with explicit labels.

Input scores must have been frozen using discovery patients only. This tool
checks declared patient identity metadata and computes ranking metrics; it
cannot audit how scores or external truth labels were produced.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any


def _patient_ids(metadata: dict, name: str) -> set[str]:
    ids = metadata.get("patient_ids")
    if not isinstance(ids, list) or not ids or any(not isinstance(x, str) or not x for x in ids):
        raise ValueError(f"{name} requires a nonempty complete patient_ids list")
    if len(ids) != len(set(ids)):
        raise ValueError(f"{name} has duplicate patient IDs")
    return set(ids)


def _scores(candidates: list[dict], key: str) -> list[tuple[float, int, str]]:
    rows = []
    for candidate in candidates:
        if not isinstance(candidate, dict):
            raise ValueError("Every candidate must be an object")
        value = candidate.get(key)
        if (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(value)
        ):
            raise ValueError(f"Every candidate requires a finite {key}")
        rows.append((float(value), candidate["label"], candidate["candidate_id"]))
    return sorted(rows, key=lambda row: (-row[0], row[2]))


def _ranking_metrics(rows: list[tuple[float, int, str]], k: int) -> dict[str, float | int]:
    positives = sum(label for _, label, _ in rows)
    seen = 0
    true_positives = 0
    average_precision = 0.0
    expected_top_k_positives = 0.0
    tie_groups = 0
    cursor = 0
    while cursor < len(rows):
        end = cursor + 1
        while end < len(rows) and rows[end][0] == rows[cursor][0]:
            end += 1
        group_size = end - cursor
        group_positives = sum(row[1] for row in rows[cursor:end])
        if group_size > 1:
            tie_groups += 1
        seen += group_size
        true_positives += group_positives
        average_precision += (group_positives / positives) * (true_positives / seen)
        in_top_k = max(0, min(k, end) - cursor)
        expected_top_k_positives += in_top_k * group_positives / group_size
        cursor = end
    return {
        "average_precision": average_precision,
        "precision_at_k": expected_top_k_positives / min(k, len(rows)),
        "tie_groups": tie_groups,
    }


def compare_heldout_rankings(
    discovery: dict[str, Any], validation: dict[str, Any], *, k: int = 5
) -> dict[str, Any]:
    """Compare frozen composite and recurrence scores on one labeled candidate set.

    Both inputs declare ``patient_id_namespace`` and complete ``patient_ids``.
    Validation additionally has ``candidates`` with unique ``candidate_id``,
    binary ``label``, ``composite_score`` and ``recurrence_score``. A shared
    patient namespace and no ID overlap are required. Labels are external
    ground truth; the function never manufactures them.
    """
    if k <= 0:
        raise ValueError("k must be positive")
    if not isinstance(discovery, dict) or not isinstance(validation, dict):
        raise ValueError("Discovery and validation must be JSON objects")
    discovery_ids = _patient_ids(discovery, "discovery")
    validation_ids = _patient_ids(validation, "validation")
    namespace = discovery.get("patient_id_namespace")
    if (
        not isinstance(namespace, str)
        or not namespace
        or namespace != validation.get("patient_id_namespace")
    ):
        raise ValueError("A shared, nonempty patient_id_namespace is required")
    overlap = discovery_ids & validation_ids
    if overlap:
        raise ValueError(f"Discovery and validation share {len(overlap)} patient IDs")
    candidates = validation.get("candidates")
    if not isinstance(candidates, list) or not candidates:
        raise ValueError("Validation requires a nonempty candidates list")
    if any(not isinstance(candidate, dict) for candidate in candidates):
        raise ValueError("Every candidate must be an object")
    ids = [candidate.get("candidate_id") for candidate in candidates]
    if any(not isinstance(x, str) or not x for x in ids) or len(ids) != len(set(ids)):
        raise ValueError("Candidate IDs must be unique nonempty strings")
    labels = [candidate.get("label") for candidate in candidates]
    if any(type(label) is not int or label not in (0, 1) for label in labels):
        raise ValueError("Every candidate requires an external binary label (0 or 1)")
    if len(set(labels)) != 2:
        raise ValueError("Validation candidates must include positive and negative labels")
    composite = _ranking_metrics(_scores(candidates, "composite_score"), k)
    recurrence = _ranking_metrics(_scores(candidates, "recurrence_score"), k)
    return {
        "patient_id_namespace": namespace,
        "discovery_patient_count": len(discovery_ids),
        "validation_patient_count": len(validation_ids),
        "patient_overlap_count": 0,
        "n_candidates": len(candidates),
        "n_positive_labels": sum(labels),
        "k": min(k, len(candidates)),
        "composite": composite,
        "recurrence": recurrence,
        "delta_average_precision": (
            composite["average_precision"] - recurrence["average_precision"]
        ),
        "delta_precision_at_k": composite["precision_at_k"] - recurrence["precision_at_k"],
        "interpretation": (
            "Descriptive ranking comparison on supplied external labels. Tied scores are "
            "evaluated as groups; precision at k allocates boundary ties fractionally. "
            "The supplied scores must have been frozen using discovery patients only; "
            "this metadata check cannot verify score derivation or external label quality."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("discovery", type=Path)
    parser.add_argument("validation", type=Path)
    parser.add_argument("--k", type=int, default=5)
    args = parser.parse_args()
    try:
        discovery = json.loads(args.discovery.read_text())
        validation = json.loads(args.validation.read_text())
        report = compare_heldout_rankings(discovery, validation, k=args.k)
    except (OSError, ValueError, TypeError, AttributeError) as exc:
        parser.error(str(exc))
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
