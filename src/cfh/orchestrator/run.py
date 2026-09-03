"""Run registered algorithms concurrently against a common input snapshot."""

from __future__ import annotations

import hashlib
import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from time import perf_counter
from typing import Any, Sequence

from cfh.algorithms.registry import get
from cfh.genes.registry import GeneConfig
from cfh.model.algorithm_result import AlgorithmResult
from cfh.model.fusion_event import FusionEvent
from cfh.model.fusion_feature import FusionFeature


def _input_fingerprint(events: Sequence[FusionEvent], features: Sequence[FusionFeature]) -> str:
    """Hash a canonical representation of the input tables."""
    payload = {
        "events": [event.model_dump(mode="json") for event in events],
        "features": [feature.model_dump(mode="json") for feature in features],
    }
    serialized = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(serialized.encode()).hexdigest()


def _params_for(
    name: str,
    params: dict[str, Any] | None,
    algorithm_names: Sequence[str] | None = None,
) -> dict[str, Any]:
    """Resolve optional per-algorithm parameters without sharing mutable state."""
    if not params:
        return {}
    named_params = params.get(name)
    if isinstance(named_params, dict):
        return dict(named_params)
    if algorithm_names and any(
        isinstance(params.get(algorithm), dict) for algorithm in algorithm_names
    ):
        return {}
    if all(isinstance(value, dict) for value in params.values()):
        return {}
    return dict(params)


def _run_one(
    name: str,
    events: list[FusionEvent],
    features: list[FusionFeature],
    gene_config: GeneConfig,
    params: dict[str, Any],
    input_fingerprint: str,
) -> AlgorithmResult:
    """Execute one plugin and convert exceptions to structured failures."""
    started_at = perf_counter()
    try:
        result = get(name)().run(events, features, gene_config, params)
        summary = dict(result.Summary or {})
        summary["Runtime_seconds"] = perf_counter() - started_at
        return result.model_copy(
            update={
                "Parameters": result.Parameters if result.Parameters is not None else params,
                "Summary": summary,
                "Warnings": result.Warnings or [],
                "Created_at": result.Created_at or datetime.now(timezone.utc),
                "Input_fingerprint": input_fingerprint,
            }
        )
    except Exception as exc:
        return AlgorithmResult(
            Algorithm=name,
            Parameters=params,
            Summary={"Runtime_seconds": perf_counter() - started_at},
            Tables={},
            Warnings=[f"Algorithm failed: {type(exc).__name__}: {exc}"],
            Created_at=datetime.now(timezone.utc),
            Input_fingerprint=input_fingerprint,
        )


def run_algorithms(
    algorithm_names: Sequence[str],
    events: list[FusionEvent],
    features: list[FusionFeature],
    gene_config: GeneConfig,
    params: dict[str, Any] | None = None,
    *,
    max_workers: int | None = None,
) -> list[AlgorithmResult]:
    """Run registered algorithms concurrently and return results in request order.

    ``params`` may contain per-algorithm dictionaries keyed by algorithm name.
    ``results_to_json`` serializes this return value as a combined JSON array.
    Exceptions from an individual plugin become that plugin's warning result,
    allowing independent algorithms in the same run to complete successfully.
    """
    input_fingerprint = _input_fingerprint(events, features)
    if not algorithm_names:
        return []

    results: list[AlgorithmResult | None] = [None] * len(algorithm_names)
    worker_count = max_workers or len(algorithm_names)
    with ThreadPoolExecutor(max_workers=worker_count) as executor:
        futures = {
            executor.submit(
                _run_one,
                name,
                [event.model_copy(deep=True) for event in events],
                [feature.model_copy(deep=True) for feature in features],
                gene_config.model_copy(deep=True),
                _params_for(name, params, algorithm_names),
                input_fingerprint,
            ): index
            for index, name in enumerate(algorithm_names)
        }
        for future in as_completed(futures):
            results[futures[future]] = future.result()

    return [result for result in results if result is not None]


def results_to_json(results: Sequence[AlgorithmResult]) -> str:
    """Serialize combined algorithm results into a JSON array for emission."""
    return json.dumps([result.model_dump(mode="json") for result in results])
