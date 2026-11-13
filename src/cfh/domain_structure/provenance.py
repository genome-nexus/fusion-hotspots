"""One-off: write INFERRED MSA provenance into the 2026-11-12 legacy checkpoints.

The 71 predictions of the committed 2026-11-13 run were made (code 3da242a)
before checkpoints recorded their MSA source, under settings id
``set_e20d3cb4e3e2``, whose key did not include the MSA server. Nothing per
job records which server produced each a3m, so their provenance can only be
*inferred* from run-level evidence (:data:`EVIDENCE_2026_11_12`). The user
chose to keep them with that provenance documented as inferred rather than
re-predict.

:func:`record_inferred_provenance` is the only way a legacy checkpoint
becomes usable. It is never called by ``predict`` or ``analyze``; it runs
from ``scripts/domain_structure_inferred_provenance.py``. It refuses unless:

* the config's MSA source (mode, server, cache) is exactly the one the
  evidence names -- it cannot relabel to any other server;
* the work directory holds exactly one MSA cache and it is the evidence's;
* the legacy settings id is the config's settings id without ``host_url``
  (every other prediction setting matches);
* every legacy job re-validates (scores, model) and its query's a3m exists in
  that cache.

Each checkpoint then gets ``msa_source``, ``msa_provenance: inferred``, the
evidence summary and the sha256 of the a3m now in the cache, and is moved to
the config's settings id. Nothing is re-predicted or deleted.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Callable

from cfh.domain_structure.predictor import (
    MsaCache,
    PredictionFailure,
    PredictorConfig,
    load_prediction,
)


@dataclass(frozen=True)
class InferredEvidence:
    legacy_settings_id: str
    msa_mode: str
    host_url: str
    msa_cache: str
    summary: tuple[str, ...]

    def msa_source(self) -> dict[str, str]:
        return {"msa_mode": self.msa_mode, "host_url": self.host_url, "msa_cache": self.msa_cache}


EVIDENCE_2026_11_12 = InferredEvidence(
    legacy_settings_id="set_e20d3cb4e3e2",
    msa_mode="mmseqs2_uniref_env",
    host_url="https://api.colabfold.com",
    msa_cache="msa_c3444ece",
    summary=(
        "single predict launch without --config (session transcript), so the default "
        "PredictorConfig applied; default host_url at launch commit 3da242a was "
        "https://api.colabfold.com",
        "single invocation in predict.log (2026-11-12 16:23:38 to 2026-11-13 05:03:03)",
        "single MSA cache msa_c3444ece = sha256('mmseqs2_uniref_env|https://api.colabfold.com')"
        "[:8], filled 2026-11-12 16:23:37-17:16:51, before the first prediction (17:17)",
        "the FLI1::ERG chimera a3m was seeded from the benchmark fetch, whose request "
        "(bench/req.json) names host_url https://api.colabfold.com",
    ),
)
"""Run-level evidence for the committed 2026-11-13 run (see docs/domain_structure.md)."""


def _settings_id_without_host(config: PredictorConfig) -> str:
    settings = {k: v for k, v in config.prediction_settings().items() if k != "host_url"}
    text = json.dumps(settings, sort_keys=True)
    return "set_" + hashlib.sha256(text.encode()).hexdigest()[:12]


def record_inferred_provenance(
    work: Path,
    config: PredictorConfig,
    evidence: InferredEvidence = EVIDENCE_2026_11_12,
    log: Callable[[str], None] = print,
) -> int:
    """Mark the legacy checkpoints as inferred-provenance and move them to the config's key.

    Returns the number of checkpoints written; raises PredictionFailure (and
    changes nothing) when any check in the module docstring fails.
    """
    if config.msa_source() != evidence.msa_source():
        raise PredictionFailure(
            f"provenance: config MSA source {config.msa_source()} is not the evidence's "
            f"{evidence.msa_source()}"
        )
    if _settings_id_without_host(config) != evidence.legacy_settings_id:
        raise PredictionFailure(
            "provenance: the config's prediction settings do not match the legacy checkpoints"
        )
    caches = sorted(path.name for path in work.glob("msa_*") if path.is_dir())
    if caches != [evidence.msa_cache]:
        raise PredictionFailure(
            f"provenance: MSA caches present {caches}; the evidence requires exactly "
            f"[{evidence.msa_cache!r}]"
        )
    msas = MsaCache(work, config)
    assert msas.root.name == evidence.msa_cache
    predictions_root = work / "predictions"
    legacy = predictions_root / evidence.legacy_settings_id
    target = predictions_root / config.settings_id()
    jobs = sorted(d for d in legacy.glob("job_*") if (d / "result.json").exists())
    if not jobs:
        raise PredictionFailure(f"provenance: no checkpoints under {legacy}")
    if any(target.glob("job_*")):
        raise PredictionFailure(f"provenance: {target.name} already has checkpoints")
    updates: list[tuple[Path, dict]] = []
    for directory in jobs:
        meta = json.loads((directory / "result.json").read_text())
        if meta.get("settings_id") != legacy.name:
            raise PredictionFailure(f"provenance: {directory.name} records another settings id")
        if meta.get("msa_source") or meta.get("msa_provenance"):
            raise PredictionFailure(f"provenance: {directory.name} already has a record")
        sequence = meta.get("sequence", "")
        load_prediction(directory, sequence)
        a3m = msas.get(sequence)
        if a3m is None:
            raise PredictionFailure(
                f"provenance: no a3m for {directory.name} in {evidence.msa_cache}"
            )
        meta.update(
            settings_id=target.name,
            msa_source=evidence.msa_source(),
            msa_provenance="inferred",
            msa_provenance_evidence={
                "legacy_settings_id": legacy.name,
                "summary": list(evidence.summary),
                "a3m_sha256_note": "sha256 of the a3m in the cache when this record was "
                "written, not recorded at prediction time",
            },
            a3m_sha256=hashlib.sha256(a3m.encode()).hexdigest(),
        )
        updates.append((directory, meta))
    target.mkdir(parents=True, exist_ok=True)
    for directory, meta in updates:
        (directory / "result.json").write_text(json.dumps(meta, indent=1))
        directory.rename(target / directory.name)
    if not any(legacy.iterdir()):
        legacy.rmdir()
    log(
        f"inferred provenance written for {len(updates)} checkpoints: "
        f"{legacy.name} -> {target.name} ({asdict(evidence)['host_url']})"
    )
    return len(updates)
