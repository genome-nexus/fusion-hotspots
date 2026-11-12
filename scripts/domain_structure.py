"""Predict and compare structures of the shortlisted chimeric domains (PR 2 of 3).

Usage::

    # long, resumable: MSAs (public ColabFold server, batched, cached) + predictions
    python scripts/domain_structure.py predict \\
        --seq-run runs/domain_function_seq_2026-11-11 --work runs/.domain_structure_work

    # fast: metrics, flags and committed outputs from the cached predictions
    python scripts/domain_structure.py analyze \\
        --seq-run runs/domain_function_seq_2026-11-11 --work runs/.domain_structure_work \\
        --out runs/domain_structure_2026-11-12

``--config predictor.json`` overrides any :class:`PredictorConfig` field (e.g. a
GPU host's ``colabfold_command``). The ColabFold env and weights live outside
the repository (``$CFH_COLABFOLD_HOME``, default ``~/.cache/cfh_colabfold``);
see ``docs/domain_structure.md``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import time
from dataclasses import asdict
from pathlib import Path

from cfh.domain_function.record import FLANK, read_records_json
from cfh.domain_structure import SCHEMA_VERSION
from cfh.domain_structure.afdb_models import AfdbCoordinates
from cfh.domain_structure.flags import RULES
from cfh.domain_structure.jobs import build_jobs, read_shortlist
from cfh.domain_structure.metrics import JUNCTION_FLANK, LOCAL_RADIUS, MIN_MAPPED_PAIRS
from cfh.domain_structure.predictor import (
    DiskSpaceError,
    MsaCache,
    Predictor,
    PredictorConfig,
    check_disk,
    colabfold_version,
    weights_manifest,
)
from cfh.domain_structure.run import (
    afdb_models,
    assess,
    export_models,
    flag_counts,
    load_predictions,
    predict_all,
    write_tsv,
)


def log(message: str) -> None:
    print(time.strftime("%Y-%m-%d %H:%M:%S"), message, flush=True)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def git_commit() -> str:
    try:
        return subprocess.run(
            ["git", "rev-parse", "HEAD"], capture_output=True, text=True, check=True
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("stage", choices=("predict", "analyze"))
    parser.add_argument("--seq-run", type=Path, required=True)
    parser.add_argument("--work", type=Path, default=Path("runs/.domain_structure_work"))
    parser.add_argument("--out", type=Path)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--limit", type=int, help="only the first N jobs (smoke test)")
    args = parser.parse_args()

    config = PredictorConfig.load(args.config)
    shortlist_path = args.seq_run / "shortlist_for_structure.tsv"
    records_path = args.seq_run / "chimera_records.json.gz"
    records = read_records_json(records_path)
    record_jobs, jobs = build_jobs(read_shortlist(shortlist_path), records)
    if args.limit:
        jobs = jobs[: args.limit]
    args.work.mkdir(parents=True, exist_ok=True)
    msas = MsaCache(args.work, config)
    predictor = Predictor(args.work / "predictions", config, msas)
    log(
        f"{len(record_jobs)} records, {len(jobs)} distinct sequences "
        f"({sum(len(j.sequence) for j in jobs)} residues); settings {config.settings_id()}"
    )

    if args.stage == "predict":
        started = time.time()
        try:
            check_disk(args.work, config.min_free_disk_gb)
            predictions, failures = predict_all(jobs, predictor, msas, log)
        except DiskSpaceError as error:
            log(f"ABORT: {error}")
            return 2
        status = {
            "finished": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "wall_seconds_this_invocation": round(time.time() - started, 1),
            "predicted": len(predictions),
            "failed": failures,
        }
        (args.work / "predict_status.json").write_text(json.dumps(status, indent=1))
        log(f"predict finished: {len(predictions)} ok, {len(failures)} failed")
        return 0

    if args.out is None:
        parser.error("analyze needs --out")
    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    predictions, failures = load_predictions(jobs, predictor)
    afdb = afdb_models(record_jobs, AfdbCoordinates(args.work / "afdb"))
    rows = assess(record_jobs, predictions, failures, afdb)
    write_tsv(rows, out / "domain_structure_metrics.tsv.gz")
    index = export_models(jobs, predictions, out / "models")
    counts = flag_counts(rows, records)
    (out / "flag_counts.json").write_text(json.dumps(counts, indent=1) + "\n")
    walls = [p.wall_seconds for p in predictions.values() if p.wall_seconds is not None]
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "code_commit": git_commit(),
        "generated": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "inputs": {
            "shortlist": {"path": str(shortlist_path), "sha256": sha256(shortlist_path)},
            "records": {"path": str(records_path), "sha256": sha256(records_path)},
        },
        "predictor": {
            "version": colabfold_version(config),
            "settings_id": config.settings_id(),
            "settings": config.prediction_settings(),
            "msa_server": config.host_url,
            "msa_user_agent": config.user_agent,
            "msa_batch_size": config.msa_batch_size,
            "weights": {
                "source": "https://storage.googleapis.com/alphafold/"
                "alphafold_params_2021-07-14.tar",
                "sha256": weights_manifest(config),
            },
            "relax": False,
            "platform": f"{sys.platform}; CPU JAX",
        },
        "windows": {"flank": FLANK, "junction_flank": JUNCTION_FLANK},
        "comparison": {"local_radius_A": LOCAL_RADIUS, "min_mapped_pairs": MIN_MAPPED_PAIRS},
        "flag_rules": RULES,
        "jobs": {
            "distinct_sequences": len(jobs),
            "residues": sum(len(j.sequence) for j in jobs),
            "predicted": len(predictions),
            "failed": failures,
        },
        "wall_time": {
            "prediction_seconds_sum": round(sum(walls), 1),
            "prediction_hours_sum": round(sum(walls) / 3600, 2),
            "per_job_seconds": {k: v["wall_seconds"] for k, v in index.items()},
        },
        "models": index,
        "flag_counts": counts,
        "config": {k: v for k, v in asdict(config).items() if k not in ("user_agent",)},
    }
    (out / "run_manifest.json").write_text(json.dumps(manifest, indent=1) + "\n")
    log(f"wrote {out}: {counts['by_flag']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
