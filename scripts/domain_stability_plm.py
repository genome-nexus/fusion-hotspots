"""Protein-language-model scores for chimeric fused domains (PR 3 of 3).

Usage::

    # long, resumable: masked PLL with ESM-2 (priority records first, then a
    # seeded random order of the rest), cached per sequence
    python scripts/domain_stability_plm.py score --max-hours 3

    # fast: metrics, disagreement flags and committed outputs from the cache
    python scripts/domain_stability_plm.py analyze --out runs/domain_stability_plm_2026-11-13

The PLM env (torch + fair-esm) and the ESM-2 weights live outside the
repository (``$CFH_PLM_HOME``, default ``~/.cache/cfh_plm``; weights under
``$CFH_PLM_HOME/torch``). See the run's REPORT.md for the method.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

from cfh.domain_function.record import FLANK, read_records_json
from cfh.domain_plm import SCHEMA_VERSION
from cfh.domain_plm.cache import ScoreCache
from cfh.domain_plm.compare import (
    DISAGREEMENT_RULES,
    RANK_FRACTION,
    add_ranks,
    disagreements,
)
from cfh.domain_plm.metrics import junction_profile
from cfh.domain_plm.run import (
    NAMED_DNA_FUSIONS,
    NAMED_FUSIONS,
    ORDER_SEED,
    collect,
    metrics_row,
    read_tsv,
    score_in_order,
    scoring_order,
)
from cfh.domain_plm.scoring import (
    JUNCTION_HALF_WINDOW,
    ScoringError,
    WorkerScorer,
    check_against_metrics,
    plan_record,
)
from cfh.domain_plm.write import write_outputs

PLM_HOME = Path(os.environ.get("CFH_PLM_HOME", Path.home() / ".cache" / "cfh_plm"))
DEFAULT_MODEL = "esm2_t12_35M_UR50D"


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
    parser.add_argument("stage", choices=("score", "analyze"))
    parser.add_argument("--seq-run", type=Path, default=Path("runs/domain_function_seq_2026-11-11"))
    parser.add_argument(
        "--structure-run", type=Path, default=Path("runs/domain_structure_2026-11-13")
    )
    parser.add_argument(
        "--extra-records",
        type=Path,
        default=Path("runs/domain_stability_plm_2026-11-13/inputs"),
        help="directory with pcawg_tcga_chimera_records.json.gz + its metrics TSV",
    )
    parser.add_argument("--work", type=Path, default=Path("runs/.domain_stability_plm_work"))
    parser.add_argument("--out", type=Path)
    parser.add_argument("--model", default=DEFAULT_MODEL)
    parser.add_argument("--device", default="mps")
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--python", default=str(PLM_HOME / "env" / "bin" / "python"))
    parser.add_argument("--max-hours", type=float, help="stop scoring after this long")
    parser.add_argument(
        "--from-end",
        action="store_true",
        help="score only the random-order records, last first (a second worker, e.g. on CPU)",
    )
    args = parser.parse_args()

    records_path = args.seq_run / "chimera_records.json.gz"
    metrics_path = args.seq_run / "domain_function_metrics.tsv.gz"
    shortlist_path = args.seq_run / "shortlist_for_structure.tsv"
    structure_path = args.structure_run / "domain_structure_metrics.tsv.gz"
    extra_records_path = args.extra_records / "pcawg_tcga_chimera_records.json.gz"
    extra_metrics_path = args.extra_records / "pcawg_tcga_domain_function_metrics.tsv.gz"

    records = read_records_json(records_path)
    pr155 = {row["record_id"]: row for row in read_tsv(metrics_path)}
    extra = read_records_json(extra_records_path) if extra_records_path.exists() else []
    extra_rows = (
        {row["record_id"]: row for row in read_tsv(extra_metrics_path)}
        if extra_metrics_path.exists()
        else {}
    )
    extra_ids = {record.record_id for record in extra} - {r.record_id for r in records}
    tier_rows = {**extra_rows, **pr155}
    all_records = records + [record for record in extra if record.record_id in extra_ids]
    plans = [plan_record(record) for record in all_records]
    shortlist_ids = [row["record_id"] for row in read_tsv(shortlist_path)]
    ordered = scoring_order(plans, shortlist_ids, NAMED_FUSIONS + NAMED_DNA_FUSIONS, extra_ids)
    groups = {plan.record.record_id: group for group, plan in ordered}
    log(
        f"{len(records)} PR #155 records + {len(extra_ids)} PCAWG DNA records; "
        f"{sum(1 for p in plans if p.problem)} cannot be scored; "
        f"order: {dict(Counter(groups.values()))}"
    )
    env = {**os.environ, "TORCH_HOME": str(PLM_HOME / "torch")}

    scorer: WorkerScorer | None = None
    model_problem = ""
    try:
        scorer = WorkerScorer(args.python, args.model, args.device, args.batch_size, env)
    except (OSError, ScoringError, KeyError, ValueError) as error:
        model_problem = f"PLM model unavailable: {type(error).__name__}: {error}"
        log(model_problem)
    cache = ScoreCache(args.work / "scores", scorer.spec) if scorer else None
    failures_path = args.work / "failures.json"

    if args.stage == "score":
        if scorer is None or cache is None:
            return 2
        log(f"model {scorer.header}; cache {cache.dir}")
        deadline = time.time() + args.max_hours * 3600 if args.max_hours else None
        started = time.time()
        todo = ordered
        if args.from_end:
            todo = [item for item in ordered if item[0] == "random_order"][::-1]
        done, failures = score_in_order(todo, scorer, cache, log, deadline)
        scorer.close()
        tag = "_from_end" if args.from_end else ""
        failures_path = args.work / f"failures{tag}.json"
        previous = json.loads(failures_path.read_text()) if failures_path.exists() else {}
        args.work.mkdir(parents=True, exist_ok=True)
        failures_path.write_text(json.dumps({**previous, **failures}, indent=1))
        (args.work / f"score_status{tag}.json").write_text(
            json.dumps(
                {
                    "finished": time.strftime("%Y-%m-%dT%H:%M:%S"),
                    "records_scored_this_invocation": done,
                    "wall_seconds_this_invocation": round(time.time() - started, 1),
                    "failures": len(failures),
                    "model": scorer.header,
                },
                indent=1,
            )
        )
        log(f"score finished: {done} records this invocation, {len(failures)} failures")
        return 0

    if args.out is None:
        parser.error("analyze needs --out")
    if scorer is not None:
        scorer.close()
    failures = {}
    for path in sorted(args.work.glob("failures*.json")):
        failures.update(json.loads(path.read_text()))
    structure = {row["record_id"]: row for row in read_tsv(structure_path)}
    rows = []
    profiles = []
    for group, plan in ordered:
        record = plan.record
        rid = record.record_id
        tier_row = tier_rows.get(rid)
        mismatch = check_against_metrics(plan, tier_row)
        row = {
            "record_id": rid,
            "sequence_id": record.sequence_id,
            "fusion": record.fusion,
            "source": "pcawg_tcga_dna" if rid in extra_ids else "pr155",
            "scoring_group": group,
            "pfam_family": record.pfam_family,
            "junction_class": record.junction_class,
            "n_samples": record.n_samples,
            "chimeric_domain_span": f"{record.domain_start}-{record.domain_end}",
            "chimeric_domain_length": len(record.chimeric_domain),
            "junction_index": record.junction_index,
            "pilot_indel_residues": ""
            if record.pilot_indel_residues is None
            else record.pilot_indel_residues,
            "pr155_tier": (tier_row or {}).get("tier", "") or "missing",
            "pr160_fold_flag": structure.get(rid, {}).get("structure_flag", "")
            or "not_in_structure_run",
            **metrics_row(plan, cache, failures, mismatch, model_problem),
        }
        rows.append(row)
        if row["plm_status"] == "scored" and cache is not None:
            scores, _ = collect(plan, cache)
            for position, residue, value in junction_profile(plan, scores["chimera"]):
                profiles.append((rid, record.fusion, position, residue, value))
    add_ranks(rows)
    for row in rows:
        row["disagreements"] = ";".join(disagreements(row))
    spec = scorer.spec if scorer else None
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "code_commit": git_commit(),
        "generated": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "inputs": {
            name: {"path": str(path), "sha256": sha256(path)}
            for name, path in (
                ("pr155_records", records_path),
                ("pr155_metrics", metrics_path),
                ("pr155_shortlist", shortlist_path),
                ("pr160_structure_metrics", structure_path),
                ("pcawg_tcga_records", extra_records_path),
                ("pcawg_tcga_metrics", extra_metrics_path),
            )
            if path.exists()
        },
        "pcawg_tcga_source": {
            "branch": "origin/polly/pcawg-tcga-sv-screen",
            "commit": "fb0cfcf7eec9880b2e3b9a59c89ffb2295d2dc1b",
            "path": "runs/pcawg_tcga_sv_screen_2026-11-13/chimera_records.json.gz",
            "blob": "488c74df28a9162a424b03b3c87bda90385a4195",
        },
        "model": {
            "spec": None if spec is None else spec.__dict__,
            "spec_id": None if spec is None else spec.spec_id(),
            "header": None if scorer is None else scorer.header,
            "problem": model_problem,
            "weights_url": f"https://dl.fbaipublicfiles.com/fair-esm/models/{args.model}.pt",
        },
        "method": {
            "score": "masked pseudo-log-likelihood: one residue masked per forward pass, "
            "log p(true residue), natural log",
            "windows": {"flank": FLANK, "junction_half_window": JUNCTION_HALF_WINDOW},
            "rank_fraction": RANK_FRACTION,
            "disagreement_rules": DISAGREEMENT_RULES,
            "random_order_seed": ORDER_SEED,
        },
        "scoring_sessions": {
            path.name: json.loads(path.read_text())
            for path in sorted(args.work.glob("score_status*.json"))
        },
        "scoring_failures": failures,
    }
    out = args.out
    write_outputs(out, rows, profiles, manifest, NAMED_FUSIONS + NAMED_DNA_FUSIONS)
    log(f"wrote {out}: {Counter(r['plm_status'] for r in rows)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
