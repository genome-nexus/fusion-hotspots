"""Driver stages: MSAs, predictions (checkpointed per job), metrics, flags and outputs."""

from __future__ import annotations

import csv
import gzip
import json
import shutil
from collections import Counter
from pathlib import Path
from typing import Any, Callable

from cfh.domain_function.record import ChimeraRecord
from cfh.domain_structure.afdb_models import AfdbCoordinates
from cfh.domain_structure.flags import assign_flag
from cfh.domain_structure.jobs import Job, RecordJobs
from cfh.domain_structure.metrics import AfdbModel, record_metrics
from cfh.domain_structure.predictor import (
    DiskSpaceError,
    MsaCache,
    Prediction,
    PredictionFailure,
    Predictor,
)

Log = Callable[[str], None]


def predict_all(
    jobs: list[Job], predictor: Predictor, msas: MsaCache, log: Log
) -> tuple[dict[str, Prediction], dict[str, str]]:
    """Fetch missing MSAs (batched), then predict every job; resumable at job granularity.

    A :class:`DiskSpaceError` propagates (the run stops cleanly); every other
    failure is recorded per job and the run moves on.
    """
    predictions: dict[str, Prediction] = {}
    failures: dict[str, str] = {}
    todo = [job for job in jobs if predictor.cached(job.job_id, job.sequence) is None]
    msa_errors = msas.fetch([job.sequence for job in todo], log=log)
    for job in jobs:
        if job.sequence in msa_errors:
            failures[job.job_id] = msa_errors[job.sequence]
            predictor._fail(predictor.job_dir(job.job_id), msa_errors[job.sequence])
            continue
        try:
            predictions[job.job_id] = predictor.predict(job.job_id, job.sequence, log=log)
        except DiskSpaceError:
            raise
        except PredictionFailure as error:
            failures[job.job_id] = str(error)
            log(f"FAILED {job.job_id}: {error}")
    return predictions, failures


def load_predictions(
    jobs: list[Job], predictor: Predictor
) -> tuple[dict[str, Prediction], dict[str, str]]:
    """Validated cached predictions; anything else gets its recorded failure reason."""
    predictions: dict[str, Prediction] = {}
    failures: dict[str, str] = {}
    for job in jobs:
        cached = predictor.cached(job.job_id, job.sequence)
        if cached is not None:
            predictions[job.job_id] = cached
            continue
        failure = predictor.job_dir(job.job_id) / "failure.json"
        failures[job.job_id] = (
            json.loads(failure.read_text())["reason"] if failure.exists() else "not_run"
        )
    return predictions, failures


def afdb_models(
    record_jobs: list[RecordJobs], client: AfdbCoordinates
) -> dict[str, AfdbModel | str]:
    out: dict[str, AfdbModel | str] = {}
    for jobs in record_jobs:
        for column in ("uniprot_5p", "uniprot_3p"):
            accession = jobs.shortlist.get(column, "")
            if accession not in out:
                out[accession] = client.model(accession)
    return out


ID_COLUMNS = (
    "record_id",
    "sequence_id",
    "fusion",
    "pfam_family",
    "junction_class",
    "cohorts",
    "n_samples",
    "shortlist_rank",
    "shortlist_reason",
    "tier",
    "chimeric_domain_span",
    "chimeric_domain_length",
    "five_prime_domain",
    "three_prime_domain",
    "junction_index",
    "aligned_inserted_sequence",
    "aligned_deleted_columns",
    "uniprot_5p",
    "uniprot_3p",
)


FLAG_COLUMNS = ("structure_flag", "structure_flag_rules", "structure_flag_missing")


def assess(
    record_jobs: list[RecordJobs],
    predictions: dict[str, Prediction],
    failures: dict[str, str],
    afdb: dict[str, AfdbModel | str],
) -> list[dict[str, Any]]:
    rows = []
    for jobs in record_jobs:
        row: dict[str, Any] = {key: jobs.shortlist.get(key, "") for key in ID_COLUMNS}
        row["pr1_tier"] = row.pop("tier")
        models = {
            "five": afdb.get(jobs.shortlist.get("uniprot_5p", ""), "no UniProt accession"),
            "three": afdb.get(jobs.shortlist.get("uniprot_3p", ""), "no UniProt accession"),
        }
        try:
            row.update(record_metrics(jobs, predictions, failures, models))
        except (ValueError, KeyError, IndexError, PredictionFailure) as error:
            row["metrics_error"] = f"{type(error).__name__}: {error}"
            for role in ("chimera", "parent5", "parent3"):
                row[f"{role}_prediction_status"] = f"metrics_error: {error}"
        flag, rules, missing = assign_flag(row)
        row["structure_flag"] = flag
        row["structure_flag_rules"] = ";".join(rules)
        row["structure_flag_missing"] = ";".join(missing)
        rows.append(row)
    return rows


def _cell(value: Any) -> Any:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    return value


def write_tsv(rows: list[dict[str, Any]], path: Path) -> None:
    columns = ["pr1_tier" if c == "tier" else c for c in ID_COLUMNS] + list(FLAG_COLUMNS)
    for row in sorted(rows, key=len, reverse=True):  # the most complete row sets the order
        for key in row:
            if key not in columns:
                columns.append(key)
    columns = [c for c in columns if any(c in row for row in rows)]
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _cell(row.get(key)) for key in columns})


def export_models(
    jobs: list[Job], predictions: dict[str, Prediction], out_dir: Path
) -> dict[str, dict[str, Any]]:
    """Top-ranked model (PDB, gzip) and compact pLDDT/PAE JSON (gzip) per prediction.

    pLDDT keeps 2 decimals and PAE 1 decimal (metrics are computed before rounding).
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    index: dict[str, dict[str, Any]] = {}
    for job in jobs:
        prediction = predictions.get(job.job_id)
        if prediction is None:
            continue
        shutil.copyfile(prediction.pdb_path, out_dir / f"{job.job_id}.pdb.gz")
        scores = {
            "job_id": job.job_id,
            "sequence": job.sequence,
            "model_name": prediction.model_name,
            "ptm": round(prediction.ptm, 4),
            "plddt": [round(x, 2) for x in prediction.plddt],
            "pae": [[round(x, 1) for x in row] for row in prediction.pae],
        }
        with gzip.open(out_dir / f"{job.job_id}.scores.json.gz", "wt") as handle:
            json.dump(scores, handle, separators=(",", ":"))
        index[job.job_id] = {
            "roles": job.roles,
            "length": len(job.sequence),
            "wall_seconds": prediction.wall_seconds,
            "model_name": prediction.model_name,
        }
    return index


def flag_counts(rows: list[dict[str, Any]], records: list[ChimeraRecord]) -> dict[str, Any]:
    by_id = {r.record_id: r for r in records}
    overall = Counter(row["structure_flag"] for row in rows)
    per_cohort: dict[str, Counter] = {}
    for row in rows:
        record = by_id.get(row["record_id"])
        for cohort in record.cohorts if record else []:
            per_cohort.setdefault(cohort, Counter())[row["structure_flag"]] += 1
    by_tier: dict[str, Counter] = {}
    for row in rows:
        by_tier.setdefault(row["pr1_tier"], Counter())[row["structure_flag"]] += 1
    rules = Counter(rule for row in rows for rule in row["structure_flag_rules"].split(";") if rule)
    return {
        "records": len(rows),
        "by_flag": dict(overall),
        "by_cohort": {k: dict(v) for k, v in sorted(per_cohort.items())},
        "by_pr1_tier": {k: dict(v) for k, v in sorted(by_tier.items())},
        "rules_fired": dict(rules),
    }


TABLE_COLUMNS = (
    ("fusion", "fusion"),
    ("pfam_family", "Pfam"),
    ("pr1_tier", "PR1 tier"),
    ("structure_flag", "flag"),
    ("structure_flag_rules", "rules"),
    ("chimera_domain_plddt", "dom pLDDT C"),
    ("parent5_domain_plddt", "5'"),
    ("parent3_domain_plddt", "3'"),
    ("delta_domain_plddt_vs_better_parent", "Δdom"),
    ("chimera_junction_plddt_mean", "jx mean C"),
    ("chimera_junction_plddt_min", "jx min C"),
    ("delta_junction_plddt_mean_vs_better_parent", "Δjx"),
    ("chimera_interhalf_pae", "half PAE C"),
    ("parent5_interhalf_pae", "5'"),
    ("parent3_interhalf_pae", "3'"),
    ("chimera_ptm", "pTM C"),
    ("vs_5p_tm", "TM vs 5'"),
    ("vs_3p_tm", "TM vs 3'"),
    ("vs_5p_local_rmsd", "jx RMSD 5'"),
    ("vs_3p_local_rmsd", "jx RMSD 3'"),
    ("vs_afdb_5p_tm", "TM AFDB 5'"),
    ("vs_afdb_3p_tm", "TM AFDB 3'"),
)


def report_table(rows: list[dict[str, Any]]) -> str:
    """Markdown table of the main per-record measures (for REPORT.md)."""
    header = "| " + " | ".join(label for _, label in TABLE_COLUMNS) + " |"
    lines = [header, "|" + "---|" * len(TABLE_COLUMNS)]
    for row in rows:
        cells = []
        for key, _ in TABLE_COLUMNS:
            value = row.get(key)
            cells.append("" if value is None else str(value).replace("|", "/"))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines) + "\n"
