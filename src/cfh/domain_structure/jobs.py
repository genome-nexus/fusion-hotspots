"""Prediction jobs for the shortlisted chimeric domains and their parent windows.

For each shortlisted record three sequences are predicted with identical
settings: the chimeric domain with the same +/- ``FLANK`` window the sequence
step used (:attr:`ChimeraRecord.chimeric_domain_flanked`) and each parent's
domain +/- ``FLANK`` window (:meth:`ParentDomain.window`). A job is keyed by
its sequence, so a sequence shared by several records (a parent window used by
two fusions, or two records with one ``sequence_id``) is predicted once.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path

from cfh.domain_function.record import ChimeraRecord
from cfh.domain_structure.predictor import sequence_hash

ROLES = ("chimera", "parent5", "parent3")


@dataclass
class Job:
    job_id: str
    sequence: str
    roles: list[str] = field(default_factory=list)
    """``<role>:<label>`` for every record use, e.g. ``chimera:FLI1::ERG``."""


@dataclass(frozen=True)
class RecordJobs:
    record: ChimeraRecord
    shortlist: dict[str, str]
    chimera: str
    """Job id of the flanked chimeric domain."""
    parent5: str
    parent3: str
    chimera_start: int
    """1-based chimera position of the first predicted residue."""
    parent5_start: int
    parent3_start: int


def job_id_for(sequence: str) -> str:
    return f"job_{sequence_hash(sequence)}"


def read_shortlist(path: Path) -> list[dict[str, str]]:
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def build_jobs(
    shortlist: list[dict[str, str]], records: list[ChimeraRecord]
) -> tuple[list[RecordJobs], list[Job]]:
    by_id = {record.record_id: record for record in records}
    jobs: dict[str, Job] = {}

    def add(sequence: str, role: str) -> str:
        job_id = job_id_for(sequence)
        job = jobs.setdefault(job_id, Job(job_id, sequence))
        if role not in job.roles:
            job.roles.append(role)
        return job_id

    out: list[RecordJobs] = []
    for row in shortlist:
        record = by_id.get(row["record_id"])
        if record is None:
            raise KeyError(f"shortlisted record {row['record_id']} is not in the records file")
        if record.sequence_id != row["sequence_id"]:
            raise ValueError(f"{row['record_id']}: sequence_id differs from the records file")
        c_start, _ = record.window()
        f_start, f_end = record.five.window()
        t_start, t_end = record.three.window()
        out.append(
            RecordJobs(
                record=record,
                shortlist=row,
                chimera=add(record.chimeric_domain_flanked, f"chimera:{record.fusion}"),
                parent5=add(
                    record.five.protein[f_start - 1 : f_end], f"parent5:{record.five.gene}"
                ),
                parent3=add(
                    record.three.protein[t_start - 1 : t_end], f"parent3:{record.three.gene}"
                ),
                chimera_start=c_start,
                parent5_start=f_start,
                parent3_start=t_start,
            )
        )
    return out, sorted(jobs.values(), key=lambda job: (len(job.sequence), job.job_id))
