"""Run the sequence-level functional-domain checks over a junction-screen table."""

from __future__ import annotations

import csv
import gzip
import json
import math
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping

from pyhmmer.plan7 import HMM

from cfh.domain_function import SCHEMA_VERSION
from cfh.domain_function.afdb import AfdbClient, ParentStructure
from cfh.domain_function.assess import assess_record
from cfh.domain_function.metrics import NOT_ASSESSED
from cfh.domain_function.pfam_hmm import HmmUnavailable, fetch_hmm, load_hmm
from cfh.domain_function.record import (
    ARTEFACT_FLAG_COLUMNS,
    FLANK,
    ChimeraRecord,
    record_to_dict,
)
from cfh.domain_function.tiers import TIER_ORDER
from cfh.domain_function.uniprot import UniProtEntry

NAMED_CASES = (
    "FLI1::ERG",
    "STAT5B::STAT5A",
    "GPC4::GPC3",
    "GPC6::GPC5",
    "GPC3::GPC4",
    "LPO::MPO",
    "PIM3::ULK1",
)
SHORTLIST_SIZE = 25

Progress = Callable[[str], None]


def read_screen_rows(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_hmms(
    families: Iterable[str], cache_dir: Path, progress: Progress | None = None
) -> tuple[dict[str, HMM | None], dict[str, str]]:
    hmms: dict[str, HMM | None] = {}
    problems: dict[str, str] = {}
    for index, family in enumerate(sorted(set(families)), start=1):
        try:
            hmms[family] = load_hmm(fetch_hmm(family, cache_dir))
        except HmmUnavailable as error:
            hmms[family] = None
            problems[family] = str(error)
        if progress and index % 50 == 0:
            progress(f"HMMs {index}")
    return hmms, problems


def load_structures(
    accessions: Iterable[str], cache_dir: Path, workers: int = 2, progress: Progress | None = None
) -> tuple[dict[str, ParentStructure | None], dict[str, str]]:
    pending = sorted(set(accessions))
    structures: dict[str, ParentStructure | None] = {}
    problems: dict[str, str] = {}

    def fetch(accession: str) -> tuple[str, ParentStructure | None, str]:
        try:
            return accession, AfdbClient(cache_dir).structure(accession), ""
        except Exception as error:  # noqa: BLE001 - recorded per accession, fail closed
            return accession, None, f"{type(error).__name__}: {error}"

    with ThreadPoolExecutor(max_workers=workers) as pool:
        for index, (accession, structure, problem) in enumerate(pool.map(fetch, pending), start=1):
            structures[accession] = structure
            if problem:
                problems[accession] = problem
            if progress and index % 100 == 0:
                progress(f"AFDB {index}/{len(pending)}")
    return structures, problems


def assess_all(
    records: list[ChimeraRecord],
    hmms: Mapping[str, HMM | None],
    entries: Mapping[str, UniProtEntry | None],
    structures: Mapping[str, ParentStructure | None],
    progress: Progress | None = None,
) -> list[dict[str, Any]]:
    out = []
    for index, record in enumerate(records, start=1):
        record_entries = {
            "five": entries.get(record.five.gene),
            "three": entries.get(record.three.gene),
        }
        record_structures = {
            which: structures.get(entry.accession) if entry else None
            for which, entry in record_entries.items()
        }
        metrics = assess_record(
            record, hmms.get(record.pfam_family), record_entries, record_structures
        )
        out.append({**summary_columns(record), **metrics})
        if progress and index % 250 == 0:
            progress(f"assessed {index}/{len(records)}")
    return out


def summary_columns(record: ChimeraRecord) -> dict[str, Any]:
    """Record fields carried into the flat metrics table."""
    flags = {
        f"{column}_counts": ",".join(
            f"{k}:{v}" for k, v in sorted(record.artefact_flag_counts.get(column, {}).items())
        )
        for column in ARTEFACT_FLAG_COLUMNS
    }
    return {
        "record_id": record.record_id,
        "fusion": record.fusion,
        "fusions": ";".join(record.fusions),
        "pfam_family": record.pfam_family,
        "junction_class": record.junction_class,
        "cohorts": ";".join(record.cohorts),
        "rows_per_cohort": ";".join(f"{k}:{v}" for k, v in sorted(record.rows_per_cohort.items())),
        "n_rows": record.n_rows,
        "n_samples": record.n_samples,
        "max_split_reads": record.max_split_reads,
        "total_split_reads": record.total_split_reads,
        "max_paired_end_reads": record.max_paired_end_reads,
        "artefact_clean_rows": record.artefact_clean_rows,
        "artefact_core_clean_rows": record.artefact_core_clean_rows,
        **flags,
        "five_prime_transcript": record.five.transcript,
        "three_prime_transcript": record.three.transcript,
        "five_prime_domain": f"{record.five.domain_start}-{record.five.domain_end}",
        "three_prime_domain": f"{record.three.domain_start}-{record.three.domain_end}",
        "five_prime_last_residue": record.five.junction_residue,
        "three_prime_first_residue": record.three.junction_residue,
        "hybrid_residue": record.hybrid_residue,
        "chimeric_domain_span": f"{record.domain_start}-{record.domain_end}",
        "junction_index": record.junction_index,
        "pilot_indel_residues": record.pilot_indel_residues,
        "pilot_indel_sequence": record.pilot_indel_sequence,
        "chimeric_domain_length": record.domain_end - record.domain_start + 1,
        "chimeric_domain_flanked": record.chimeric_domain_flanked,
    }


def schema_e_max(row: Mapping[str, Any]) -> float:
    values = [row.get("schema_E_5p"), row.get("schema_E_3p")]
    numbers = [v for v in values if isinstance(v, int) and not isinstance(v, bool)]
    return float(max(numbers)) if numbers else math.inf


def _flagged_columns(row: Mapping[str, Any]) -> list[str]:
    """Artefact flags with at least one ``flagged`` row for this record."""
    flagged = []
    for column in ARTEFACT_FLAG_COLUMNS:
        counts = str(row.get(f"{column}_counts") or "")
        if any(item.startswith("flagged:") for item in counts.split(",")):
            flagged.append(column.removesuffix("_flag"))
    return flagged


def _as_int(value: Any) -> int:
    try:
        return int(value or 0)
    except (TypeError, ValueError):
        return 0


def shortlist(rows: list[dict[str, Any]], size: int = SHORTLIST_SIZE) -> list[dict[str, Any]]:
    """Structure-prediction shortlist: one record per fusion and Pfam family.

    Order within each stage: tier, then lowest SCHEMA E (max over the parent
    structures; missing last), then highest GA margin, then max split reads and
    sample count. Stages, up to ``size`` rows in total:

    1. ``ranked``: records with at least one row passing every artefact flag;
    2. ``fill``: records with a row passing every flag except the near-universal
       context flags (:data:`~cfh.domain_function.record.CONTEXT_FLAG_COLUMNS`);
       the flags that fired are named in ``shortlist_reason``.

    The named cases are then appended (best record per family) when absent,
    marked ``named_case`` with their flagged columns.
    """

    def key(row: Mapping[str, Any]) -> tuple[Any, ...]:
        margin = row.get("chimera_ga_margin")
        return (
            TIER_ORDER.index(row["tier"]),
            schema_e_max(row),
            -(margin if isinstance(margin, (int, float)) else -math.inf),
            -_as_int(row.get("max_split_reads")),
            -_as_int(row.get("n_samples")),
            row["record_id"],
        )

    def best_per_group(candidates: Iterable[dict[str, Any]]) -> list[dict[str, Any]]:
        best: dict[tuple[str, str], dict[str, Any]] = {}
        for row in sorted(candidates, key=key):
            best.setdefault((row["fusion"], row["pfam_family"]), row)
        return sorted(best.values(), key=key)

    def note(stage: str, row: Mapping[str, Any]) -> str:
        flagged = _flagged_columns(row)
        return f"{stage};flagged:{','.join(flagged)}" if flagged else stage

    chosen: list[dict[str, Any]] = []
    groups: set[tuple[str, str]] = set()
    stages = (
        ("ranked", lambda row: _as_int(row.get("artefact_clean_rows")) > 0),
        ("fill", lambda row: _as_int(row.get("artefact_core_clean_rows")) > 0),
    )
    for stage, eligible in stages:
        for row in best_per_group(r for r in rows if eligible(r)):
            group = (row["fusion"], row["pfam_family"])
            if len(chosen) >= size or group in groups:
                continue
            chosen.append(dict(row, shortlist_reason=note(stage, row)))
            groups.add(group)
    for name in NAMED_CASES:
        for row in best_per_group(r for r in rows if r["fusion"] == name):
            group = (row["fusion"], row["pfam_family"])
            if group not in groups:
                chosen.append(dict(row, shortlist_reason=note("named_case", row)))
                groups.add(group)
    for rank, row in enumerate(chosen, start=1):
        row["shortlist_rank"] = rank
        row["schema_E_max"] = None if math.isinf(schema_e_max(row)) else schema_e_max(row)
    return chosen


def tier_counts(rows: list[dict[str, Any]], records: list[ChimeraRecord]) -> dict[str, Any]:
    """Tier counts per cohort and junction class, by record and by screen row."""
    by_id = {record.record_id: record for record in records}
    record_counts: dict[str, Counter[str]] = defaultdict(Counter)
    row_counts: dict[str, Counter[str]] = defaultdict(Counter)
    overall: Counter[str] = Counter()
    for row in rows:
        record = by_id[row["record_id"]]
        overall[row["tier"]] += 1
        for cohort, n in record.rows_per_cohort.items():
            group = f"{cohort}|{record.junction_class}"
            record_counts[group][row["tier"]] += 1
            row_counts[group][row["tier"]] += n
    return {
        "overall_records": dict(overall),
        "records": {k: dict(v) for k, v in sorted(record_counts.items())},
        "rows": {k: dict(v) for k, v in sorted(row_counts.items())},
    }


def _cell(value: Any) -> Any:
    if value is None:
        return ""
    if isinstance(value, float) and math.isinf(value):
        return ""
    return value


def write_tsv(rows: list[dict[str, Any]], path: Path) -> None:
    columns: list[str] = []
    for row in rows:
        for column in row:
            if column not in columns:
                columns.append(column)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: _cell(row.get(k, NOT_ASSESSED)) for k in columns})


def write_records(records: list[ChimeraRecord], path: Path, provenance: Mapping[str, Any]) -> None:
    payload = {
        "schema_version": SCHEMA_VERSION,
        "flank": FLANK,
        **provenance,
        "records": [record_to_dict(record) for record in records],
    }
    with gzip.open(path, "wt") as handle:
        json.dump(payload, handle, separators=(",", ":"))
