"""Flag domain-screen rows from several cohorts and write the pilot outputs.

Usage::

    python -I scripts/domain_screen_pilot_flags.py \\
        --screen COHORT=BUILD=DIR ...           # each DIR holds domain_screen.tsv
        --tcga-gtex NORMALIZED_TCGA_GTEX.tsv \\
        --repeats BUILD=SUPERDUPS.txt.gz,RMSK.txt.gz ...
        --reference BUILD=GENOME.2bit ...        # local reference windows (else Ensembl)
        --cache-dir FRESH_FLAG_CACHE --output-dir runs/<pilot>

Inputs are read as data only. For every cohort it writes
``<cohort>.domain_screen.tsv`` (every shared-family row, with the flag
columns), ``combined_candidates.tsv`` (all cohorts, rows that reached the
junction-translation step) and ``flag_summary.json``.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any

from cfh.cohort.artefact_flags import (
    FLAG_COLUMNS,
    EnsemblSequenceSource,
    RepeatTracks,
    TwoBitReference,
    annotate_rows,
)
from cfh.cohort.domain_screen import SCREEN_COLUMNS
from cfh.ingestion.open_fusion_sv_adapter import OPEN_PASSTHROUGH_COLUMNS

csv.field_size_limit(sys.maxsize)

NOT_TRANSLATED = {"no_shared_family", "no_transcript_model"}


def _read(path: Path, cohort: str) -> list[dict[str, Any]]:
    with path.open(newline="") as handle:
        return [{"cohort": cohort, **row} for row in csv.DictReader(handle, delimiter="\t")]


def _write(path: Path, rows: list[dict[str, Any]], columns: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--screen", action="append", required=True)
    parser.add_argument("--tcga-gtex", type=Path, required=True)
    parser.add_argument("--repeats", action="append", default=[])
    parser.add_argument("--reference", action="append", default=[])
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)

    cohorts: list[tuple[str, str, Path]] = []
    for spec in args.screen:
        cohort, build, directory = spec.split("=", 2)
        cohorts.append((cohort, build, Path(directory)))
    builds = sorted({build for _, build, _ in cohorts})

    repeat_tracks: dict[str, RepeatTracks | None] = {build: None for build in builds}
    unavailable = {
        build: f"no UCSC segdup/RepeatMasker tables given for {build}" for build in builds
    }
    for spec in args.repeats:
        build, files = spec.split("=", 1)
        superdups, rmsk = (Path(part) for part in files.split(","))
        repeat_tracks[build] = RepeatTracks(build, superdups, rmsk)

    references = {
        build: TwoBitReference(path)
        for build, path in (spec.split("=", 1) for spec in args.reference)
    }
    sequences = {
        build: EnsemblSequenceSource(
            build, args.cache_dir / f"flags_{build}", reference=references.get(build)
        )
        for build in builds
    }
    liftover = sequences["GRCh37"].liftover_to_grch38 if "GRCh37" in sequences else None

    rows: list[dict[str, Any]] = []
    summaries: dict[str, Any] = {}
    for cohort, _build, directory in cohorts:
        rows.extend(_read(directory / "domain_screen.tsv", cohort))
        summaries[cohort] = json.loads((directory / "domain_screen_summary.json").read_text())

    annotated = annotate_rows(
        rows,
        sequences=sequences,
        repeat_tracks=repeat_tracks,
        repeat_unavailable=unavailable,
        tcga_gtex_table=args.tcga_gtex,
        liftover=liftover,
        progress=lambda line: print(line, file=sys.stderr, flush=True),
        prefetch_workers=args.workers,
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    columns = ["cohort", *SCREEN_COLUMNS, *OPEN_PASSTHROUGH_COLUMNS, *FLAG_COLUMNS]
    flag_counts: dict[str, Any] = {}
    for cohort, _build, _directory in cohorts:
        own = [row for row in annotated if row["cohort"] == cohort]
        _write(args.output_dir / f"{cohort}.domain_screen.tsv", own, columns)
        flag_counts[cohort] = {
            flag: dict(Counter(str(row[flag]) for row in own))
            for flag in (
                "paralog_sequence_flag",
                "genomic_similarity_flag",
                "breakpoint_repeat_flag",
                "sample_burden_flag",
                "read_evidence_flag",
                "normal_recurrence_flag",
            )
        }
    translated = [row for row in annotated if row["reconstitution_class"] not in NOT_TRANSLATED]
    _write(args.output_dir / "combined_candidates.tsv", translated, columns)
    (args.output_dir / "flag_summary.json").write_text(
        json.dumps(
            {"screen_summaries": summaries, "flag_counts": flag_counts}, indent=2, sort_keys=True
        )
        + "\n"
    )
    print(f"flagged {len(annotated)} rows; {len(translated)} reached junction translation")


if __name__ == "__main__":
    main()
