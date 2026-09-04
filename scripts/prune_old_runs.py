#!/usr/bin/env python3
"""Keep only the most recent committed run artifact per (type, study) group.

Directories under ``runs/`` are named ``<type>_<ISO8601-timestamp>``, e.g.
``braf_msk-impact-50k-2026_20260904T172738Z`` or
``cohort-scan_msk_impact_50k_2026_20260904T144201Z``. Successive fix rounds
regenerate the same run type, and each regeneration used to get committed
alongside every earlier one -- bloating the repo with redundant historical
copies. This script groups run directories by everything before the trailing
timestamp, and removes every directory in a group except the one with the
latest timestamp.

It never touches anything that isn't a ``runs/<type>_<timestamp>/`` directory
(standalone files directly under ``runs/`` are left alone), and it never
hardcodes a gene name, study id, or run type -- the grouping is derived
purely from directory names.

Usage:
    python scripts/prune_old_runs.py            # dry run: report only
    python scripts/prune_old_runs.py --apply     # actually delete

Run this before committing a new run artifact under ``runs/`` so only the
latest run per (gene-or-scan-kind, study_id) combination is ever committed.
"""

from __future__ import annotations

import argparse
import re
import shutil
import sys
from pathlib import Path

# ``<anything>_<8 digits>T<6 digits>Z`` e.g. ``..._20260904T172738Z``.
_RUN_DIR_PATTERN = re.compile(r"^(?P<run_type>.+)_(?P<timestamp>\d{8}T\d{6}Z)$")


def group_runs_by_type(runs_dir: Path) -> dict[str, list[tuple[str, Path]]]:
    """Group timestamped run directories under ``runs_dir`` by their type.

    Returns a mapping of run type -> list of (timestamp, path) pairs, each
    list sorted oldest-first (ISO8601 timestamps sort lexicographically).
    Non-directories and directories that don't match the
    ``<type>_<timestamp>`` naming convention are ignored entirely.
    """
    groups: dict[str, list[tuple[str, Path]]] = {}
    for entry in sorted(runs_dir.iterdir()):
        if not entry.is_dir():
            continue
        match = _RUN_DIR_PATTERN.match(entry.name)
        if not match:
            continue
        groups.setdefault(match.group("run_type"), []).append((match.group("timestamp"), entry))
    for entries in groups.values():
        entries.sort(key=lambda pair: pair[0])
    return groups


def find_stale_runs(runs_dir: Path) -> list[Path]:
    """Return every run directory that is not the latest in its group."""
    stale: list[Path] = []
    for entries in group_runs_by_type(runs_dir).values():
        stale.extend(path for _timestamp, path in entries[:-1])
    return stale


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--runs-dir",
        type=Path,
        default=Path("runs"),
        help="Directory containing run artifacts (default: runs)",
    )
    parser.add_argument(
        "--apply",
        action="store_true",
        help="Actually delete the stale run directories (default: dry run)",
    )
    args = parser.parse_args(argv)

    runs_dir: Path = args.runs_dir
    if not runs_dir.is_dir():
        print(f"error: {runs_dir} is not a directory", file=sys.stderr)
        return 1

    stale = find_stale_runs(runs_dir)
    if not stale:
        print(f"Nothing to prune under {runs_dir}: every run type has a single latest copy.")
        return 0

    verb = "Removing" if args.apply else "Would remove (pass --apply to delete)"
    for path in stale:
        print(f"{verb}: {path}")
        if args.apply:
            shutil.rmtree(path)

    print(
        f"\n{len(stale)} redundant run director{'y' if len(stale) == 1 else 'ies'} "
        f"{'removed' if args.apply else 'found'}."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
