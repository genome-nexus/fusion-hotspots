#!/usr/bin/env python3
"""Assemble a static GitHub Pages site from committed ``cfh`` run artifacts.

For every directory under ``runs/`` this (re)generates the backend-less HTML
viewer bundle via :func:`cfh.reporting.html_viewer.build_run_viewer` -- the
same offline, no-network function exercised by
``tests/test_html_viewer.py`` against the real committed BRAF and
cohort-scan runs -- then copies each bundle plus a searchable gene landing
page (:mod:`cfh.reporting.site_index`) into a single output directory
ready to hand to ``actions/upload-pages-artifact``.

This script only reads already-committed run artifacts (``results.json``,
``cohort_scan/summary.json``, ``*.svg``) and writes the generated
``viewer/`` bundle into each run directory (mirroring what
``cfh analyze``/``cfh cohort-scan --html`` does locally). It never
re-runs ingestion, mapping, or algorithm code, and it makes no network
calls.

Usage:
    python scripts/build_pages_site.py --runs-dir runs --output-dir _site
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

from cfh.reporting.html_viewer import build_run_viewer
from cfh.reporting.site_index import collect_site_index, render_site_index


def build_site(runs_dir: Path, output_dir: Path) -> list[str]:
    """Generate a viewer bundle for every eligible run under ``runs_dir``
    and copy it into ``output_dir``.

    Returns the sorted list of run directory names that got a published
    viewer (a run directory with neither a single-gene ``results.json`` nor
    a cohort-scan ``cohort_scan/summary.json`` is silently skipped, matching
    :func:`build_run_viewer`).
    """
    if output_dir.exists():
        shutil.rmtree(output_dir)
    (output_dir / "runs").mkdir(parents=True)

    published: list[str] = []
    for run_dir in sorted(p for p in runs_dir.iterdir() if p.is_dir()):
        if build_run_viewer(run_dir) is None:
            continue
        shutil.copytree(run_dir / "viewer", output_dir / "runs" / run_dir.name)
        published.append(run_dir.name)

    index = collect_site_index(runs_dir, published)
    (output_dir / "index.html").write_text(render_site_index(index))
    return published


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runs-dir", type=Path, default=Path("runs"))
    parser.add_argument("--output-dir", type=Path, default=Path("_site"))
    args = parser.parse_args(argv)

    runs_dir: Path = args.runs_dir
    if not runs_dir.is_dir():
        print(f"error: {runs_dir} is not a directory", file=sys.stderr)
        return 1

    published = build_site(runs_dir, args.output_dir)
    if not published:
        print(f"error: no run under {runs_dir} produced a viewer bundle", file=sys.stderr)
        return 1

    print(f"Published {len(published)} run viewer(s) to {args.output_dir}:")
    for name in published:
        print(f"  - {name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
