"""Validate the domain-retention track (``domain_retention_outliers.svg``)
against real BRAF and RET artifacts from the committed cohort scan.

No network access, no synthetic fixtures: reads the exact
``results.json``/``visualizations/domain_retention_outliers.svg`` already
committed to the repo and checks that the committed SVG carries the real
domain name/boundaries and real exon numbers from that gene's own
canonical transcript -- not merely that the file exists or draws an
unlabeled highlight block.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

REPO_ROOT = Path(__file__).parent.parent
GENE_REPORTS = (
    REPO_ROOT
    / "runs"
    / "cohort-scan_msk_impact_50k_2026_20260904T144201Z"
    / "cohort_scan"
    / "gene_reports"
)
RUN_DIRS = {
    "BRAF": GENE_REPORTS / "braf",
    "RET": GENE_REPORTS / "ret",
}


def _payload(run_dir: Path) -> dict:
    return json.loads((run_dir / "results.json").read_text())


def _svg(run_dir: Path) -> str:
    return (run_dir / "visualizations" / "domain_retention_outliers.svg").read_text()


def test_real_runs_have_gene_track_with_exon_boundaries_and_domains():
    for gene, run_dir in RUN_DIRS.items():
        payload = _payload(run_dir)
        assert payload["gene_track"] is not None, gene
        assert payload["gene_track"]["exon_boundaries_aa"], gene
        assert payload["gene_track"]["domains"], gene


def test_domain_retention_track_shows_real_domain_name_and_boundaries():
    for gene, run_dir in RUN_DIRS.items():
        payload = _payload(run_dir)
        accession = payload["summary"]["domain_accession"]
        start = payload["summary"]["domain_start_aa"]
        end = payload["summary"]["domain_end_aa"]
        name = next(
            d["name"] for d in payload["gene_track"]["domains"] if d["accession"] == accession
        )
        svg = _svg(run_dir)
        assert f"{name} ({start}-{end})" in svg, gene


def test_domain_retention_track_shows_real_exon_numbers():
    for gene, run_dir in RUN_DIRS.items():
        payload = _payload(run_dir)
        svg = _svg(run_dir)
        exon_ranks = {b["exon_rank"] for b in payload["gene_track"]["exon_boundaries_aa"]}
        assert exon_ranks, gene
        for rank in exon_ranks:
            assert f">E{rank}<" in svg, (gene, rank)


def test_domain_retention_track_x_axis_spans_the_real_protein_length():
    for gene, run_dir in RUN_DIRS.items():
        payload = _payload(run_dir)
        svg = _svg(run_dir)
        protein_length = payload["gene_track"]["protein_length"]
        assert f">{protein_length} aa<" in svg, gene
        assert ">1<" in svg, gene
        assert ">5' / N-terminus<" in svg, gene
        assert ">3' / C-terminus<" in svg, gene
        # At least one interior hundred-aa tick between 0 and the protein
        # length actually appears (not just the two endpoints).
        interior_hundreds = [str(hundred) for hundred in range(100, protein_length, 100)]
        assert any(f">{value}<" in svg for value in interior_hundreds), gene


def test_domain_retention_track_svg_is_well_formed():
    for gene, run_dir in RUN_DIRS.items():
        svg = _svg(run_dir)
        assert svg.strip().startswith("<svg")
        assert svg.strip().endswith("</svg>")
        match = re.search(r'width="(\d+)" height="(\d+)"', svg)
        assert match is not None, gene
        assert int(match.group(1)) > 0 and int(match.group(2)) > 0
