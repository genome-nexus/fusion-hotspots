"""Tests that a real, valid, non-trivial PDF is produced from a committed run.

This exercises the full PDF-rendering path (``cfh.reporting.pdf``) against
the already-committed real BRAF benchmark artifact under
``runs/braf_msk-impact-2017_20260905T012645Z/`` -- no network access, no
new fixtures -- and asserts the extracted text contains the actual numbers
from that run's results.json, not just that "a PDF was produced".
"""

from __future__ import annotations

import json
from pathlib import Path

from pypdf import PdfReader

from cfh.reporting.pdf import render_pdf_report
from cfh.reporting.text import format_stat

REPO_ROOT = Path(__file__).parent.parent
BRAF_RUN_DIR = REPO_ROOT / "runs" / "braf_msk-impact-2017_20260905T012645Z"


def test_pdf_report_generated_from_real_braf_run_contains_actual_numbers(tmp_path):
    payload = json.loads((BRAF_RUN_DIR / "results.json").read_text())

    output_path = tmp_path / "report.pdf"
    result_path = render_pdf_report(
        payload,
        output_path,
        results_tsv_path=BRAF_RUN_DIR / "results.tsv",
        visualizations_dir=BRAF_RUN_DIR / "visualizations",
    )

    assert result_path == output_path
    assert output_path.exists()
    assert output_path.stat().st_size > 5_000  # not a trivially-empty PDF

    reader = PdfReader(str(output_path))
    assert len(reader.pages) > 1  # abstract/summary + tables + figures span pages

    text = "".join(page.extract_text() or "" for page in reader.pages)

    summary = payload["summary"]
    assert payload["gene_symbol"] in text
    assert payload["study_id"] in text
    assert str(summary["total_fusions"]) in text
    assert f"{summary['in_frame_percent']:.1f}%" in text
    assert f"{summary['kinase_retained_percent']:.1f}%" in text
    assert summary["domain_accession"] in text
    assert format_stat(summary["fisher_p_value"]) in text

    # A real partner-gene name from the run's frequency table, embedded via
    # a real PDF table (not linked/omitted).
    assert "SND1" in text

    # The abstract and results-summary section headings are present.
    assert "Abstract" in text
    assert "Results summary" in text
    assert "Fusion partner frequency" in text
    assert "Domain retention" in text

    # Figures section: the run's SVG visualizations were embedded, not just
    # linked -- their captions appear as rendered figure titles.
    assert "Figures" in text
    assert "domain retention outliers" in text
    assert "reference comparison" in text


def test_pdf_report_is_deterministic_across_repeated_renders(tmp_path):
    payload = json.loads((BRAF_RUN_DIR / "results.json").read_text())

    first = tmp_path / "first.pdf"
    second = tmp_path / "second.pdf"
    render_pdf_report(
        payload,
        first,
        results_tsv_path=BRAF_RUN_DIR / "results.tsv",
        visualizations_dir=BRAF_RUN_DIR / "visualizations",
    )
    render_pdf_report(
        payload,
        second,
        results_tsv_path=BRAF_RUN_DIR / "results.tsv",
        visualizations_dir=BRAF_RUN_DIR / "visualizations",
    )

    first_text = "".join(page.extract_text() or "" for page in PdfReader(str(first)).pages)
    second_text = "".join(page.extract_text() or "" for page in PdfReader(str(second)).pages)
    assert first_text == second_text
