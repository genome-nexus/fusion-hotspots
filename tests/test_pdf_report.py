"""Tests that a real, valid, non-trivial PDF is produced from a committed run.

This exercises the full PDF-rendering path (``cfh.reporting.pdf``) against
the already-committed real, latest ``braf_msk-impact-2017_*`` benchmark
artifact under ``runs/`` -- no network access, no new fixtures -- and
asserts the extracted text contains the actual numbers from that run's
results.json, not just that "a PDF was produced".
"""

from __future__ import annotations

import csv
import json

from pypdf import PdfReader

from cfh.reporting.pdf import _format_cell_value, render_pdf_report
from cfh.reporting.text import format_stat
from conftest import latest_run_dir

BRAF_RUN_DIR = latest_run_dir("braf_msk-impact-2017")


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
    # The report prefers a configured domain's human-readable name over its
    # bare accession when one is available (see ``cfh.reporting.domain_names``);
    # this run's ``summary.key_domains`` supplies "Protein kinase domain" for
    # PF07714, so that's what actually appears in the rendered text.
    domain_name = next(
        (
            domain["name"]
            for domain in summary.get("key_domains") or []
            if domain.get("accession") == summary["domain_accession"]
        ),
        summary["domain_accession"],
    )
    assert domain_name in text
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


def test_format_cell_value_summarizes_long_lists_and_passes_through_scalars():
    assert _format_cell_value("plain") == "plain"
    assert _format_cell_value(42) == 42
    assert _format_cell_value(("a", "b")) == "a, b"
    long_value = tuple(f"event_{i}" for i in range(200))
    formatted = _format_cell_value(long_value)
    assert isinstance(formatted, str)
    assert "200 total" in formatted
    assert formatted.count(",") < 200


def test_pdf_report_does_not_overflow_on_a_table_with_a_long_list_column(tmp_path):
    """Regression test: a real algorithm table (e.g. window_detection's
    ``window_scan``/``top_windows``) can carry an ``event_ids_inside``
    column holding dozens to hundreds of recurrent event ids per row. Before
    ``_format_cell_value`` existed, rendering that column verbatim produced
    one gigantic multi-thousand-point-tall table cell that overflowed the
    landscape page layout with a ``LayoutError``, aborting the whole report.
    """
    payload = json.loads((BRAF_RUN_DIR / "results.json").read_text())
    payload = dict(payload)
    payload["algorithm_results"] = [
        *payload["algorithm_results"],
        {
            "Algorithm": "window_detection",
            "Tables": {
                "window_scan": [
                    {
                        "start_aa": 100,
                        "end_aa": 200,
                        "event_ids_inside": [f"EVT-{i:04d}" for i in range(150)],
                    }
                ]
            },
        },
    ]

    output_path = tmp_path / "report.pdf"
    render_pdf_report(
        payload,
        output_path,
        results_tsv_path=BRAF_RUN_DIR / "results.tsv",
        visualizations_dir=BRAF_RUN_DIR / "visualizations",
    )

    assert output_path.exists()
    reader = PdfReader(str(output_path))
    text = "".join(page.extract_text() or "" for page in reader.pages)
    assert "150 total" in text


def test_pdf_report_paginates_many_tall_narrow_results_rows(tmp_path):
    """Long 24-column result tables do not enter ReportLab's split failure.

    A tall source annotation is still individually shorter than a landscape
    page, but its repeated header made ReportLab fail when it appeared late in
    a large table.  This is deliberately synthetic so it protects every gene,
    rather than just the RET/FGFR2 artifacts that exposed the issue.

    The annotation text is a single unbroken 300-char run (no spaces) rather
    than repeated short words: ReportLab can only line-wrap at word
    boundaries, so a long *space-separated* value wraps into many short lines
    and never actually gets tall enough to reproduce the original
    ``LayoutError`` -- confirmed by running this exact scenario against
    ``pdf.py`` as of commit 4f9db8c (the merge-base before this pagination
    fix existed): with space-separated text it passes even on that
    unpatched code, i.e. it protects nothing. The unbroken variant below
    does reproduce the crash on 4f9db8c and only passes with the fix.
    """
    payload = json.loads((BRAF_RUN_DIR / "results.json").read_text())
    payload = {**payload, "gene_symbol": "PAGINATION_TEST", "algorithm_results": []}
    header = [f"column_{index}" for index in range(24)]
    header[22] = "source_annotation_text"
    annotation = "syntheticannotation" + "x" * (300 - len("syntheticannotation"))
    rows = [header]
    for index in range(180):
        row = [f"value-{index}" for _ in header]
        row[0] = f"event-{index:03d}"
        row[22] = annotation if index in {31, 87, 151} else "short annotation"
        rows.append(row)

    tsv_path = tmp_path / "results.tsv"
    with tsv_path.open("w", newline="") as handle:
        csv.writer(handle, delimiter="\t").writerows(rows)

    output_path = tmp_path / "report.pdf"
    render_pdf_report(payload, output_path, results_tsv_path=tsv_path, visualizations_dir=tmp_path)

    reader = PdfReader(str(output_path))
    text = "".join(page.extract_text() or "" for page in reader.pages)
    assert len(reader.pages) > 10
    assert "PAGINATION_TEST fusion-hotspot benchmark report" in text
    assert "event-179" in "".join(text.split())
    assert "syntheticannotation" in "".join(text.split())


def test_pdf_report_paginates_table_with_long_unbroken_header_text(tmp_path):
    """A long, unbroken header value must not overflow the page either.

    ``_generic_table_flowable`` bounds body-cell height, but before this
    fix the header row (row 0) only ever went through the old
    character-count-only ``_truncate_cell_text`` -- unbounded in rendered
    height. ``_page_sized_table_flowables`` also never validated a bare
    ``header + first row`` candidate before accepting it (its overflow
    guard was skipped whenever ``current_rows`` was still empty), so an
    over-tall header combined with even one sufficiently tall body row was
    silently handed to ReportLab as a raw ``Table`` and crashed with a
    ``LayoutError`` -- reproduced live against commit 81b7e32 (this PR's
    prior state) using the exact scenario below.
    """
    payload = json.loads((BRAF_RUN_DIR / "results.json").read_text())
    payload = {**payload, "gene_symbol": "PAGINATION_TEST", "algorithm_results": []}
    header = [f"column_{index}" for index in range(24)]
    long_header = "fusion_annotation_crosscheck_protein_position"
    header[22] = long_header + "z" * (300 - len(long_header))
    annotation = "tallfirstrowannotation" + "w" * (300 - len("tallfirstrowannotation"))
    rows = [header]
    for index in range(180):
        row = [f"value-{index}" for _ in header]
        row[0] = f"event-{index:03d}"
        row[23] = annotation if index == 0 else "short"
        rows.append(row)

    tsv_path = tmp_path / "results.tsv"
    with tsv_path.open("w", newline="") as handle:
        csv.writer(handle, delimiter="\t").writerows(rows)

    output_path = tmp_path / "report.pdf"
    render_pdf_report(payload, output_path, results_tsv_path=tsv_path, visualizations_dir=tmp_path)

    reader = PdfReader(str(output_path))
    text = "".join(page.extract_text() or "" for page in reader.pages)
    assert len(reader.pages) > 5
    assert long_header in "".join(text.split())
    assert "event-179" in "".join(text.split())
