"""Assemble the self-contained, human-reviewer-facing ``report.pdf``.

Consumes exactly the artifacts already written by
``cfh.real_benchmark.write_outputs`` for one run directory (``results.json``,
``results.tsv``, ``visualizations/*.svg``) plus the deterministic templated
text from :mod:`cfh.reporting.text`. Layout only -- no numbers are computed
here; every figure comes from ``text.py`` or is read verbatim from the run's
own tables.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.pagesizes import LETTER, landscape
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.platypus import (
    NextPageTemplate,
    PageBreak,
    PageTemplate,
    Paragraph,
    Spacer,
    Table,
    TableStyle,
)
from reportlab.platypus.doctemplate import BaseDocTemplate
from reportlab.platypus.frames import Frame
from svglib.svglib import svg2rlg

from cfh.reporting.text import render_abstract, render_results_summary

_PORTRAIT_TEMPLATE = "portrait"
_LANDSCAPE_TEMPLATE = "landscape"

_MAX_TABLE_ROWS = 500
"""Hard cap on rendered TSV rows so a pathologically large run stays a
readable, boundedly-sized PDF instead of a runaway multi-thousand-page
document; a note is appended when rows are truncated."""

_MAX_LIST_CELL_ITEMS = 5
"""Cap on how many list/tuple elements are spelled out inline in one
rendered table cell (e.g. an algorithm's ``event_ids_inside`` column). A
cell holding a long list (dozens to hundreds of recurrent-event ids, as a
window-scan candidate can) would otherwise render as one enormous
multi-thousand-point-tall cell that overflows the page layout; the full
list is still available verbatim in ``results.json``/``results.tsv``, this
is a display-only summary."""


def _format_cell_value(value: object) -> object:
    """Summarize a list/tuple cell value for display; pass everything else through."""
    if isinstance(value, (list, tuple)):
        if len(value) <= _MAX_LIST_CELL_ITEMS:
            return ", ".join(str(item) for item in value)
        shown = ", ".join(str(item) for item in value[:_MAX_LIST_CELL_ITEMS])
        return f"{shown}, ... ({len(value)} total)"
    return value


def _styles() -> dict[str, ParagraphStyle]:
    base = getSampleStyleSheet()
    return {
        "Title": base["Title"],
        "Heading1": base["Heading1"],
        "Heading2": base["Heading2"],
        "Body": base["BodyText"],
        "Cell": ParagraphStyle("Cell", parent=base["BodyText"], fontSize=6, leading=7.5),
        "CellHeader": ParagraphStyle(
            "CellHeader", parent=base["BodyText"], fontSize=6.5, leading=8, textColor=colors.white
        ),
        "Caption": ParagraphStyle("Caption", parent=base["Italic"], fontSize=9),
    }


def _load_svg_drawing(path: Path, max_width: float):
    drawing = svg2rlg(str(path))
    if drawing is None or not drawing.width:
        return None
    scale = min(1.0, max_width / drawing.width)
    if scale < 1.0:
        drawing.width *= scale
        drawing.height *= scale
        drawing.scale(scale, scale)
    return drawing


def _generic_table_flowable(rows: list[list], styles: dict, header: bool = True) -> Table:
    cell_style = styles["Cell"]
    header_style = styles["CellHeader"]
    formatted = []
    for row_index, row in enumerate(rows):
        style = header_style if header and row_index == 0 else cell_style
        formatted.append([Paragraph(str(value), style) for value in row])
    table = Table(formatted, repeatRows=1 if header else 0)
    table_style = [
        ("GRID", (0, 0), (-1, -1), 0.4, colors.grey),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("ROWBACKGROUNDS", (0, 1 if header else 0), (-1, -1), [colors.white, colors.whitesmoke]),
    ]
    if header:
        table_style.append(("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#2878b5")))
    table.setStyle(TableStyle(table_style))
    return table


def _results_tsv_table(tsv_path: Path, styles: dict) -> list:
    if not tsv_path.exists():
        return []
    with tsv_path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        rows = list(reader)
    if not rows:
        return []
    flowables: list = [Paragraph("Per-event results (results.tsv)", styles["Heading2"])]
    body_rows = rows[1:]
    truncated = len(body_rows) > _MAX_TABLE_ROWS
    display_rows = [rows[0]] + body_rows[:_MAX_TABLE_ROWS]
    flowables.append(_generic_table_flowable(display_rows, styles))
    if truncated:
        flowables.append(
            Paragraph(
                f"Showing the first {_MAX_TABLE_ROWS} of {len(body_rows)} rows; "
                "see results.tsv for the complete table.",
                styles["Caption"],
            )
        )
    return flowables


def _algorithm_tables(payload: dict, styles: dict) -> list:
    """Render each algorithm's own ``Tables`` entries (contingency tables,
    partner-gene counts, cutpoint scan, and -- forward-compatibly -- a
    ``composite_score`` ranked table if one is ever present) as real tables.
    """
    flowables: list = []
    for result in payload.get("algorithm_results") or []:
        tables = result.get("Tables") or {}
        rendered_any = False
        for table_name, value in tables.items():
            if isinstance(value, dict) and value.get("omitted_from_artifact"):
                continue
            if isinstance(value, list) and value and isinstance(value[0], dict):
                columns = list(value[0].keys())
                rows = [columns] + [
                    [_format_cell_value(row.get(col)) for col in columns] for row in value
                ]
            elif isinstance(value, list) and value and isinstance(value[0], list):
                rows = value
            else:
                continue
            if not rendered_any:
                flowables.append(
                    Paragraph(
                        f"{result.get('Algorithm')} tables",
                        styles["Heading2"],
                    )
                )
                rendered_any = True
            flowables.append(Paragraph(table_name, styles["Caption"]))
            display_rows = rows[: _MAX_TABLE_ROWS + 1]
            flowables.append(
                _generic_table_flowable(display_rows, styles, header=isinstance(value[0], dict))
            )
            if len(rows) > _MAX_TABLE_ROWS + 1:
                flowables.append(
                    Paragraph(
                        f"Showing the first {_MAX_TABLE_ROWS} of {len(rows) - 1} rows.",
                        styles["Caption"],
                    )
                )
            flowables.append(Spacer(1, 0.15 * inch))
    return flowables


def _figures(visualization_dir: Path, styles: dict, max_width: float) -> list:
    flowables: list = []
    if not visualization_dir.exists():
        return flowables
    svg_paths = sorted(visualization_dir.glob("*.svg"))
    if not svg_paths:
        return flowables
    flowables.append(Paragraph("Figures", styles["Heading1"]))
    for svg_path in svg_paths:
        drawing = _load_svg_drawing(svg_path, max_width)
        if drawing is None:
            continue
        flowables.append(Paragraph(svg_path.stem.replace("_", " "), styles["Heading2"]))
        flowables.append(drawing)
        flowables.append(Spacer(1, 0.2 * inch))
    return flowables


def render_pdf_report(
    payload: dict,
    output_path: str | Path,
    *,
    results_tsv_path: str | Path | None = None,
    visualizations_dir: str | Path | None = None,
) -> Path:
    """Render one run's PDF report to ``output_path`` and return that path.

    ``payload`` is the same dict shape written to (and read back from)
    ``results.json``. ``results_tsv_path``/``visualizations_dir`` default to
    ``results.tsv``/``visualizations`` next to ``output_path`` -- the
    standard ``runs/<run_id>/`` layout -- but can be overridden (e.g. in
    tests using synthetic fixtures that live elsewhere).
    """
    output_path = Path(output_path)
    run_dir = output_path.parent
    tsv_path = Path(results_tsv_path) if results_tsv_path else run_dir / "results.tsv"
    viz_dir = Path(visualizations_dir) if visualizations_dir else run_dir / "visualizations"

    styles = _styles()
    portrait_width, portrait_height = LETTER
    landscape_size = landscape(LETTER)

    doc = BaseDocTemplate(
        str(output_path),
        pagesize=LETTER,
        leftMargin=0.75 * inch,
        rightMargin=0.75 * inch,
        topMargin=0.75 * inch,
        bottomMargin=0.75 * inch,
    )
    portrait_frame = Frame(
        doc.leftMargin,
        doc.bottomMargin,
        portrait_width - doc.leftMargin - doc.rightMargin,
        portrait_height - doc.topMargin - doc.bottomMargin,
        id="portrait",
    )
    landscape_frame = Frame(
        doc.leftMargin,
        doc.bottomMargin,
        landscape_size[0] - doc.leftMargin - doc.rightMargin,
        landscape_size[1] - doc.topMargin - doc.bottomMargin,
        id="landscape",
    )
    doc.addPageTemplates(
        [
            PageTemplate(id=_PORTRAIT_TEMPLATE, frames=[portrait_frame], pagesize=LETTER),
            PageTemplate(id=_LANDSCAPE_TEMPLATE, frames=[landscape_frame], pagesize=landscape_size),
        ]
    )

    gene = payload.get("gene_symbol") or "Unknown gene"
    study = payload.get("study_id") or "unknown study"

    story: list = [
        Paragraph(f"{gene} fusion-hotspot benchmark report", styles["Title"]),
        Paragraph(f"Study: {study}", styles["Body"]),
        Spacer(1, 0.2 * inch),
        Paragraph("Abstract", styles["Heading1"]),
        Paragraph(render_abstract(payload), styles["Body"]),
        Spacer(1, 0.25 * inch),
        Paragraph("Results summary", styles["Heading1"]),
    ]
    for section in render_results_summary(payload):
        story.append(Paragraph(section["heading"], styles["Heading2"]))
        story.append(Paragraph(section["paragraph"], styles["Body"]))
        story.append(Spacer(1, 0.1 * inch))

    story.append(NextPageTemplate(_LANDSCAPE_TEMPLATE))
    story.append(PageBreak())
    story.append(Paragraph("Tables", styles["Heading1"]))
    story.extend(_algorithm_tables(payload, styles))
    story.extend(_results_tsv_table(tsv_path, styles))

    story.append(NextPageTemplate(_PORTRAIT_TEMPLATE))
    story.append(PageBreak())
    story.extend(_figures(viz_dir, styles, portrait_width - doc.leftMargin - doc.rightMargin))

    doc.build(story)
    return output_path


def render_cohort_summary_pdf(
    output_path: str | Path,
    *,
    title: str,
    subtitle: str,
    notes: list[str],
    rows: list[list],
    extra_tables: list[dict] | None = None,
    figures_dir: str | Path | None = None,
) -> Path:
    """Render a simple, landscape, one-table PDF summarizing every gene in a
    cohort scan, reusing the same table-flowable styling as the per-gene
    ``report.pdf`` (see :func:`_generic_table_flowable`). ``rows`` is a
    header row followed by one row per scanned gene, already
    string-formatted by the caller (no numbers are computed here).

    ``extra_tables``, if given, is a list of ``{"heading", "note", "rows"}``
    dicts rendered as additional sections after the main table -- e.g. the
    "honorable mentions" highly ranked non-FDR-significant tier -- each with
    its own heading, an italic caption, and a header-row-first table exactly
    like the main one.

    ``figures_dir``, if given, is rendered exactly like the per-gene
    report's own figures section (see :func:`_figures`) -- e.g. the
    genome-wide Manhattan/volcano summary SVG written alongside
    ``summary.pdf`` in the same ``cohort_scan/`` directory.
    """
    output_path = Path(output_path)
    styles = _styles()
    landscape_size = landscape(LETTER)

    doc = BaseDocTemplate(
        str(output_path),
        pagesize=landscape_size,
        leftMargin=0.5 * inch,
        rightMargin=0.5 * inch,
        topMargin=0.6 * inch,
        bottomMargin=0.6 * inch,
    )
    frame = Frame(
        doc.leftMargin,
        doc.bottomMargin,
        landscape_size[0] - doc.leftMargin - doc.rightMargin,
        landscape_size[1] - doc.topMargin - doc.bottomMargin,
        id="landscape",
    )
    doc.addPageTemplates(
        [PageTemplate(id=_LANDSCAPE_TEMPLATE, frames=[frame], pagesize=landscape_size)]
    )

    story: list = [
        Paragraph(title, styles["Title"]),
        Paragraph(subtitle, styles["Body"]),
        Spacer(1, 0.15 * inch),
    ]
    for note in notes:
        story.append(Paragraph(note, styles["Caption"]))
    story.append(Spacer(1, 0.2 * inch))

    truncated = len(rows) - 1 > _MAX_TABLE_ROWS
    display_rows = [rows[0]] + rows[1 : _MAX_TABLE_ROWS + 1] if rows else rows
    if display_rows:
        story.append(_generic_table_flowable(display_rows, styles))
    if truncated:
        story.append(
            Paragraph(
                f"Showing the first {_MAX_TABLE_ROWS} of {len(rows) - 1} scanned genes.",
                styles["Caption"],
            )
        )

    for extra_table in extra_tables or []:
        story.append(Spacer(1, 0.25 * inch))
        story.append(Paragraph(extra_table["heading"], styles["Heading1"]))
        if extra_table.get("note"):
            story.append(Paragraph(extra_table["note"], styles["Caption"]))
        story.append(Spacer(1, 0.1 * inch))
        extra_rows = extra_table["rows"]
        if extra_rows:
            story.append(_generic_table_flowable(extra_rows, styles))

    if figures_dir is not None:
        story.append(PageBreak())
        story.extend(
            _figures(
                Path(figures_dir),
                styles,
                landscape_size[0] - doc.leftMargin - doc.rightMargin,
            )
        )

    doc.build(story)
    return output_path


def render_manuscript_pdf(
    output_path: str | Path,
    *,
    title: str,
    abstract: str,
    methods: str,
    manhattan_svg_path: str | Path | None,
    manhattan_caption: str,
    results_table_rows: list[list],
    gene_highlights: list[dict],
    discussion_bullets: list[str],
    appendix_rows: list[list],
) -> Path:
    """Render the cross-gene manuscript-style synthesis report
    (``paper.pdf``) to ``output_path`` and return that path.

    Layout only -- every string/number here is already computed by
    :mod:`cfh.reporting.manuscript_text` and :mod:`cfh.cohort.outputs`, and
    every figure is a pre-existing SVG on disk (the cohort scan's own
    ``manhattan.svg`` and, per highlighted gene, that gene's own already
    -generated ``report.pdf`` figure) -- nothing is regenerated here.

    ``gene_highlights`` is a list of ``{"heading", "paragraph",
    "figure_path", "figure_caption", "report_note"}`` dicts, one per
    highlighted gene (FDR-significant, honorable-mention, and hand-curated
    genes); ``figure_path``/``figure_caption``/``report_note`` may be
    ``None``. ``results_table_rows`` and ``appendix_rows`` are each a header
    row followed by data rows, already string-formatted by the caller.
    """
    output_path = Path(output_path)
    styles = _styles()
    portrait_width, portrait_height = LETTER
    landscape_size = landscape(LETTER)

    doc = BaseDocTemplate(
        str(output_path),
        pagesize=LETTER,
        leftMargin=0.75 * inch,
        rightMargin=0.75 * inch,
        topMargin=0.75 * inch,
        bottomMargin=0.75 * inch,
    )
    portrait_frame = Frame(
        doc.leftMargin,
        doc.bottomMargin,
        portrait_width - doc.leftMargin - doc.rightMargin,
        portrait_height - doc.topMargin - doc.bottomMargin,
        id="portrait",
    )
    landscape_frame = Frame(
        doc.leftMargin,
        doc.bottomMargin,
        landscape_size[0] - doc.leftMargin - doc.rightMargin,
        landscape_size[1] - doc.topMargin - doc.bottomMargin,
        id="landscape",
    )
    doc.addPageTemplates(
        [
            PageTemplate(id=_PORTRAIT_TEMPLATE, frames=[portrait_frame], pagesize=LETTER),
            PageTemplate(id=_LANDSCAPE_TEMPLATE, frames=[landscape_frame], pagesize=landscape_size),
        ]
    )

    content_width = portrait_width - doc.leftMargin - doc.rightMargin

    story: list = [
        Paragraph(title, styles["Title"]),
        Spacer(1, 0.2 * inch),
        Paragraph("Abstract", styles["Heading1"]),
        Paragraph(abstract, styles["Body"]),
        Spacer(1, 0.2 * inch),
        Paragraph("Methods", styles["Heading1"]),
        Paragraph(methods, styles["Body"]),
        Spacer(1, 0.25 * inch),
        Paragraph("Results", styles["Heading1"]),
        Paragraph("Genome-wide summary", styles["Heading2"]),
        Paragraph(manhattan_caption, styles["Caption"]),
    ]
    if manhattan_svg_path is not None and Path(manhattan_svg_path).exists():
        drawing = _load_svg_drawing(Path(manhattan_svg_path), content_width)
        if drawing is not None:
            story.append(drawing)
    story.append(Spacer(1, 0.2 * inch))

    if results_table_rows:
        story.append(NextPageTemplate(_LANDSCAPE_TEMPLATE))
        story.append(PageBreak())
        story.append(Paragraph("FDR-significant and honorable-mention genes", styles["Heading2"]))
        story.append(_generic_table_flowable(results_table_rows, styles))
        story.append(NextPageTemplate(_PORTRAIT_TEMPLATE))
        story.append(PageBreak())

    story.append(Paragraph("Gene highlights", styles["Heading2"]))
    for highlight in gene_highlights:
        story.append(Paragraph(highlight["heading"], styles["Heading2"]))
        story.append(Paragraph(highlight["paragraph"], styles["Body"]))
        figure_path = highlight.get("figure_path")
        if figure_path is not None and Path(figure_path).exists():
            drawing = _load_svg_drawing(Path(figure_path), content_width)
            if drawing is not None:
                if highlight.get("figure_caption"):
                    story.append(Paragraph(highlight["figure_caption"], styles["Caption"]))
                story.append(drawing)
        if highlight.get("report_note"):
            story.append(Paragraph(highlight["report_note"], styles["Caption"]))
        story.append(Spacer(1, 0.15 * inch))

    story.append(Paragraph("Discussion", styles["Heading1"]))
    for bullet in discussion_bullets:
        story.append(Paragraph(f"&bull; {bullet}", styles["Body"]))
        story.append(Spacer(1, 0.05 * inch))

    if appendix_rows:
        story.append(NextPageTemplate(_LANDSCAPE_TEMPLATE))
        story.append(PageBreak())
        story.append(Paragraph("Appendix: per-gene report index", styles["Heading1"]))
        story.append(_generic_table_flowable(appendix_rows, styles))

    doc.build(story)
    return output_path


def render_pdf_report_for_run_dir(run_dir: str | Path) -> Path:
    """Convenience wrapper: render ``report.pdf`` for an on-disk run directory
    that already has ``results.json`` (and, if present, ``results.tsv`` /
    ``visualizations/``) written by ``write_outputs``.
    """
    run_dir = Path(run_dir)
    payload = json.loads((run_dir / "results.json").read_text())
    return render_pdf_report(payload, run_dir / "report.pdf")
