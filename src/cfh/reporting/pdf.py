"""Assemble the self-contained, human-reviewer-facing ``report.pdf`` and the
cross-gene ``paper.pdf``/``summary.pdf``.

Consumes exactly the artifacts already written by
``cfh.real_benchmark.write_outputs`` / ``cfh.cohort.outputs`` for one run
(``results.json``, ``results.tsv``, ``visualizations/*.svg``) plus the
deterministic templated text from :mod:`cfh.reporting.text` and
:mod:`cfh.reporting.manuscript_text`. Layout only -- no numbers are computed
here; every sentence comes from those two modules and every figure is read
verbatim from the run's own SVGs.

``render_pdf_report`` (the per-gene ``report.pdf``) and
``render_manuscript_pdf`` (the cohort-level manuscript ``paper.pdf``) render
a real two-column academic-paper-style LaTeX document, compiled with
Tectonic (see :mod:`cfh.reporting.latex` for why Tectonic and how escaping/
compilation work). ``render_cohort_summary_pdf`` (the plain one-table-per-
gene ``summary.pdf``) is unaffected -- it keeps using reportlab/platypus
directly, as before.
"""

from __future__ import annotations

import csv
import json
import tempfile
from pathlib import Path
from typing import Any

from reportlab.lib import colors
from reportlab.lib.pagesizes import LETTER, landscape
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.platypus import (
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

from cfh.reporting.latex import compile_latex, escape_latex, latex_long_table, svg_to_pdf
from cfh.reporting.text import render_abstract, render_results_summary

_PORTRAIT_TEMPLATE = "portrait"
_LANDSCAPE_TEMPLATE = "landscape"

_MAX_TABLE_ROWS = 500
"""Hard cap on rendered TSV rows so a pathologically large run stays a
readable, boundedly-sized PDF instead of a runaway multi-thousand-page
document; a note is appended when rows are truncated."""

_LATEX_PREAMBLE = r"""\documentclass[9pt,twocolumn]{article}
\usepackage[letterpaper,margin=0.65in]{geometry}
\usepackage[T1]{fontenc}
\usepackage{lmodern}
\usepackage{booktabs}
\usepackage{graphicx}
\usepackage{grffile}
\usepackage{pdflscape}
\usepackage{xltabular}
\usepackage{caption}
\setlength{\parindent}{0pt}
\setlength{\parskip}{0.55em}
\captionsetup{font=small,labelfont=bf}
\pagestyle{plain}
"""


class _FigureCollector:
    """Converts each SVG handed to :meth:`add` to a real standalone PDF
    (via :func:`cfh.reporting.latex.svg_to_pdf`, reusing this project's
    existing svglib/reportlab SVG-loading path) under a plain,
    LaTeX-special-character-free relative filename (``fig0.pdf``,
    ``fig1.pdf``, ...) inside ``scratch_dir``, and returns that filename for
    use in an ``\\includegraphics{...}`` call. Returns ``None`` -- adding
    nothing -- when the SVG has no renderable content, mirroring
    ``svg_to_pdf``'s own None-on-empty behavior.
    """

    def __init__(self, scratch_dir: Path) -> None:
        self._scratch_dir = scratch_dir
        self._next_index = 0
        self.figures: dict[str, Path] = {}

    def add(self, svg_path: Path) -> str | None:
        name = f"fig{self._next_index}.pdf"
        destination = self._scratch_dir / name
        if svg_to_pdf(svg_path, destination) is None:
            return None
        self._next_index += 1
        self.figures[name] = destination
        return name


def _title_and_abstract_tex(title: str, subtitle: str, abstract: str) -> str:
    """The shared ``\\twocolumn[...]`` title-page idiom: a full-width title
    and single-column abstract, followed by the two-column body -- see
    :mod:`cfh.reporting.latex`'s module docstring for why this toolchain/
    layout was verified against a real compile before being adopted here.
    """
    lines = [
        "\\title{" + escape_latex(title) + "}",
        "\\author{" + escape_latex(subtitle) + "}",
        "\\date{}",
        "\\begin{document}",
        "\\twocolumn[",
        "  \\begin{@twocolumnfalse}",
        "  \\maketitle",
        "  \\begin{abstract}",
        escape_latex(abstract),
        "  \\end{abstract}",
        "  \\vspace{1em}",
        "  \\end{@twocolumnfalse}",
        "]",
        "",
    ]
    return "\n".join(lines)


def _figure_block(name: str, caption: str) -> str:
    lines = [
        "\\begin{figure*}[htbp]",
        "\\centering",
        f"\\includegraphics[width=0.92\\textwidth]{{{name}}}",
    ]
    if caption:
        lines.append("\\caption{" + escape_latex(caption) + "}")
    lines.append("\\end{figure*}")
    return "\n".join(lines) + "\n"


def _landscape_table_section(heading: str, table_tex: str) -> str:
    return (
        "\\clearpage\n\\onecolumn\n\\begin{landscape}\n"
        + ("\\section{" + escape_latex(heading) + "}\n" if heading else "")
        + table_tex
        + "\\end{landscape}\n\\clearpage\n\\twocolumn\n"
    )


def _algorithm_tables_tex(payload: dict) -> str:
    """Render each algorithm's own ``Tables`` entries (contingency tables,
    partner-gene counts, cutpoint scan, and -- forward-compatibly -- a
    ``composite_score`` ranked table if one is ever present) as real
    booktabs/xltabular tables. Reused verbatim data shape from the previous
    reportlab implementation -- only the rendering target changed.
    """
    parts: list[str] = []
    for result in payload.get("algorithm_results") or []:
        algorithm = result.get("Algorithm")
        for table_name, value in (result.get("Tables") or {}).items():
            if isinstance(value, dict) and value.get("omitted_from_artifact"):
                continue
            if isinstance(value, list) and value and isinstance(value[0], dict):
                columns = list(value[0].keys())
                rows: list[list[Any]] = [columns] + [
                    [row.get(col) for col in columns] for row in value
                ]
            elif isinstance(value, list) and value and isinstance(value[0], list):
                rows = value
            else:
                continue
            heading = f"{algorithm} \u2014 {table_name}" if algorithm else table_name
            parts.append("\\subsection{" + escape_latex(heading) + "}\n")
            parts.append(latex_long_table(rows, header=True))
            parts.append("\n\\vspace{0.6em}\n")
    return "".join(parts)


def _results_tsv_tex(tsv_path: Path) -> str:
    if not tsv_path.exists():
        return ""
    with tsv_path.open(newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if not rows:
        return ""
    body_rows = rows[1:]
    truncated = len(body_rows) > _MAX_TABLE_ROWS
    display_rows = [rows[0]] + body_rows[:_MAX_TABLE_ROWS]
    parts = [
        "\\subsection{Per-event results (results.tsv)}\n",
        latex_long_table(display_rows, header=True),
    ]
    if truncated:
        parts.append(
            "\n\\par\\textit{Showing the first "
            f"{_MAX_TABLE_ROWS} of {len(body_rows)} rows; "
            "see results.tsv for the complete table.}\n"
        )
    return "".join(parts)


def _figures_section_tex(visualization_dir: Path, collector: _FigureCollector) -> str:
    if not visualization_dir.exists():
        return ""
    svg_paths = sorted(visualization_dir.glob("*.svg"))
    if not svg_paths:
        return ""
    parts = ["\\section{Figures}\n"]
    for svg_path in svg_paths:
        name = collector.add(svg_path)
        if name is None:
            continue
        parts.append(_figure_block(name, svg_path.stem.replace("_", " ")))
    return "".join(parts)


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

    Renders a real two-column academic-paper-style LaTeX document (title,
    single-column abstract, two-column Results-summary body, a landscape
    booktabs/xltabular Tables section, and a Figures section), compiled with
    Tectonic -- see :mod:`cfh.reporting.latex`. Every sentence comes from
    :func:`cfh.reporting.text.render_abstract` /
    :func:`cfh.reporting.text.render_results_summary`; nothing is computed
    here, and every piece of data-derived text is escaped via
    :func:`cfh.reporting.latex.escape_latex` before being embedded.
    """
    output_path = Path(output_path)
    run_dir = output_path.parent
    tsv_path = Path(results_tsv_path) if results_tsv_path else run_dir / "results.tsv"
    viz_dir = Path(visualizations_dir) if visualizations_dir else run_dir / "visualizations"

    gene = payload.get("gene_symbol") or "Unknown gene"
    study = payload.get("study_id") or "unknown study"

    with tempfile.TemporaryDirectory(prefix="cfh-latex-figures-") as scratch:
        collector = _FigureCollector(Path(scratch))

        parts = [
            _LATEX_PREAMBLE,
            _title_and_abstract_tex(
                f"{gene} fusion-hotspot benchmark report",
                f"Study: {study}",
                render_abstract(payload),
            ),
            "\\section{Results summary}\n",
        ]
        for section in render_results_summary(payload):
            parts.append("\\subsection{" + escape_latex(section["heading"]) + "}\n")
            parts.append(escape_latex(section["paragraph"]) + "\n")

        tables_tex = _algorithm_tables_tex(payload) + _results_tsv_tex(tsv_path)
        if tables_tex:
            parts.append(_landscape_table_section("Tables", tables_tex))

        parts.append(_figures_section_tex(viz_dir, collector))
        parts.append("\\end{document}\n")

        tex_source = "\n".join(parts)
        return compile_latex(tex_source, output_path, figures=collector.figures)


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

    Renders the same real two-column academic-paper-style LaTeX document
    (title, single-column abstract, two-column body, landscape
    booktabs/xltabular tables, figures) as :func:`render_pdf_report` --
    compiled with Tectonic, see :mod:`cfh.reporting.latex`. Nothing is
    computed here; every piece of data-derived text is escaped via
    :func:`cfh.reporting.latex.escape_latex` before being embedded.
    """
    output_path = Path(output_path)

    with tempfile.TemporaryDirectory(prefix="cfh-latex-figures-") as scratch:
        collector = _FigureCollector(Path(scratch))

        parts = [
            _LATEX_PREAMBLE,
            _title_and_abstract_tex(title, "", abstract),
            "\\section{Methods}\n",
            escape_latex(methods) + "\n",
            "\\section{Results}\n",
            "\\subsection{Genome-wide summary}\n",
            escape_latex(manhattan_caption) + "\n",
        ]
        if manhattan_svg_path is not None and Path(manhattan_svg_path).exists():
            name = collector.add(Path(manhattan_svg_path))
            if name is not None:
                parts.append(_figure_block(name, manhattan_caption))

        if results_table_rows:
            parts.append(
                _landscape_table_section(
                    "FDR-significant and honorable-mention genes",
                    latex_long_table(results_table_rows, header=True),
                )
            )

        parts.append("\\subsection{Gene highlights}\n")
        for highlight in gene_highlights:
            parts.append("\\subsubsection{" + escape_latex(highlight["heading"]) + "}\n")
            parts.append(escape_latex(highlight["paragraph"]) + "\n")
            figure_path = highlight.get("figure_path")
            if figure_path is not None and Path(figure_path).exists():
                name = collector.add(Path(figure_path))
                if name is not None:
                    parts.append(_figure_block(name, highlight.get("figure_caption") or ""))
            if highlight.get("report_note"):
                parts.append("\\par\\textit{" + escape_latex(highlight["report_note"]) + "}\n")

        parts.append("\\section{Discussion}\n")
        parts.append("\\begin{itemize}\n")
        for bullet in discussion_bullets:
            parts.append("\\item " + escape_latex(bullet) + "\n")
        parts.append("\\end{itemize}\n")

        if appendix_rows:
            parts.append(
                _landscape_table_section(
                    "Appendix: per-gene report index",
                    latex_long_table(appendix_rows, header=True),
                )
            )

        parts.append("\\end{document}\n")

        tex_source = "\n".join(parts)
        return compile_latex(tex_source, output_path, figures=collector.figures)


def render_pdf_report_for_run_dir(run_dir: str | Path) -> Path:
    """Convenience wrapper: render ``report.pdf`` for an on-disk run directory
    that already has ``results.json`` (and, if present, ``results.tsv`` /
    ``visualizations/``) written by ``write_outputs``.
    """
    run_dir = Path(run_dir)
    payload = json.loads((run_dir / "results.json").read_text())
    return render_pdf_report(payload, run_dir / "report.pdf")
