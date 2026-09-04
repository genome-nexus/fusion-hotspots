"""LaTeX document assembly and Tectonic-based PDF compilation.

This module owns exactly two things: (1) escaping arbitrary data-derived
text for safe inclusion in ``.tex`` source, and (2) invoking the Tectonic
LaTeX engine to turn ``.tex`` source (plus any figure PDFs) into a real
compiled PDF. It computes no report content itself -- every sentence, table
cell, and number rendered through this module is produced elsewhere (see
:mod:`cfh.reporting.text` and :mod:`cfh.reporting.manuscript_text`) and
handed to :func:`escape_latex`/:func:`latex_table` here only for safe
typesetting.

Tectonic (https://tectonic-typesetting.io/) was chosen over a traditional
TeX Live install because it is a single self-contained, statically linked
binary (no multi-GB package install) that fetches only the LaTeX packages a
given document actually uses, on demand, from its own bundle server, and
caches them locally -- this was verified directly (not assumed) against a
real two-column ``article``-class document using ``booktabs``, ``graphicx``,
``pdflscape``, and ``xltabular`` before being adopted here.

If the ``tectonic`` binary is not on ``PATH`` at render time, or the
document fails to compile, this module fails loudly with an actionable
:class:`LatexRenderError` rather than silently falling back to a
different-looking PDF -- a report that renders successfully but with the
wrong layout is a worse failure mode than a build that stops and says so.
"""

from __future__ import annotations

import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any

from reportlab.graphics import renderPDF
from svglib.svglib import svg2rlg

_TECTONIC_BINARY = "tectonic"

_TECTONIC_INSTALL_HINT = (
    "Install Tectonic (a self-contained LaTeX engine -- no TeX Live required) "
    "to render PDF reports:\n"
    "  - macOS (Homebrew): brew install tectonic\n"
    "  - Linux/CI: download a prebuilt static binary from "
    "https://github.com/tectonic-typesetting/tectonic/releases "
    "(the *-x86_64-unknown-linux-musl.tar.gz asset needs no shared libraries) "
    "and put it on PATH\n"
    "  - or: cargo install tectonic\n"
    "See https://tectonic-typesetting.github.io/en-US/install.html for full "
    "instructions."
)


class LatexRenderError(RuntimeError):
    """Raised when the LaTeX toolchain is unavailable or a document fails to
    compile -- deliberately fatal (see the module docstring): a PDF report
    that silently degrades to a different renderer/layout when Tectonic is
    missing would be far more confusing to a reviewer than a loud failure
    that says exactly what to install.
    """


# Every character LaTeX treats as markup outside of math mode. Backslash
# MUST be handled in the same single regex pass as the rest (not a separate
# ``str.replace``) -- ``re.sub`` scans only the ORIGINAL string once, so a
# backslash produced by escaping (e.g. ``\%``) is never re-escaped, which a
# naive sequential ``str.replace("%", "\\%").replace("\\", ...)`` would get
# wrong in the opposite order.
_ESCAPE_MAP = {
    "\\": r"\textbackslash{}",
    "&": r"\&",
    "%": r"\%",
    "$": r"\$",
    "#": r"\#",
    "_": r"\_",
    "{": r"\{",
    "}": r"\}",
    "~": r"\textasciitilde{}",
    "^": r"\textasciicircum{}",
}
_ESCAPE_RE = re.compile("|".join(re.escape(char) for char in _ESCAPE_MAP))


def escape_latex(value: Any) -> str:
    """Escape a plain data value for safe inclusion in LaTeX source.

    ``None`` becomes ``""``; everything else is coerced with ``str()`` and
    has every LaTeX-special character (``\\ & % $ # _ { } ~ ^``) replaced
    with its literal-text escape. This is the ONLY place in the LaTeX
    rendering path that is allowed to touch raw data strings -- every gene
    name, note, warning, and formatted number that reaches a ``.tex`` file
    must pass through here first.
    """
    if value is None:
        return ""
    return _ESCAPE_RE.sub(lambda match: _ESCAPE_MAP[match.group()], str(value))


def latex_table(rows: list[list[Any]], *, header: bool = True, align: str | None = None) -> str:
    """Render ``rows`` (a header row followed by data rows, when
    ``header=True``) as a compact booktabs-style ``tabular``. Every cell is
    escaped via :func:`escape_latex`; no numbers are computed or reformatted
    here -- ``rows`` is already the exact, final display data (the same
    already-string-formatted rows the Markdown table renderer consumes).
    """
    if not rows:
        return ""
    n_cols = max(len(row) for row in rows)
    col_spec = align or ("l" * n_cols)
    lines = [f"\\begin{{tabular}}{{{col_spec}}}", "\\toprule"]
    for index, row in enumerate(rows):
        cells = [escape_latex(cell) for cell in row]
        cells += [""] * (n_cols - len(cells))
        lines.append(" & ".join(cells) + r" \\")
        if header and index == 0:
            lines.append("\\midrule")
    lines.append("\\bottomrule")
    lines.append("\\end{tabular}")
    return "\n".join(lines)


def latex_long_table(rows: list[list[Any]], *, header: bool = True, font_size: str = "tiny") -> str:
    """Render ``rows`` as a wide, auto-wrapping, multi-page ``xltabular``
    (equal-width ``X`` columns that wrap long cell text, repeating the
    header row on every page) -- used for the full per-event ``results.tsv``
    table and other wide/many-row tables that a fixed-width single-page
    ``tabular`` cannot hold. Every cell is escaped via :func:`escape_latex`.
    """
    if not rows:
        return ""
    n_cols = max(len(row) for row in rows)
    col_spec = "|" + "X|" * n_cols
    lines = [
        f"\\begingroup\\{font_size}",
        "\\setlength{\\tabcolsep}{2pt}",
        f"\\begin{{xltabular}}{{\\linewidth}}{{{col_spec}}}",
        "\\toprule",
    ]
    data_rows = rows
    if header:
        header_cells = [escape_latex(cell) for cell in rows[0]]
        header_cells += [""] * (n_cols - len(header_cells))
        lines.append(" & ".join(header_cells) + r" \\")
        lines.append("\\midrule")
        lines.append("\\endhead")
        lines.append("\\bottomrule")
        lines.append("\\endfoot")
        data_rows = rows[1:]
    else:
        lines.append("\\endhead")
    for row in data_rows:
        cells = [escape_latex(cell) for cell in row]
        cells += [""] * (n_cols - len(cells))
        lines.append(" & ".join(cells) + r" \\")
    lines.append("\\end{xltabular}")
    lines.append("\\endgroup")
    return "\n".join(lines)


def svg_to_pdf(svg_path: str | Path, pdf_path: str | Path) -> Path | None:
    """Convert one SVG file to a real, standalone PDF file.

    Reuses this project's existing svglib/reportlab SVG-loading path
    (``svglib.svg2rlg`` -> a reportlab ``Drawing``) and simply exports that
    same in-memory ``Drawing`` straight to a PDF file with
    ``reportlab.graphics.renderPDF.drawToFile`` instead of adding a new
    SVG-rendering dependency (inkscape/cairosvg). Returns ``pdf_path`` on
    success, or ``None`` if the SVG has no renderable content (mirrors
    ``cfh.reporting.pdf._load_svg_drawing``'s own None-on-empty behavior) --
    the caller is responsible for omitting that figure from the document
    rather than referencing a PDF that was never written.
    """
    drawing = svg2rlg(str(svg_path))
    if drawing is None or not drawing.width or not drawing.height:
        return None
    pdf_path = Path(pdf_path)
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    renderPDF.drawToFile(drawing, str(pdf_path), showBoundary=False)
    return pdf_path


def _tectonic_binary() -> str:
    binary = shutil.which(_TECTONIC_BINARY)
    if binary is None:
        raise LatexRenderError(
            "LaTeX rendering requires the 'tectonic' binary, which was not found "
            f"on PATH.\n{_TECTONIC_INSTALL_HINT}"
        )
    return binary


def compile_latex(
    tex_source: str,
    output_pdf_path: str | Path,
    *,
    figures: dict[str, Path] | None = None,
) -> Path:
    """Compile ``tex_source`` with Tectonic and write the resulting PDF to
    ``output_pdf_path``.

    ``figures`` maps each plain relative filename referenced by an
    ``\\includegraphics{...}`` in ``tex_source`` (e.g. ``"fig0.pdf"`` -- no
    directory component, no LaTeX-special characters, so it never needs
    ``grffile`` or path escaping) to an already-rendered PDF file on disk
    (see :func:`svg_to_pdf`); each is copied into the compile working
    directory before Tectonic runs, so every asset Tectonic sees is a
    self-contained relative filename rather than an absolute path (Tectonic
    warns that absolute-path inputs make a build non-reproducible).

    Raises :class:`LatexRenderError` -- never silently degrades -- if the
    ``tectonic`` binary is unavailable or the document fails to compile.
    """
    binary = _tectonic_binary()
    output_pdf_path = Path(output_pdf_path)
    output_pdf_path.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="cfh-latex-") as tmp:
        work_dir = Path(tmp)
        for relative_name, source_pdf in (figures or {}).items():
            shutil.copyfile(source_pdf, work_dir / relative_name)

        tex_path = work_dir / "report.tex"
        tex_path.write_text(tex_source)
        out_dir = work_dir / "out"
        out_dir.mkdir()

        result = subprocess.run(
            [binary, "report.tex", "--outdir", "out"],
            cwd=work_dir,
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            tail = "\n".join((result.stdout + result.stderr).splitlines()[-80:])
            raise LatexRenderError(
                "Tectonic failed to compile the LaTeX report "
                f"(exit code {result.returncode}). Last output:\n{tail}"
            )

        produced = out_dir / "report.pdf"
        if not produced.exists():
            raise LatexRenderError(
                "Tectonic exited successfully but did not produce report.pdf. "
                f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
            )
        shutil.copyfile(produced, output_pdf_path)

    return output_pdf_path
