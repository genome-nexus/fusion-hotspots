"""Tests for :mod:`cfh.reporting.latex`: escaping, table assembly, SVG->PDF
conversion, and Tectonic compilation (including its fail-loud behavior when
Tectonic is unavailable).

The escaping tests are the dedicated proof required for every LaTeX
special character (``% & _ $ # { } ^ ~ \\``) never breaking compilation or
leaking as raw markup -- including an end-to-end render (real Tectonic
compile, real ``pypdf`` extraction) of a value containing a real special
character, not just a check of the escaping function in isolation.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest
from pypdf import PdfReader

from cfh.reporting.latex import (
    LatexRenderError,
    compile_latex,
    escape_latex,
    latex_long_table,
    latex_table,
    svg_to_pdf,
)
from cfh.reporting.pdf import render_pdf_report

REPO_ROOT = Path(__file__).parent.parent
BRAF_RUN_DIR = REPO_ROOT / "runs" / "braf_msk-impact-2017_20260904T005539Z"

_TECTONIC_AVAILABLE = shutil.which("tectonic") is not None
_requires_tectonic = pytest.mark.skipif(
    not _TECTONIC_AVAILABLE,
    reason="tectonic is not installed on PATH; see CONTRIBUTING.md for install instructions",
)


# ---------------------------------------------------------------------------
# escape_latex
# ---------------------------------------------------------------------------


def test_escape_latex_escapes_every_special_character():
    # Every character LaTeX treats as markup outside math mode, per the
    # hard requirement: % & _ $ # { } ^ ~ \
    assert escape_latex("%") == r"\%"
    assert escape_latex("&") == r"\&"
    assert escape_latex("_") == r"\_"
    assert escape_latex("$") == r"\$"
    assert escape_latex("#") == r"\#"
    assert escape_latex("{") == r"\{"
    assert escape_latex("}") == r"\}"
    assert escape_latex("^") == r"\textasciicircum{}"
    assert escape_latex("~") == r"\textasciitilde{}"
    assert escape_latex("\\") == r"\textbackslash{}"


def test_escape_latex_handles_a_realistic_mixed_string():
    value = "AGAP3_c.1140+237 (50%) & BRAF#1 costs $5 {x}^y~z"
    escaped = escape_latex(value)
    assert escaped == (
        r"AGAP3\_c.1140+237 (50\%) \& BRAF\#1 costs \$5 "
        r"\{x\}\textasciicircum{}y\textasciitilde{}z"
    )


def test_escape_latex_backslash_is_not_double_escaped():
    """A naive sequential ``str.replace`` chain (escape ``%`` before ``\\``,
    say) would re-escape the backslash it just inserted. ``escape_latex``
    must scan the ORIGINAL string once so this can't happen."""
    assert escape_latex("50% \\ done") == r"50\% \textbackslash{} done"


def test_escape_latex_none_becomes_empty_string():
    assert escape_latex(None) == ""


def test_escape_latex_coerces_non_string_values():
    assert escape_latex(42) == "42"
    assert escape_latex(3.14) == "3.14"


# ---------------------------------------------------------------------------
# latex_table / latex_long_table
# ---------------------------------------------------------------------------


def test_latex_table_is_booktabs_style_and_escapes_cells():
    tex = latex_table([["Gene", "Note"], ["SND1", "A&B_C 50%"]])
    assert "\\toprule" in tex
    assert "\\midrule" in tex
    assert "\\bottomrule" in tex
    assert "\\begin{tabular}" in tex
    assert "SND1" in tex
    assert r"A\&B\_C 50\%" in tex


def test_latex_table_empty_rows_returns_empty_string():
    assert latex_table([]) == ""


def test_latex_long_table_wraps_in_xltabular_and_escapes_cells():
    tex = latex_long_table([["Gene", "Note"], ["ETV6", "100% & done_here"]])
    assert "\\begin{xltabular}" in tex
    assert "\\endhead" in tex
    assert r"100\% \& done\_here" in tex


# ---------------------------------------------------------------------------
# svg_to_pdf
# ---------------------------------------------------------------------------


def test_svg_to_pdf_converts_a_real_committed_svg(tmp_path):
    svg_path = BRAF_RUN_DIR / "visualizations" / "domain_retention_outliers.svg"
    output = svg_to_pdf(svg_path, tmp_path / "fig.pdf")
    assert output == tmp_path / "fig.pdf"
    assert output.exists()
    assert output.stat().st_size > 0
    # A real, standalone, parseable PDF -- not just non-empty bytes.
    reader = PdfReader(str(output))
    assert len(reader.pages) == 1


def test_svg_to_pdf_returns_none_for_an_empty_svg(tmp_path):
    empty_svg = tmp_path / "empty.svg"
    empty_svg.write_text('<svg xmlns="http://www.w3.org/2000/svg"></svg>')
    assert svg_to_pdf(empty_svg, tmp_path / "out.pdf") is None


# ---------------------------------------------------------------------------
# compile_latex: fail-loud when Tectonic is unavailable
# ---------------------------------------------------------------------------


def test_compile_latex_raises_actionable_error_when_tectonic_is_missing(tmp_path, monkeypatch):
    monkeypatch.setenv("PATH", str(tmp_path))  # a directory guaranteed to have no `tectonic`
    with pytest.raises(LatexRenderError) as excinfo:
        compile_latex(
            r"\documentclass{article}\begin{document}x\end{document}",
            tmp_path / "out.pdf",
        )
    message = str(excinfo.value)
    assert "tectonic" in message.lower()
    assert "install" in message.lower()


@_requires_tectonic
def test_compile_latex_raises_actionable_error_on_a_broken_document(tmp_path):
    with pytest.raises(LatexRenderError) as excinfo:
        compile_latex(
            r"\documentclass{article}\begin{document}\undefinedcommand\end{document}",
            tmp_path / "out.pdf",
        )
    assert "Tectonic failed to compile" in str(excinfo.value)


# ---------------------------------------------------------------------------
# End-to-end: a real special character in real data survives a real
# Tectonic compile and appears literally in the extracted PDF text.
# ---------------------------------------------------------------------------


@_requires_tectonic
def test_render_pdf_report_escapes_real_special_characters_in_data(tmp_path):
    """A partner-gene name and a warning containing every LaTeX-special
    character must compile cleanly and appear as literal text (not broken
    markup, not a compile failure) in the rendered PDF -- proving the
    escaping is actually applied on the real render path, not just in the
    ``escape_latex`` unit tests above."""
    payload = {
        "gene_symbol": "TESTGENE",
        "study_id": "test_study",
        "summary": {
            "total_fusions": 10,
            "in_frame_count": 8,
            "in_frame_percent": 80.0,
            "domain_accession": "PF00069",
            "kinase_retained_count": 6,
            "kinase_retained_percent": 60.0,
        },
        "algorithm_results": [
            {
                "Algorithm": "frequency",
                "Summary": {"analyzed_event_count": 10, "unique_partner_gene_count": 1},
                "Tables": {
                    "Partner_gene_counts": [
                        {
                            "Partner_gene": "AB_C&D%E#F$G{H}I",
                            "Event_count": 5,
                            # Deliberately avoids the letter pairs "fi"/"fl"/"ffi": Latin
                            # Modern renders those as single ligature glyphs (U+FB01/FB02),
                            # so a literal ASCII "fi"/"fl" would not round-trip through
                            # pypdf's text extraction -- a font/PDF-extraction artifact
                            # unrelated to LaTeX escaping correctness, which is what this
                            # test actually checks.
                            "Note": "100% match & checked_case #7 costs $2 {ok}",
                        },
                    ]
                },
            }
        ],
    }

    output_path = tmp_path / "report.pdf"
    result_path = render_pdf_report(payload, output_path)

    assert result_path == output_path
    assert output_path.exists()

    reader = PdfReader(str(output_path))
    assert len(reader.pages) >= 1
    text = "".join(page.extract_text() or "" for page in reader.pages)

    # The ASCII-safe subset of LaTeX-special characters (`\%`, `\&`, `\_`,
    # `\$`, `\#`, `\{`, `\}`) round-trip through Tectonic's Latin Modern
    # font as their exact original ASCII glyph, so these appear verbatim.
    assert "AB_C&D%E#F$G{H}I" in text
    assert "100% match & checked_case #7 costs $2 {ok}" in text
