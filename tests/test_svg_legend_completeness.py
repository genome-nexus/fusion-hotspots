"""Every semantic color drawn anywhere in a renderer's SVG body must have a
documented meaning in that same SVG (a legend swatch+label, or -- for a
color whose exact per-instance shade is arbitrary, like a per-partner or
extra-key-domain highlight -- an explanatory note covering the whole
scheme). See ``cfh.reporting.palette``'s module docstring and
``SEMANTIC_COLORS``/``CHROME_COLORS``.

Uses the real, already-committed BRAF and RET benchmark artifacts under
``runs/`` (no network access, no synthetic fixtures) so this checks the
renderers' actual real-world output, not a fixture crafted to pass.
"""

from __future__ import annotations

import re
from pathlib import Path

from cfh.reporting import palette as palette_module
from cfh.reporting.palette import CHROME_COLORS, SEMANTIC_COLORS

REPO_ROOT = Path(__file__).parent.parent
RUN_DIRS = {
    "BRAF": REPO_ROOT / "runs" / "braf_msk-impact-50k-2026_20260905T011754Z",
    "RET": REPO_ROOT / "runs" / "ret_msk-impact-50k-2026_20260905T011503Z",
}

_SVG_NAMES = [
    "domain_retention_outliers.svg",
    "reference_comparison.svg",
    "fusion_schematic.svg",
    "intragenic_deletion_schematic.svg",
]

_COLOR_ATTR_RE = re.compile(r'(?:fill|stroke)="(#[0-9a-fA-F]{6})"')
_SWATCH_WITH_LABEL_RE = re.compile(
    r'<(?:rect|circle)[^>]*(?:fill|stroke)="(#[0-9a-fA-F]{6})"[^>]*/>\s*'
    r"<text[^>]*>[^<]+</text>"
)


def _committed_svgs() -> dict[str, str]:
    svgs = {}
    for gene, run_dir in RUN_DIRS.items():
        for name in _SVG_NAMES:
            path = run_dir / "visualizations" / name
            if path.exists():
                svgs[f"{gene}/{name}"] = path.read_text()
    return svgs


def _used_colors(svg: str) -> set[str]:
    return {match.group(1) for match in _COLOR_ATTR_RE.finditer(svg)}


def _documented_colors(svg: str) -> set[str]:
    """Colors that appear as a swatch (``<rect>``/``<circle>``) immediately
    followed by a text label -- true both of a dedicated legend entry and
    of an inline per-item label (e.g. a domain-highlight rect followed by
    its own name/range text)."""
    return {match.group(1) for match in _SWATCH_WITH_LABEL_RE.finditer(svg)}


# fusion_schematic.svg draws one arbitrary hash-derived shade per fusion
# partner (see partner_color) -- by construction not any fixed hex value,
# so it can't be listed in SEMANTIC_COLORS/CHROME_COLORS. The renderer
# documents that whole *scheme* via an explanatory note (checked in
# test_semantic_colors_used_in_fusion_schematic_are_all_legended), not a
# legend entry per shade, so it's exempted from the "every color is a
# known palette constant" check below.
_HASH_COLOR_SVG_NAMES = {"fusion_schematic.svg"}


def test_every_committed_svg_has_at_least_one_swatch_and_uses_only_known_colors():
    svgs = _committed_svgs()
    assert len(svgs) >= 6, "expected both BRAF and RET to contribute several SVGs"
    known = set(SEMANTIC_COLORS) | set(CHROME_COLORS)
    for name, svg in svgs.items():
        used = _used_colors(svg)
        assert used, name
        if name.split("/", 1)[1] in _HASH_COLOR_SVG_NAMES:
            continue
        # Every color drawn is a documented palette constant -- never an
        # unrelated stray hex literal.
        unknown = used - known
        assert not unknown, (name, unknown)


def test_semantic_colors_used_in_domain_retention_track_are_all_legended():
    for gene, run_dir in RUN_DIRS.items():
        svg = (run_dir / "visualizations" / "domain_retention_outliers.svg").read_text()
        used_semantic = _used_colors(svg) & set(SEMANTIC_COLORS)
        assert used_semantic, gene  # sanity: this run does use semantic colors
        documented = _documented_colors(svg)
        assert used_semantic <= documented, (gene, used_semantic - documented)


def test_semantic_colors_used_in_reference_comparison_are_all_legended():
    for gene, run_dir in RUN_DIRS.items():
        svg = (run_dir / "visualizations" / "reference_comparison.svg").read_text()
        used_semantic = _used_colors(svg) & set(SEMANTIC_COLORS)
        assert used_semantic, gene
        documented = _documented_colors(svg)
        assert used_semantic <= documented, (gene, used_semantic - documented)


def test_semantic_colors_used_in_fusion_schematic_are_all_legended():
    for gene, run_dir in RUN_DIRS.items():
        path = run_dir / "visualizations" / "fusion_schematic.svg"
        if not path.exists():
            continue
        svg = path.read_text()
        used_semantic = _used_colors(svg) & set(SEMANTIC_COLORS)
        assert used_semantic, gene
        documented = _documented_colors(svg)
        assert used_semantic <= documented, (gene, used_semantic - documented)
        # Partner colors are real colors drawn in the body that are NOT in
        # SEMANTIC_COLORS/CHROME_COLORS (arbitrary per-partner hash) --
        # confirm the renderer explains that decorative scheme rather than
        # silently leaving them unaccounted for.
        assert "decorative" in svg


def test_semantic_colors_used_in_intragenic_deletion_schematic_are_all_legended():
    for gene, run_dir in RUN_DIRS.items():
        path = run_dir / "visualizations" / "intragenic_deletion_schematic.svg"
        if not path.exists():
            continue
        svg = path.read_text()
        used_semantic = _used_colors(svg) & set(SEMANTIC_COLORS)
        assert used_semantic, gene
        documented = _documented_colors(svg)
        assert used_semantic <= documented, (gene, used_semantic - documented)
        # The connector line's meaning is given by the schematic's own
        # explanatory subtitle rather than a legend swatch.
        assert "plain connector" in svg


def test_palette_categorizes_every_color_constant_as_semantic_or_chrome():
    """Regression guard: a newly added color constant in palette.py must be
    filed into SEMANTIC_COLORS or CHROME_COLORS, or this test (and the
    completeness tests above) can silently stop covering it."""
    all_hex_constants = {
        value
        for name, value in vars(palette_module).items()
        if name.isupper() and isinstance(value, str) and value.startswith("#")
    }
    categorized = set(SEMANTIC_COLORS) | set(CHROME_COLORS)
    assert all_hex_constants <= categorized, all_hex_constants - categorized
