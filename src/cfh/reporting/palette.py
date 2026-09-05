"""Single source of truth for every reporting renderer's colors.

The semantic colors below require an in-figure legend or explanatory key in
every renderer that uses them. Neutral chrome colors are self-evident
structure (axes, labels, connectors, or decorative report styling) and do
not need legend entries -- even where a chrome color happens to share a hex
value with a semantic one (coincidence, not shared meaning).

Every SVG/PDF renderer in :mod:`cfh.reporting` and
:mod:`cfh.real_benchmark` draws from this palette rather than hardcoding
its own hex literals, so no two renderers can silently drift apart on what
a color means.
"""

from __future__ import annotations

import colorsys
import zlib

RETAINED_COLOR = "#2878b5"
"""A domain (or domain segment) that is fully retained."""

TRUNCATED_COLOR = "#f2a93b"
"""A domain (or domain segment) that is partially retained/truncated."""

LOST_COLOR = "#777777"
"""A domain that is fully lost."""

BREAKPOINT_COLOR = "#d62728"
"""Breakpoint marker; also used as the lollipop track's
reference-discrepancy outline color."""

DOMAIN_HIGHLIGHT_COLOR = "#62b36f"
"""Configured key-domain span in the domain-retention lollipop track."""

BACKBONE_COLOR = "#e2e2e2"
"""Neutral protein backbone outside a drawn domain."""

REFERENCE_BAR_COLOR = "#999999"
"""Reference-cohort bar in a reference-versus-run comparison."""

AXIS_COLOR = "#444444"
"""Structural axes and protein backbones."""

GRID_COLOR = "#888888"
"""Structural threshold/grid lines and placeholder text."""

TEXT_COLOR = "#111111"
"""Primary text over a colored shape."""

SECONDARY_TEXT_COLOR = "#555555"
"""Secondary explanatory text and exon-boundary strokes."""

MUTED_TEXT_COLOR = "#666666"
"""Muted annotation text."""

EXON_TICK_COLOR = "#6b6b6b"
"""Structural exon-boundary tick marks. Deliberately its own shade rather
than reusing ``LOST_COLOR``'s grey, so the same hex value is never
ambiguous between "a lost domain" and "an exon boundary tick"."""

CONNECTOR_COLOR = "#a3a3a3"
"""Structural connector across an intragenic deleted span. Deliberately
its own shade rather than reusing ``REFERENCE_BAR_COLOR``'s grey, so the
same hex value is never ambiguous between "a reference-cohort bar" and "a
deleted-span connector"."""

TABLE_HEADER_COLOR = "#2878b5"
"""Decorative PDF report table header background. Coincidentally the same
blue as ``RETAINED_COLOR`` (for visual consistency with the SVG figures),
but purely decorative report styling -- not a domain-retention indicator,
so it needs no legend."""

SEMANTIC_COLORS = (
    RETAINED_COLOR,
    TRUNCATED_COLOR,
    LOST_COLOR,
    BREAKPOINT_COLOR,
    DOMAIN_HIGHLIGHT_COLOR,
    BACKBONE_COLOR,
    REFERENCE_BAR_COLOR,
)
"""Every color that requires a legend or explanatory key wherever it's
used. Test coverage (see ``tests/test_svg_legend_completeness.py``) checks
that any of these appearing in a rendered SVG is matched by a legend entry
in that same SVG. A gene-hash-derived per-domain highlight shade (see
:func:`cfh.real_benchmark._domain_highlight_color`) or per-partner shade
(:func:`cfh.reporting.fusion_schematic.partner_color`) is deliberately not
listed here -- individually arbitrary, but the *scheme* they belong to is
documented once by the corresponding legend/note text, not per hex value."""

CHROME_COLORS = (
    AXIS_COLOR,
    GRID_COLOR,
    TEXT_COLOR,
    SECONDARY_TEXT_COLOR,
    MUTED_TEXT_COLOR,
    EXON_TICK_COLOR,
    CONNECTOR_COLOR,
    TABLE_HEADER_COLOR,
)
"""Structural/decorative colors that need no legend entry (see the module
docstring). Listed explicitly, alongside ``SEMANTIC_COLORS``, so the two
sets together account for every color constant this module defines."""


def deterministic_color(label: str, *, lightness: float = 0.55, saturation: float = 0.55) -> str:
    """Deterministic, arbitrary-but-stable hex color for an arbitrary
    string label: the same label always gets the same color within a run
    and across runs (a pure function of the label), so a reader can
    visually track one entity across a diagram. Uses a CRC32 hash into
    hue space rather than Python's salted ``hash()``, which is randomized
    per-process and would make the same label render a different color on
    every regeneration.

    Used by :func:`cfh.reporting.fusion_schematic.partner_color` (each
    fusion partner gene gets its own stable, decorative color -- see that
    function's docstring for why partner colors need no legend) and by
    :func:`cfh.real_benchmark._domain_highlight_color` (a second or later
    configured key domain gets its own stable shade), so the two can't
    silently diverge on how a stable color is derived from a name.
    """
    hue = (zlib.crc32(label.encode("utf-8")) % 360) / 360.0
    r, g, b = colorsys.hls_to_rgb(hue, lightness, saturation)
    return f"#{int(r * 255):02x}{int(g * 255):02x}{int(b * 255):02x}"
