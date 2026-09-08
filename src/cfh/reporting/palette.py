"""Single source of truth for the SVG color conventions shared by every
renderer that visualizes them: domain-retention status, and the fixed
decorative fill used for fusion partner-gene blocks.

Two renderers currently draw from this palette:
:func:`cfh.real_benchmark._domain_track_svg` (the per-event domain-retention
lollipop track) and :mod:`cfh.reporting.fusion_schematic` (the
fusion-transcript schematic). Both must import these constants rather than
hardcoding their own hex literals, so the two visualizations can never
silently drift apart on what a color means.
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

PARTNER_COLOR = "#8064a2"
"""Fixed, non-semantic fill for every fusion partner-gene block in
:mod:`cfh.reporting.fusion_schematic`. Every partner-gene block renders in
this one color regardless of which gene it is -- it exists only to set the
partner block visually apart from the domain-colored target-gene block
sharing its row, not to distinguish one partner from another (row labels
do that, and there is deliberately no legend entry mapping this color to
any partner identity)."""


def deterministic_color(label: str, *, lightness: float = 0.55, saturation: float = 0.55) -> str:
    """Deterministic, arbitrary-but-stable hex color for an arbitrary
    string label: the same label always gets the same color within a run
    and across runs (a pure function of the label), so a reader can
    visually track one entity across a diagram. Uses a CRC32 hash into
    hue space rather than Python's salted ``hash()``, which is randomized
    per-process and would make the same label render a different color on
    every regeneration.

    Used by :func:`cfh.real_benchmark._domain_track_svg` for per-domain
    highlight coloring when a gene configures more than one key domain
    (each such domain gets its own legend entry, so distinct hues are the
    point). Fusion partner-gene blocks do *not* use this function -- see
    :func:`cfh.reporting.fusion_schematic.partner_color`, which is a fixed
    color instead, deliberately not derived from the partner's name.
    """
    hue = (zlib.crc32(label.encode("utf-8")) % 360) / 360.0
    r, g, b = colorsys.hls_to_rgb(hue, lightness, saturation)
    return f"#{int(r * 255):02x}{int(g * 255):02x}{int(b * 255):02x}"
