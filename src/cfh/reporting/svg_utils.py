"""Small shared helpers for hand-built inline SVG markup.

Currently just XML text-escaping, used by every renderer that embeds
data-derived text (gene symbols, sample IDs, partner names, ...) into SVG
``<text>``/``<title>`` element content -- :mod:`cfh.reporting.manhattan`,
:mod:`cfh.reporting.fusion_schematic`, and
:func:`cfh.real_benchmark._domain_track_svg`. A single shared function keeps
all of them escaping the same way rather than each carrying its own
possibly-drifting copy.
"""

from __future__ import annotations


def escape_xml_text(text: str) -> str:
    """Escape ``&``/``<``/``>`` for safe use as SVG element text content
    (e.g. inside ``<title>...</title>``).

    Quotes are intentionally left unescaped: this is for text *content*,
    not attribute values, so ``"``/``'`` are not XML-significant here.
    """
    return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
