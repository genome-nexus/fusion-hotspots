"""Tests for the shared SVG text-escaping helper (:mod:`cfh.reporting.svg_utils`),
used by every hand-built SVG renderer that embeds data-derived text (gene
names, sample IDs, ...) into ``<title>``/``<text>`` element content.
"""

from __future__ import annotations

from cfh.reporting.svg_utils import escape_xml_text


def test_escapes_ampersand_and_angle_brackets():
    assert escape_xml_text("A&B<C>D") == "A&amp;B&lt;C&gt;D"


def test_ampersand_escaped_before_reintroducing_entities():
    # A naive multi-pass replace could double-escape; this input would
    # produce "&amp;lt;" if "&" were escaped after "<" instead of before.
    assert escape_xml_text("<") == "&lt;"
    assert escape_xml_text("&lt;") == "&amp;lt;"


def test_leaves_ordinary_text_and_quotes_unchanged():
    assert escape_xml_text("BRAF p.V600E (sample 'S1')") == "BRAF p.V600E (sample 'S1')"
    assert escape_xml_text('say "hi"') == 'say "hi"'
