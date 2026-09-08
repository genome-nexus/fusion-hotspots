"""Unit tests for cfh.reporting.exon_labels."""

from __future__ import annotations

from cfh.reporting.exon_labels import (
    exon_label_for_protein_position,
    format_breakpoint_exon,
)


def test_format_breakpoint_exon_states_unavailable_when_missing():
    assert format_breakpoint_exon(None) == "exon unavailable"
    assert format_breakpoint_exon(None, is_intronic=True) == "exon unavailable"


def test_format_breakpoint_exon_plain_exon_when_not_intronic():
    assert format_breakpoint_exon(10) == "exon 10"
    assert format_breakpoint_exon(10, is_intronic=False) == "exon 10"
    assert format_breakpoint_exon(10, is_intronic=None) == "exon 10"


def test_format_breakpoint_exon_flags_intronic_approximation():
    assert format_breakpoint_exon(8, is_intronic=True) == "intronic breakpoint, nearest exon 8"


BRAF_EXON_BOUNDARIES = [
    {"exon_rank": 7, "start_aa": 287, "end_aa": 327},
    {"exon_rank": 8, "start_aa": 327, "end_aa": 380},
    {"exon_rank": 9, "start_aa": 380, "end_aa": 421},
]


def test_exon_label_for_protein_position_inside_one_exon():
    assert exon_label_for_protein_position(BRAF_EXON_BOUNDARIES, 350) == "exon 8"


def test_exon_label_for_protein_position_at_shared_boundary_names_both_exons():
    # BRAF's real committed run: aa 327 is exactly exon 7's end and exon 8's
    # start (a clamped intronic-breakpoint estimate artifact).
    assert exon_label_for_protein_position(BRAF_EXON_BOUNDARIES, 327) == "exon 7/8 boundary"
    assert exon_label_for_protein_position(BRAF_EXON_BOUNDARIES, 380) == "exon 8/9 boundary"


def test_exon_label_for_protein_position_unavailable_cases():
    unavailable = "exon position unavailable"
    assert exon_label_for_protein_position(None, 100) == unavailable
    assert exon_label_for_protein_position([], 100) == unavailable
    assert exon_label_for_protein_position(BRAF_EXON_BOUNDARIES, None) == unavailable
    # Position outside every configured exon span (e.g. downstream of the
    # last exon's end).
    assert exon_label_for_protein_position(BRAF_EXON_BOUNDARIES, 9999) == unavailable
