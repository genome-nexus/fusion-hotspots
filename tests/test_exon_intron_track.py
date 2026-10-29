import json
import re

import pytest

from cfh.reporting.exon_intron_track import _kept, _slot_key, render_exon_intron_track_svg
from cfh.reporting.html_viewer import build_run_viewer
from conftest import latest_run_dir

RANKS = [2, 3, 4, 5]


def _event(role, exon, intronic, status="retained", **extra):
    return {
        "event_id": f"{role}-{exon}-{intronic}",
        "sample_id": "S",
        "partner_gene": "P",
        "target_role": role,
        "breakpoint_exon": exon,
        "is_intronic_breakpoint": intronic,
        "domain_status": status,
        **extra,
    }


@pytest.mark.parametrize(
    ("role", "exon", "intronic", "slot"),
    [
        # 5' gene snaps an intronic break to the last exon it keeps: the
        # break is in the intron after that exon.
        ("five_prime", 3, True, ("intron", 3)),
        # 3' gene snaps to the first exon it keeps: the intron before it.
        ("three_prime", 4, True, ("intron", 3)),
        ("three_prime", 2, True, ("intron", 0)),  # upstream of the first coding exon
        ("five_prime", 5, True, ("intron", 5)),  # downstream of the last
        ("five_prime", 3, False, ("exon", 3)),
        ("three_prime", 3, False, ("exon", 3)),
        ("other", 3, True, None),
        ("five_prime", None, True, None),
    ],
)
def test_breakpoints_land_in_the_right_exon_or_intron(role, exon, intronic, slot):
    assert _slot_key(_event(role, exon, intronic), RANKS) == slot


@pytest.mark.parametrize(
    ("role", "exon", "intronic", "whole", "partial"),
    [
        ("five_prime", 3, True, {2, 3}, set()),
        ("five_prime", 3, False, {2}, {3}),
        ("three_prime", 4, True, {4, 5}, set()),
        ("three_prime", 4, False, {5}, {4}),
    ],
)
def test_kept_exons_follow_the_gene_role(role, exon, intronic, whole, partial):
    assert _kept(_event(role, exon, intronic), RANKS) == (whole, partial)


def _payload(events):
    return {
        "gene_symbol": "TEST",
        "summary": {"key_domains": [{"name": "Kinase", "start_aa": 60, "end_aa": 140}]},
        "gene_track": {
            "exon_boundaries_aa": [
                {"exon_rank": rank, "start_aa": 1 + 50 * i, "end_aa": 50 * (i + 1)}
                for i, rank in enumerate(RANKS)
            ]
        },
        "events": events,
    }


def test_svg_counts_breaks_and_exon_inclusion():
    events = [
        _event("three_prime", 4, True),
        _event("three_prime", 4, True),
        _event("three_prime", 4, False, status="disrupted"),
        _event("five_prime", 3, True, status="lost", partner_gene="Q"),
    ]
    svg = render_exon_intron_track_svg(_payload(events))

    assert svg is not None and svg.startswith("<svg")
    titles = re.findall(r"<title>([^<]*)</title>", svg)
    # Two 3' breaks before exon 4 and one 5' break after exon 3 share intron 3.
    assert (
        "intron 3: 3 breakpoint(s) (retained 2, lost 1; as 5' gene 1, as 3' gene 2); "
        "top partners P 2, Q 1"
    ) in titles
    assert any(t.startswith("exon 4 (mid-exon): 1 breakpoint(s) (truncated 1") for t in titles)
    # Exon 5 is kept whole by all three 3' fusions; exon 4 by two, in part by one.
    assert "exon 5: kept whole by 3 of 4 fusions (75%)" in titles
    assert "exon 4: kept whole by 2 of 4 fusions (50%), in part by 1" in titles
    assert "exon 2: kept whole by 1 of 4 fusions (25%)" in titles
    assert "Kinase" in svg


def test_svg_is_omitted_without_exon_boundaries_or_mappable_breaks():
    assert render_exon_intron_track_svg({"events": [_event("five_prime", 3, True)]}) is None
    assert render_exon_intron_track_svg(_payload([])) is None
    assert render_exon_intron_track_svg(_payload([_event("other", 3, True)])) is None


def _write_run(tmp_path, payload, *, protein_svg=True):
    run_dir = tmp_path / "run"
    (run_dir / "visualizations").mkdir(parents=True)
    (run_dir / "results.json").write_text(json.dumps(payload))
    if protein_svg:
        (run_dir / "visualizations" / "domain_retention_outliers.svg").write_text(
            '<svg xmlns="http://www.w3.org/2000/svg"><circle data-event-id="E1" r="3"/></svg>'
        )
    return run_dir


def test_viewer_switches_between_protein_and_exon_views(tmp_path):
    html = build_run_viewer(_write_run(tmp_path, _payload([_event("three_prime", 4, True)])))
    html = html.read_text()

    assert re.findall(r'data-track-view="(\w+)"', html) == ["protein", "exon"]
    assert '<div data-track-panel="protein">' in html
    assert '<div data-track-panel="exon" hidden>' in html
    # The protein view keeps its id so click-to-event wiring still finds it.
    assert 'id="domain-track-svg"' in html
    assert 'id="exon-intron-svg"' in html
    assert 'querySelectorAll("[data-track-panel]")' in html


def test_viewer_shows_the_exon_view_alone_without_a_protein_track(tmp_path):
    payload = _payload([_event("three_prime", 4, True)])
    html = build_run_viewer(_write_run(tmp_path, payload, protein_svg=False)).read_text()
    assert '<button type="button" role="tab"' not in html
    assert 'id="exon-intron-svg"' in html


@pytest.mark.parametrize(
    ("gene", "slot_title"),
    [
        ("alk", "intron 19: 236 breakpoint(s)"),
        ("fgfr2", "intron 17: 125 breakpoint(s)"),
    ],
)
def test_committed_run_puts_textbook_junctions_in_the_right_intron(gene, slot_title):
    """ALK's kinase-retaining break is in intron 19 and FGFR2's in intron 17."""
    reports = latest_run_dir("cohort-scan_msk_impact_50k_2026") / "cohort_scan" / "gene_reports"
    path = reports / gene / "results.json"
    if not path.exists():
        pytest.skip("committed cohort scan has no report for this gene")
    svg = render_exon_intron_track_svg(json.loads(path.read_text()))
    assert svg is not None
    peak = max(
        (title for title in re.findall(r"<title>([^<]*)</title>", svg) if "breakpoint(s)" in title),
        key=lambda title: int(re.search(r": (\d+) breakpoint", title).group(1)),
    )
    assert peak.startswith(slot_title)


def _with_genomic(payload, strand=1):
    # Exons 50 bp apart except a 10 kb intron after exon 3 (ranks 2-5).
    starts = {2: 1000, 3: 1200, 4: 11400, 5: 11600}
    exons = [{"exon_rank": r, "start": s, "end": s + 149} for r, s in starts.items()]
    if strand == -1:
        exons = [
            {"exon_rank": e["exon_rank"], "start": -e["end"], "end": -e["start"]} for e in exons
        ]
    payload["gene_track"]["strand"] = strand
    payload["gene_track"]["exons_genomic"] = exons
    return payload


@pytest.mark.parametrize("strand", [1, -1])
def test_intron_lengths_come_from_exon_coordinates(strand):
    from cfh.reporting.exon_intron_track import _intron_lengths

    track = _with_genomic(_payload([]), strand)["gene_track"]
    assert _intron_lengths(track, RANKS) == {2: 50, 3: 10050, 4: 50}
    assert _intron_lengths({}, RANKS) == {}


def _intron_rects(svg):
    """Width of each dashed intron band in the gene-structure row."""
    pattern = (
        r'<g><title>(intron [^<:,]+)[^<]*</title><rect x="[\d.]+" y="[\d.]+" '
        r'width="([\d.]+)"[^>]*stroke-dasharray'
    )
    return {m.group(1): float(m.group(2)) for m in re.finditer(pattern, svg)}


def test_introns_are_drawn_and_sized_by_log_length_when_known():
    events = [_event("three_prime", 4, True), _event("three_prime", 4, True)]
    svg = render_exon_intron_track_svg(_with_genomic(_payload(events)))
    widths = _intron_rects(svg)

    assert set(widths) == {"intron 2", "intron 3", "intron 4"}
    assert widths["intron 3"] > widths["intron 2"] == widths["intron 4"]
    assert "intron 3, 10.1 kb: 2 breakpoint(s)" in svg
    assert "intron 2, 50 bp: no breakpoints" in svg
    assert ">i3</text>" in svg  # the intron holding breakpoints is named
    assert "intron width follows log genomic length" in svg


def test_introns_without_coordinates_get_a_visible_fixed_width():
    svg = render_exon_intron_track_svg(_payload([_event("three_prime", 4, True)]))
    widths = _intron_rects(svg)
    assert set(widths) == {"intron 2", "intron 3", "intron 4"}
    assert len(set(widths.values())) == 1 and next(iter(widths.values())) >= 10
    assert "intron lengths not recorded in this run" in svg


def test_gene_track_records_genomic_exons_and_strand(
    genome_nexus_canonical_transcript_fixture_path,
):
    from cfh.mapping.genome_nexus_source import (
        gene_track_from_canonical_transcript,
        parse_canonical_transcript,
    )

    canonical = parse_canonical_transcript(
        json.loads(genome_nexus_canonical_transcript_fixture_path.read_text())
    )
    track = gene_track_from_canonical_transcript(canonical)
    assert track["strand"] == canonical.exons[0].strand
    assert [e["exon_rank"] for e in track["exons_genomic"]] == sorted(
        e.rank for e in canonical.exons
    )
    first = min(canonical.exons, key=lambda e: e.rank)
    assert track["exons_genomic"][0] == {
        "exon_rank": first.rank,
        "start": first.start,
        "end": first.end,
    }
