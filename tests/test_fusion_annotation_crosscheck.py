import json

import pytest

from cfh.mapping import fusion_annotation_crosscheck as fac
from cfh.mapping import genome_nexus_source as gns


def _braf_canonical(genome_nexus_canonical_transcript_fixture_path):
    payload = json.loads(genome_nexus_canonical_transcript_fixture_path.read_text())
    return gns.parse_canonical_transcript(payload)


def test_interior_coding_exon_breakpoint_agrees_with_our_own_resolver(
    genome_nexus_canonical_transcript_fixture_path,
):
    """Sanity baseline: an interior (non-boundary, non-intronic) coding-exon
    breakpoint must produce the same residue from both the existing resolver
    and the ported fusion-annotation arithmetic, for at least one fusion
    role -- proving the two independently-implemented strand/CDS-offset
    arithmetics genuinely agree on ordinary input, not just on error paths.
    """
    canonical = _braf_canonical(genome_nexus_canonical_transcript_fixture_path)
    exon5 = next(e for e in canonical.exons if e.rank == 5)
    breakpoint_genomic = exon5.start + 62
    strand = canonical.exons[0].strand
    cds_min, cds_max = gns.cds_bounds_from_utrs(canonical.utrs)

    ours = gns.map_genomic_breakpoint_to_protein_position(
        canonical.exons,
        breakpoint_genomic,
        strand,
        cds_min_genomic=cds_min,
        cds_max_genomic=cds_max,
    )
    assert ours.is_intronic is False
    assert ours.protein_position == 217

    result = fac.crosscheck_breakpoint_protein_position(
        canonical.exons,
        breakpoint_genomic,
        strand,
        "three_prime",
        ours.protein_position,
        cds_min_genomic=cds_min,
        cds_max_genomic=cds_max,
    )

    assert result.error is None
    assert result.agrees is True
    assert result.protein_position == 217
    assert result.note == "exact match"
    assert result.is_hybrid_codon is False
    assert result.exon_rank == 5


def test_hybrid_codon_breakpoint_is_a_known_convention_difference_not_a_disagreement(
    genome_nexus_canonical_transcript_fixture_path,
):
    """The same breakpoint from the five_prime role's perspective lands on a
    hybrid (non-codon-boundary) CDS position, so fusion-annotation's "last
    complete 5' residue" is exactly one less than this project's
    hybrid-junction-inclusive residue -- expected, and must be reported as
    agreeing, not as a disagreement.
    """
    canonical = _braf_canonical(genome_nexus_canonical_transcript_fixture_path)
    exon5 = next(e for e in canonical.exons if e.rank == 5)
    breakpoint_genomic = exon5.start + 62
    strand = canonical.exons[0].strand
    cds_min, cds_max = gns.cds_bounds_from_utrs(canonical.utrs)

    result = fac.crosscheck_breakpoint_protein_position(
        canonical.exons,
        breakpoint_genomic,
        strand,
        "five_prime",
        217,
        cds_min_genomic=cds_min,
        cds_max_genomic=cds_max,
    )

    assert result.error is None
    assert result.agrees is True
    assert result.protein_position == 216
    assert result.is_hybrid_codon is True
    assert "hybrid" in result.note


def test_genuine_disagreement_is_reported_as_such(genome_nexus_canonical_transcript_fixture_path):
    canonical = _braf_canonical(genome_nexus_canonical_transcript_fixture_path)
    exon5 = next(e for e in canonical.exons if e.rank == 5)
    breakpoint_genomic = exon5.start + 62
    strand = canonical.exons[0].strand
    cds_min, cds_max = gns.cds_bounds_from_utrs(canonical.utrs)

    result = fac.crosscheck_breakpoint_protein_position(
        canonical.exons,
        breakpoint_genomic,
        strand,
        "three_prime",
        999,
        cds_min_genomic=cds_min,
        cds_max_genomic=cds_max,
    )

    assert result.error is None
    assert result.agrees is False
    assert "disagreement" in result.note


def test_missing_utr_bounds_falls_back_to_outermost_exon_coordinate(
    genome_nexus_canonical_transcript_fixture_path,
):
    """A transcript can have no annotated UTR of one/both types at all (this
    project's own cds_bounds_from_utrs() then reports that bound as
    ``None``, e.g. NTRK1's canonical transcript has no three_prime_UTR).
    The ported fusion-annotation arithmetic needs concrete integer bounds,
    so the cross-check falls back to the transcript's outermost exon
    coordinate on that side -- the same fallback fusion-annotation's own
    gn_provider.py documents for this exact case -- rather than refusing to
    run at all.
    """
    canonical = _braf_canonical(genome_nexus_canonical_transcript_fixture_path)
    exon5 = next(e for e in canonical.exons if e.rank == 5)
    breakpoint_genomic = exon5.start + 62
    strand = canonical.exons[0].strand

    result = fac.crosscheck_breakpoint_protein_position(
        canonical.exons,
        breakpoint_genomic,
        strand,
        "three_prime",
        217,
        cds_min_genomic=None,
        cds_max_genomic=None,
    )

    # With no UTR bounds to clip against at all, the fallback uses the raw
    # exon extents -- so this necessarily produces a *different* (UTR-naive)
    # CDS position than the properly-clipped 217, exactly the kind of gap
    # this project's own cds_bounds_from_utrs() exists to avoid when real
    # UTR data is available. The point of this test is that the cross-check
    # still runs and returns a real, comparable value instead of an error.
    assert result.error is None
    assert result.agrees is not None
    assert result.protein_position is not None


def test_breakpoint_outside_the_cds_entirely_degrades_gracefully(
    genome_nexus_canonical_transcript_fixture_path,
):
    """A three_prime breakpoint past the far end of every coding exon (here,
    genomic position 1 -- far 3' of BRAF's minus-strand transcript) cannot
    map into the CDS at all; the ported cds_coord_at_genomic raises
    ValueError for exactly this case, which must be caught and reported via
    ``error``, never propagated.
    """
    canonical = _braf_canonical(genome_nexus_canonical_transcript_fixture_path)
    strand = canonical.exons[0].strand
    cds_min, cds_max = gns.cds_bounds_from_utrs(canonical.utrs)

    result = fac.crosscheck_breakpoint_protein_position(
        canonical.exons,
        1,
        strand,
        "three_prime",
        1,
        cds_min_genomic=cds_min,
        cds_max_genomic=cds_max,
    )

    assert result.agrees is None
    assert result.error is not None
    assert result.protein_position is None


def test_unsupported_role_degrades_gracefully(genome_nexus_canonical_transcript_fixture_path):
    canonical = _braf_canonical(genome_nexus_canonical_transcript_fixture_path)
    cds_min, cds_max = gns.cds_bounds_from_utrs(canonical.utrs)

    result = fac.crosscheck_breakpoint_protein_position(
        canonical.exons,
        140_000_000,
        canonical.exons[0].strand,
        "not_a_real_role",
        123,
        cds_min_genomic=cds_min,
        cds_max_genomic=cds_max,
    )

    assert result.agrees is None
    assert result.error is not None
    assert "unsupported fusion role" in result.error


def test_no_exon_data_degrades_gracefully():
    result = fac.crosscheck_breakpoint_protein_position(
        [],
        140_000_000,
        -1,
        "five_prime",
        123,
        cds_min_genomic=100,
        cds_max_genomic=200,
    )

    assert result.agrees is None
    assert result.error is not None


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        (117324415, 117324415),
        ("117324415", 117324415),
        ("chr6:117324415", 117324415),
        ("6:117324415", 117324415),
        ("g.117324415", 117324415),
        ("chr6:g.117324415", 117324415),
        ("117,324,415", 117324415),
    ],
)
def test_parse_genomic_breakpoint_accepts_documented_formats(value, expected):
    assert fac.parse_genomic_breakpoint(value) == expected


def test_parse_genomic_breakpoint_rejects_unparseable_value():
    with pytest.raises(ValueError):
        fac.parse_genomic_breakpoint("not-a-position")


def test_build_exon_genomic_map_and_cds_map_agree_on_exon_ordering(
    genome_nexus_canonical_transcript_fixture_path,
):
    """BRAF is minus-strand; both maps must independently order exons
    5'->3' by descending genomic start, matching this project's own
    rank-based exon ordering.
    """
    canonical = _braf_canonical(genome_nexus_canonical_transcript_fixture_path)
    strand = canonical.exons[0].strand
    cds_min, cds_max = gns.cds_bounds_from_utrs(canonical.utrs)
    exon_dicts = [{"start": e.start, "end": e.end} for e in canonical.exons]

    genomic_map = fac.build_exon_genomic_map(strand, exon_dicts)
    cds_map = fac.build_exon_cds_map(strand, exon_dicts, cds_min, cds_max)

    assert len(genomic_map) == len(exon_dicts) == len(cds_map)
    by_rank = sorted(canonical.exons, key=lambda e: e.rank)
    expected_order = sorted((e.start, e.end) for e in by_rank)[::-1]
    assert genomic_map == expected_order
