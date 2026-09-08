"""Tests for the Genome Nexus canonical-transcript client and the
genomic-breakpoint-to-protein-position CDS-offset math.

Fixtures frozen from real responses:
  GET https://www.genomenexus.org/ensembl/canonical-transcript/hgnc/BRAF?isoformOverrideSource=mskcc
    -> tests/fixtures/genome_nexus/canonical_transcript_braf.json
  GET https://www.genomenexus.org/ensembl/transcript/ENST00000288602
    -> tests/fixtures/genome_nexus/transcript_ENST00000288602.json
  GET https://www.genomenexus.org/ensembl/canonical-transcript/hgnc/PIK3CA?isoformOverrideSource=mskcc
    -> tests/fixtures/genome_nexus/canonical_transcript_pik3ca.json
(all fetched 2026-09-02; the BRAF endpoints happen to return the same
transcript record shape for that gene, including exon coordinates. PIK3CA
is plus-strand with two separate 5' UTR segments -- exon 1
(178866311-178866391) is entirely 5' UTR, and exon 2 (178916538-178916965)
starts with a second 5' UTR segment (178916538-178916613) before its
coding portion begins -- used to test multi-segment UTR aggregation.)
"""

import json
import math
from unittest.mock import MagicMock

import pytest

from cfh.mapping import genome_nexus_source as gns
from cfh.mapping.domain_source import UniProtDomainSource


def _load(path):
    return json.loads(path.read_text())


def test_parse_canonical_transcript_from_real_fixture(
    genome_nexus_canonical_transcript_fixture_path,
):
    payload = _load(genome_nexus_canonical_transcript_fixture_path)
    canonical = gns.parse_canonical_transcript(payload)

    assert canonical.transcript_id == "ENST00000288602"
    assert canonical.refseq_mrna_id == "NM_004333"

    kinase_domains = [d for d in canonical.pfam_domains if d.pfam_id == "PF07714"]
    assert len(kinase_domains) == 1
    assert kinase_domains[0].start_aa == 458
    assert kinase_domains[0].end_aa == 712


def test_parse_canonical_transcript_skips_malformed_pfam_domain_entries():
    """Regression: real Genome Nexus canonical-transcript payloads for some
    genes (observed live on a genome-wide gene set, e.g. PALB2/BBC3/INPP4B)
    include an empty {} entry in pfamDomains alongside otherwise well-formed
    ones. Parsing must skip that entry, not raise KeyError and crash the
    whole gene."""
    payload = {
        "transcriptId": "ENST00000TEST",
        "refseqMrnaId": "NM_TEST",
        "proteinId": "ENSP00000TEST",
        "proteinLength": 100,
        "pfamDomains": [
            {},
            {"pfamDomainId": "PF00001", "pfamDomainStart": 10, "pfamDomainEnd": 50},
        ],
        "exons": [],
        "utrs": [],
    }

    canonical = gns.parse_canonical_transcript(payload)

    assert len(canonical.pfam_domains) == 1
    assert canonical.pfam_domains[0].pfam_id == "PF00001"


def test_parse_canonical_transcript_retains_utrs_from_payload(
    genome_nexus_canonical_transcript_fixture_path,
):
    """Regression: utrs must not be dropped -- they're required to derive
    correct CDS bounds (see cds_bounds_from_utrs)."""
    payload = _load(genome_nexus_canonical_transcript_fixture_path)
    canonical = gns.parse_canonical_transcript(payload)

    assert len(canonical.utrs) == len(payload["utrs"]) == 2
    utr_types = {u.utr_type for u in canonical.utrs}
    assert utr_types == {"five_prime_UTR", "three_prime_UTR"}


def test_cds_bounds_from_utrs_matches_independently_derived_bounds(
    genome_nexus_canonical_transcript_fixture_path,
):
    payload = _load(genome_nexus_canonical_transcript_fixture_path)
    canonical = gns.parse_canonical_transcript(payload)

    # Independently derive expected bounds straight from the raw fixture
    # (minus-strand gene: 5' UTR at the high end, 3' UTR at the low end).
    five_prime_utr = next(u for u in payload["utrs"] if u["type"] == "five_prime_UTR")
    three_prime_utr = next(u for u in payload["utrs"] if u["type"] == "three_prime_UTR")
    expected_cds_max = five_prime_utr["start"] - 1
    expected_cds_min = three_prime_utr["end"] + 1

    cds_min, cds_max = gns.cds_bounds_from_utrs(canonical.utrs)

    assert cds_min == expected_cds_min
    assert cds_max == expected_cds_max


def test_cds_bounds_from_utrs_returns_none_when_utr_missing():
    utrs = [gns.UtrRecord(utr_type="five_prime_UTR", start=100, end=200, strand=1)]
    cds_min, cds_max = gns.cds_bounds_from_utrs(utrs)
    assert cds_min == 201
    assert cds_max is None


def test_cds_bounds_from_utrs_aggregates_all_segments_of_a_type_plus_strand(
    genome_nexus_canonical_transcript_pik3ca_fixture_path,
):
    """Regression: PIK3CA (plus strand) has TWO separate 5' UTR segments
    (exon 1 is entirely 5' UTR; exon 2 starts with a second 5' UTR segment
    before its coding portion begins). Picking only the first segment
    encountered would derive cds_min from the wrong (inner) segment and
    miscount the second segment's 76nt as coding sequence.
    """
    payload = _load(genome_nexus_canonical_transcript_pik3ca_fixture_path)
    canonical = gns.parse_canonical_transcript(payload)

    five_prime_utrs = [u for u in payload["utrs"] if u["type"] == "five_prime_UTR"]
    three_prime_utrs = [u for u in payload["utrs"] if u["type"] == "three_prime_UTR"]
    assert len(five_prime_utrs) == 2  # confirms this fixture actually exercises the bug
    assert canonical.exons[0].strand == 1

    # Correct: adjacent to the OUTERMOST (highest-end) 5' UTR segment, and
    # the single 3' UTR segment.
    expected_cds_min = max(u["end"] for u in five_prime_utrs) + 1
    expected_cds_max = min(u["start"] for u in three_prime_utrs) - 1

    # The bug this guards against: deriving cds_min from only the FIRST
    # 5' UTR segment encountered in the (unordered) payload list.
    buggy_cds_min_from_first_segment_only = five_prime_utrs[0]["end"] + 1
    assert expected_cds_min != buggy_cds_min_from_first_segment_only

    cds_min, cds_max = gns.cds_bounds_from_utrs(canonical.utrs)

    assert cds_min == expected_cds_min == 178916614
    assert cds_max == expected_cds_max == 178952152
    assert cds_min != buggy_cds_min_from_first_segment_only


def test_pfam_domains_convert_to_protein_domain_shape(
    genome_nexus_canonical_transcript_fixture_path,
):
    payload = _load(genome_nexus_canonical_transcript_fixture_path)
    canonical = gns.parse_canonical_transcript(payload)
    domains = gns.pfam_domains_to_protein_domains(canonical)

    kinase = next(d for d in domains if d.name == "PF07714")
    assert kinase.start_aa == 458
    assert kinase.end_aa == 712
    assert kinase.source == "genome_nexus"


def test_exon_protein_boundaries_are_ordered_and_consistent_with_breakpoint_mapping(
    genome_nexus_canonical_transcript_fixture_path,
):
    """Every real breakpoint mapped inside one exon (via
    :func:`map_genomic_breakpoint_to_protein_position`) should land inside
    that same exon's reported protein-residue range from
    :func:`exon_protein_boundaries` -- both derive from the same
    coding-interval/codon-rounding arithmetic, so they must agree."""
    payload = _load(genome_nexus_canonical_transcript_fixture_path)
    canonical = gns.parse_canonical_transcript(payload)
    cds_min, cds_max = gns.cds_bounds_from_utrs(canonical.utrs)
    boundaries = gns.exon_protein_boundaries(
        canonical.exons, cds_min_genomic=cds_min, cds_max_genomic=cds_max
    )

    ranks = [b.exon_rank for b in boundaries]
    assert ranks == sorted(ranks)
    for boundary in boundaries:
        assert boundary.start_aa <= boundary.end_aa

    strand = canonical.exons[0].strand
    for exon in canonical.exons:
        interval = gns._coding_interval(exon, cds_min, cds_max)
        if interval is None:
            continue
        genomic_midpoint = (interval[0] + interval[1]) // 2
        mapped = gns.map_genomic_breakpoint_to_protein_position(
            canonical.exons,
            genomic_midpoint,
            strand,
            cds_min_genomic=cds_min,
            cds_max_genomic=cds_max,
        )
        boundary = next(b for b in boundaries if b.exon_rank == exon.rank)
        assert boundary.start_aa <= mapped.protein_position <= boundary.end_aa


def test_gene_track_from_canonical_transcript_matches_fixture(
    genome_nexus_canonical_transcript_fixture_path,
):
    payload = _load(genome_nexus_canonical_transcript_fixture_path)
    canonical = gns.parse_canonical_transcript(payload)
    track = gns.gene_track_from_canonical_transcript(canonical)

    assert track["protein_length"] == 766
    kinase = next(d for d in track["domains"] if d["accession"] == "PF07714")
    assert (kinase["start_aa"], kinase["end_aa"]) == (458, 712)
    assert len(track["exon_boundaries_aa"]) == len(canonical.exons)
    assert all(b["start_aa"] >= 1 for b in track["exon_boundaries_aa"])


def test_cds_offset_breakpoint_mapping_matches_hand_computed_arithmetic(
    genome_nexus_transcript_fixture_path,
):
    payload = _load(genome_nexus_transcript_fixture_path)
    canonical = gns.parse_canonical_transcript(payload)

    # Independently derive CDS genomic bounds from the fixture's own UTR
    # entries (not copy-pasted from the implementation): the 5' UTR sits at
    # the high-coordinate end of the first exon and the 3' UTR at the
    # low-coordinate end of the last exon (minus-strand gene).
    five_prime_utr = next(u for u in payload["utrs"] if u["type"] == "five_prime_UTR")
    three_prime_utr = next(u for u in payload["utrs"] if u["type"] == "three_prime_UTR")
    cds_max_genomic = five_prime_utr["start"] - 1
    cds_min_genomic = three_prime_utr["end"] + 1

    # Pick a breakpoint inside an interior (fully-coding, no-UTR) exon.
    exon5 = next(e for e in canonical.exons if e.rank == 5)
    breakpoint_genomic = exon5.start + 62  # well inside the exon

    # Hand-compute the expected protein position independently, directly
    # from the fixture's raw exon list.
    ordered = sorted(payload["exons"], key=lambda e: e["rank"])

    def coding_len(exon):
        start = max(exon["exonStart"], cds_min_genomic)
        end = min(exon["exonEnd"], cds_max_genomic)
        return max(0, end - start + 1)

    preceding = sum(coding_len(e) for e in ordered if e["rank"] < 5)
    offset = exon5.end - breakpoint_genomic  # minus strand: counts down from exon end
    expected_cds_nt = preceding + offset + 1
    expected_protein_position = math.ceil(expected_cds_nt / 3)

    result = gns.map_genomic_breakpoint_to_protein_position(
        canonical.exons,
        breakpoint_genomic,
        strand=canonical.exons[0].strand,
        cds_min_genomic=cds_min_genomic,
        cds_max_genomic=cds_max_genomic,
    )

    assert result.exon_rank == 5
    assert result.is_intronic is False
    assert result.protein_position == expected_protein_position


def test_cds_offset_handles_intronic_breakpoint_by_clamping_to_nearest_exon(
    genome_nexus_transcript_fixture_path,
):
    payload = _load(genome_nexus_transcript_fixture_path)
    canonical = gns.parse_canonical_transcript(payload)

    exon5 = next(e for e in canonical.exons if e.rank == 5)
    exon6 = next(e for e in canonical.exons if e.rank == 6)
    # A position strictly between exon 5 and exon 6 (intronic on this
    # minus-strand transcript, since rank 6 has lower genomic coordinates
    # than rank 5).
    intronic_position = (exon5.start + exon6.end) // 2
    assert exon6.end < intronic_position < exon5.start

    result = gns.map_genomic_breakpoint_to_protein_position(
        canonical.exons, intronic_position, strand=canonical.exons[0].strand
    )

    assert result.is_intronic is True
    assert result.exon_rank in (5, 6)


def _synthetic_plus_strand_exons():
    # Two coding exons on a plus-strand transcript with an intron between
    # them: exon 1 spans 100-200, exon 2 spans 203-400. Gene-agnostic,
    # synthetic coordinates -- not tied to any real gene.
    return [
        gns.ExonRecord(exon_id="synthetic-exon-1", start=100, end=200, rank=1, strand=1),
        gns.ExonRecord(exon_id="synthetic-exon-2", start=203, end=400, rank=2, strand=1),
    ]


def test_intronic_breakpoint_role_overrides_geometrically_nearer_wrong_side_exon():
    """An intronic breakpoint must never snap to whichever exon is merely
    closer in raw base pairs when the target gene's 5'/3' role in the
    fusion is known -- it must snap to the correct side of the splice
    junction for that role, even when that side is geometrically farther.
    """
    exons = _synthetic_plus_strand_exons()

    # Breakpoint at 202: nearest by raw genomic distance is exon 2's start
    # (203, distance 1) over exon 1's end (200, distance 2). A 5'-partner
    # gene retains everything through the end of exon 1 and must resolve
    # there regardless -- never into the downstream exon just because it
    # is one base pair closer.
    five_prime_result = gns.map_genomic_breakpoint_to_protein_position(
        exons, breakpoint_genomic=202, strand=1, role="five_prime"
    )
    assert five_prime_result.is_intronic is True
    assert five_prime_result.exon_rank == 1

    # Confirm the naive nearest-distance fallback (no role supplied) really
    # would have picked the other, wrong-for-five_prime exon here -- this
    # is what makes the case above a genuine semantic gap, not a no-op.
    no_role_result = gns.map_genomic_breakpoint_to_protein_position(
        exons, breakpoint_genomic=202, strand=1
    )
    assert no_role_result.exon_rank == 2

    # Mirror case, breakpoint at 201: nearest by raw distance is now exon
    # 1's end (200, distance 1) over exon 2's start (203, distance 2). A
    # 3'-partner gene retains everything from the start of exon 2 onward
    # and must resolve there regardless -- never into the upstream exon.
    three_prime_result = gns.map_genomic_breakpoint_to_protein_position(
        exons, breakpoint_genomic=201, strand=1, role="three_prime"
    )
    assert three_prime_result.is_intronic is True
    assert three_prime_result.exon_rank == 2

    no_role_mirror_result = gns.map_genomic_breakpoint_to_protein_position(
        exons, breakpoint_genomic=201, strand=1
    )
    assert no_role_mirror_result.exon_rank == 1


def _synthetic_minus_strand_exons():
    # Two coding exons on a minus-strand transcript: transcript order runs
    # from high to low genomic coordinate, so rank 1 (the 5' start) is the
    # *higher*-coordinate exon (203-400) and rank 2 is the lower-coordinate
    # exon (100-200), with a 2bp intronic gap between them (201-202).
    # Gene-agnostic, synthetic coordinates -- not tied to any real gene.
    return [
        gns.ExonRecord(exon_id="synthetic-exon-1", start=203, end=400, rank=1, strand=-1),
        gns.ExonRecord(exon_id="synthetic-exon-2", start=100, end=200, rank=2, strand=-1),
    ]


def test_intronic_breakpoint_role_overrides_nearer_exon_on_minus_strand():
    """Same semantic gap as the plus-strand case above, but on a minus
    strand transcript, where transcript order runs opposite to increasing
    genomic coordinate -- the upstream/downstream side of an intronic
    breakpoint flips accordingly, and the fix must still get it right.
    """
    exons = _synthetic_minus_strand_exons()

    # Breakpoint at 201: nearest by raw genomic distance is exon 2's end
    # (200, distance 1) over exon 1's start (203, distance 2). But exon 2
    # is the *downstream* (3'-side) exon on this minus-strand transcript, so
    # a 5'-partner gene must still resolve to exon 1 (upstream), not the
    # geometrically nearer exon 2.
    five_prime_result = gns.map_genomic_breakpoint_to_protein_position(
        exons, breakpoint_genomic=201, strand=-1, role="five_prime"
    )
    assert five_prime_result.is_intronic is True
    assert five_prime_result.exon_rank == 1

    no_role_result = gns.map_genomic_breakpoint_to_protein_position(
        exons, breakpoint_genomic=201, strand=-1
    )
    assert no_role_result.exon_rank == 2

    # Mirror case, breakpoint at 202: nearest by raw distance is now exon
    # 1's start (203, distance 1) over exon 2's end (200, distance 2), but
    # exon 1 is the *upstream* (5'-side) exon here, so a 3'-partner gene
    # must resolve to exon 2 instead.
    three_prime_result = gns.map_genomic_breakpoint_to_protein_position(
        exons, breakpoint_genomic=202, strand=-1, role="three_prime"
    )
    assert three_prime_result.is_intronic is True
    assert three_prime_result.exon_rank == 2

    no_role_mirror_result = gns.map_genomic_breakpoint_to_protein_position(
        exons, breakpoint_genomic=202, strand=-1
    )
    assert no_role_mirror_result.exon_rank == 1


def test_resolve_domains_falls_back_to_uniprot_when_gene_not_in_genome_nexus(
    uniprot_fixture_path,
):
    mock_gn_client = MagicMock()
    mock_gn_client.fetch_canonical_transcript.side_effect = gns.GenomeNexusGeneNotFound("not found")

    uniprot_source = UniProtDomainSource()
    uniprot_payload = json.loads(uniprot_fixture_path.read_text())
    uniprot_source._cache["P15056"] = uniprot_source.parse(uniprot_payload)

    domains = gns.resolve_domains(
        "SOME_GENE_NOT_IN_GENOME_NEXUS",
        "P15056",
        genome_nexus_client=mock_gn_client,
        uniprot_source=uniprot_source,
    )

    assert any("kinase" in d.name.lower() for d in domains)
    mock_gn_client.fetch_canonical_transcript.assert_called_once()


def test_resolve_domains_uses_genome_nexus_when_available(
    genome_nexus_canonical_transcript_fixture_path,
):
    mock_gn_client = MagicMock()
    mock_gn_client.fetch_canonical_transcript.return_value = _load(
        genome_nexus_canonical_transcript_fixture_path
    )
    uniprot_source = MagicMock()

    domains = gns.resolve_domains(
        "SOME_GENE",
        "SOME_ACCESSION",
        genome_nexus_client=mock_gn_client,
        uniprot_source=uniprot_source,
    )

    assert any(d.name == "PF07714" for d in domains)
    uniprot_source.fetch.assert_not_called()


def test_fetch_canonical_transcript_raises_gene_not_found_on_404():
    mock_session = MagicMock()
    mock_session.get.return_value.status_code = 404
    client = gns.GenomeNexusClient(session=mock_session)

    with pytest.raises(gns.GenomeNexusGeneNotFound):
        client.fetch_canonical_transcript("NOT_A_REAL_GENE")


def test_fetch_canonical_transcript_caches_after_first_call(
    genome_nexus_canonical_transcript_fixture_path,
):
    mock_session = MagicMock()
    mock_session.get.return_value.status_code = 200
    mock_session.get.return_value.json.return_value = _load(
        genome_nexus_canonical_transcript_fixture_path
    )
    client = gns.GenomeNexusClient(session=mock_session)

    first = client.fetch_canonical_transcript("SOME_GENE")
    second = client.fetch_canonical_transcript("SOME_GENE")

    assert first == second
    mock_session.get.assert_called_once()


def test_get_retries_on_5xx_then_succeeds(monkeypatch):
    mock_session = MagicMock()
    failing_response = MagicMock(status_code=503)
    ok_response = MagicMock(status_code=200)
    mock_session.get.side_effect = [failing_response, ok_response]

    monkeypatch.setattr(gns.time, "sleep", lambda *_: None)
    client = gns.GenomeNexusClient(session=mock_session, max_retries=3, backoff_seconds=0.01)

    result = client._get("/ensembl/canonical-transcript/hgnc/SOME_GENE")

    assert result is ok_response
    assert mock_session.get.call_count == 2


def test_get_gives_up_after_max_retries(monkeypatch):
    mock_session = MagicMock()
    always_failing = MagicMock(status_code=500)
    mock_session.get.return_value = always_failing

    monkeypatch.setattr(gns.time, "sleep", lambda *_: None)
    client = gns.GenomeNexusClient(session=mock_session, max_retries=2, backoff_seconds=0.01)

    result = client._get("/ensembl/canonical-transcript/hgnc/SOME_GENE")

    assert result is always_failing
    assert mock_session.get.call_count == 3  # initial attempt + 2 retries


@pytest.mark.network
def test_fetch_canonical_transcript_real_network_call():
    """Excluded from default `pytest -m "not network"` runs."""
    client = gns.GenomeNexusClient()
    payload = client.fetch_canonical_transcript("BRAF")
    assert payload["refseqMrnaId"] == "NM_004333"
