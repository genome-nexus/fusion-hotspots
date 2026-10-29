import json
from pathlib import Path
from unittest.mock import MagicMock

import pytest
import requests

from cfh.mapping.domain_reconstitution import (
    COMPLEMENTARY,
    GAP,
    OVERLAPPING,
    PARTNER_DOMAIN_INTACT,
    PARTNER_DOMAIN_LOST,
    PARTNER_LACKS_DOMAIN,
    PARTNER_UNMAPPED,
    TARGET_DOMAIN_NOT_IN_PFAM,
    assess_reconstitution,
    kept_span,
    pfam_hit_for_interval,
    reconstitution_columns,
)
from cfh.mapping.genome_nexus_source import (
    GenomeNexusClient,
    GenomeNexusGeneNotFound,
    PfamDomain,
    parse_canonical_transcript,
)
from cfh.real_benchmark import _domain_reconstitution

FIXTURES = Path(__file__).parent / "fixtures" / "genome_nexus"
ETS = "PF00178"
# Two 100-residue copies of one family at different absolute positions.
TARGET = PfamDomain(ETS, 201, 300)
PARTNER = PfamDomain(ETS, 501, 600)


def test_kept_span_matches_retention_boundaries():
    # The breakpoint residue stays with the surviving fragment.
    assert kept_span(TARGET, 250, "five_prime") == (0.0, 0.5)
    assert kept_span(TARGET, 250, "three_prime") == (0.49, 1.0)
    assert kept_span(TARGET, 400, "five_prime") == (0.0, 1.0)
    assert kept_span(TARGET, 100, "five_prime") is None
    assert kept_span(TARGET, 301, "three_prime") is None


@pytest.mark.parametrize(
    ("partner_breakpoint", "status", "gap_residues"),
    [
        (551, COMPLEMENTARY, 0.0),  # 5' keeps 0-0.50, 3' resumes at 0.50
        (545, COMPLEMENTARY, -6.0),  # small overlap within tolerance
        (531, OVERLAPPING, -20.0),
        (571, GAP, 20.0),
        (450, PARTNER_DOMAIN_INTACT, -50.0),
    ],
)
def test_assess_reconstitution_classifies_the_junction(partner_breakpoint, status, gap_residues):
    result = assess_reconstitution(
        target_domain=TARGET,
        target_breakpoint_aa=250,
        target_role="five_prime",
        partner_pfam_domains=[PfamDomain("PF00001", 10, 90), PARTNER],
        partner_breakpoint_aa=partner_breakpoint,
    )
    assert result.status == status
    assert result.junction_gap_residues == gap_residues
    assert result.five_prime_span == (0.0, 0.5)
    assert (result.partner_domain_start_aa, result.partner_domain_end_aa) == (501, 600)


def test_assess_reconstitution_is_symmetric_in_role():
    """The 3' gene's view of a junction reports the same spans and gap."""
    as_five = assess_reconstitution(
        target_domain=TARGET,
        target_breakpoint_aa=250,
        target_role="five_prime",
        partner_pfam_domains=[PARTNER],
        partner_breakpoint_aa=551,
    )
    as_three = assess_reconstitution(
        target_domain=PARTNER,
        target_breakpoint_aa=551,
        target_role="three_prime",
        partner_pfam_domains=[TARGET],
        partner_breakpoint_aa=250,
    )
    assert as_five.status == as_three.status == COMPLEMENTARY
    assert as_five.five_prime_span == as_three.five_prime_span
    assert as_five.three_prime_span == as_three.three_prime_span
    assert as_five.junction_gap_fraction == as_three.junction_gap_fraction


def test_partner_without_the_family_or_with_it_cut_away():
    lacks = assess_reconstitution(
        target_domain=TARGET,
        target_breakpoint_aa=250,
        target_role="five_prime",
        partner_pfam_domains=[PfamDomain("PF00001", 10, 90)],
        partner_breakpoint_aa=551,
    )
    assert lacks.status == PARTNER_LACKS_DOMAIN
    assert lacks.five_prime_span is None
    # The 3' partner keeps only residues after its domain.
    lost = assess_reconstitution(
        target_domain=TARGET,
        target_breakpoint_aa=250,
        target_role="five_prime",
        partner_pfam_domains=[PARTNER],
        partner_breakpoint_aa=650,
    )
    assert lost.status == PARTNER_DOMAIN_LOST


def test_repeated_partner_copies_judge_the_closest_complement():
    result = assess_reconstitution(
        target_domain=TARGET,
        target_breakpoint_aa=250,
        target_role="five_prime",
        partner_pfam_domains=[PfamDomain(ETS, 401, 500), PARTNER],
        partner_breakpoint_aa=551,
    )
    assert result.status == COMPLEMENTARY
    assert result.partner_domain_start_aa == 501


def test_pfam_hit_for_interval_picks_the_largest_overlap():
    hits = [PfamDomain("PF1", 100, 160), PfamDomain("PF2", 150, 260)]
    assert pfam_hit_for_interval(hits, 140, 250).pfam_id == "PF2"
    assert pfam_hit_for_interval(hits, 300, 400) is None
    assert pfam_hit_for_interval(hits, None, 400) is None


def _client(*genes):
    payloads = {
        gene: json.loads((FIXTURES / f"canonical_transcript_{gene.lower()}.json").read_text())
        for gene in genes
    }
    client = MagicMock(spec=GenomeNexusClient)

    def fetch(symbol):
        if symbol not in payloads:
            raise GenomeNexusGeneNotFound(symbol)
        return payloads[symbol]

    client.fetch_canonical_transcript.side_effect = fetch
    return client, {gene: parse_canonical_transcript(payload) for gene, payload in payloads.items()}


def test_fli1_erg_reconstitutes_the_ets_domain_from_either_side():
    """P-0053737-T01-IM6, MSK-IMPACT 50k. The source row labels site1 FLI1
    but gives it ERG's chr21 coordinate, so the partner breakpoint must be
    chosen by locus, as the target's is."""
    row = {
        "Site1_Hugo_Symbol": "FLI1",
        "Site1_Position": 39755655,
        "Site2_Hugo_Symbol": "ERG",
        "Site2_Position": 128680564,
    }
    client, canonical = _client("ERG", "FLI1")
    erg_domain = next(d for d in canonical["ERG"].pfam_domains if d.pfam_id == ETS)
    fli1_domain = next(d for d in canonical["FLI1"].pfam_domains if d.pfam_id == ETS)

    from_erg = _domain_reconstitution(
        row, "FLI1", "three_prime", 370, erg_domain.start_aa, erg_domain.end_aa,
        canonical["ERG"], client,
    )  # fmt: skip
    from_fli1 = _domain_reconstitution(
        row, "ERG", "five_prime", 347, fli1_domain.start_aa, fli1_domain.end_aa,
        canonical["FLI1"], client,
    )  # fmt: skip

    for result in (from_erg, from_fli1):
        assert result.status == COMPLEMENTARY
        assert result.pfam_accession == ETS
        assert result.five_prime_span == pytest.approx((0.0, 0.8193), abs=1e-4)
        assert result.three_prime_span == pytest.approx((0.7229, 1.0), abs=1e-4)
        assert result.junction_gap_residues == -8.0
        assert result.partner_breakpoint_mid_exon is True
    assert from_erg.partner_breakpoint_protein_position == 347
    assert from_fli1.partner_breakpoint_protein_position == 370


def test_gabpa_erg_pieces_overlap_rather_than_complement():
    """P-0023818-T01-IM6: GABPA breaks in an intron after its exon 9, ERG
    mid-exon 10, so the two kept pieces share about 51 residues of ETS."""
    row = {
        "Site1_Hugo_Symbol": "GABPA",
        "Site1_Position": 27139889,
        "Site2_Hugo_Symbol": "ERG",
        "Site2_Position": 39755807,
    }
    client, canonical = _client("ERG", "GABPA")
    erg_domain = next(d for d in canonical["ERG"].pfam_domains if d.pfam_id == ETS)

    result = _domain_reconstitution(
        row, "GABPA", "three_prime", 320, erg_domain.start_aa, erg_domain.end_aa,
        canonical["ERG"], client,
    )  # fmt: skip

    assert result.status == OVERLAPPING
    assert result.junction_gap_residues == -51.0
    assert result.partner_breakpoint_protein_position == 379
    assert result.partner_breakpoint_mid_exon is False


def test_unresolvable_partner_or_target_domain_is_reported_not_raised():
    row = {"Site1_Position": 39755655, "Site2_Position": 128680564}
    client, canonical = _client("ERG")
    erg_domain = next(d for d in canonical["ERG"].pfam_domains if d.pfam_id == ETS)

    unknown_gene = _domain_reconstitution(
        row, "NOT_A_GENE", "three_prime", 370, erg_domain.start_aa, erg_domain.end_aa,
        canonical["ERG"], client,
    )  # fmt: skip
    assert unknown_gene.status == PARTNER_UNMAPPED

    no_pfam = _domain_reconstitution(
        row, "FLI1", "three_prime", 370, 1, 5, canonical["ERG"], client
    )
    assert no_pfam.status == TARGET_DOMAIN_NOT_IN_PFAM


def test_partner_lookup_network_failure_propagates():
    """A flaky lookup must fail the gene (and be retried), never quietly
    record a different answer."""
    _, canonical = _client("ERG")
    client = MagicMock(spec=GenomeNexusClient)
    client.fetch_canonical_transcript.side_effect = requests.ConnectionError("DNS failed")
    erg_domain = next(d for d in canonical["ERG"].pfam_domains if d.pfam_id == ETS)
    with pytest.raises(requests.ConnectionError):
        _domain_reconstitution(
            {"Site1_Position": 39755655, "Site2_Position": 128680564},
            "FLI1", "three_prime", 370, erg_domain.start_aa, erg_domain.end_aa,
            canonical["ERG"], client,
        )  # fmt: skip


def test_reconstitution_columns_keep_one_column_set():
    result = assess_reconstitution(
        target_domain=TARGET,
        target_breakpoint_aa=250,
        target_role="five_prime",
        partner_pfam_domains=[PARTNER],
        partner_breakpoint_aa=551,
    )
    filled = reconstitution_columns(result)
    blank = reconstitution_columns(None)
    assert filled.keys() == blank.keys()
    assert "reconstitution_status" not in filled
    assert filled["reconstitution_five_prime_span"] == "0.0000-0.5000"
    assert filled["reconstitution_three_prime_span"] == "0.5000-1.0000"
    assert filled["reconstitution_junction_gap_residues"] == 0.0
    assert set(blank.values()) == {None}


def test_report_section_lists_only_partners_sharing_the_family():
    from types import SimpleNamespace

    from cfh.real_benchmark import _domain_reconstitution_lines

    def row(fusion, partner, status, **extra):
        return {
            "fusion_name": fusion,
            "sample_id": "S1",
            "partner_gene": partner,
            "frame_status": None,
            "is_intronic_breakpoint": False,
            "domain_reconstitution": status,
            **reconstitution_columns(None),
            **extra,
        }

    run = SimpleNamespace(
        gene_symbol="ERG",
        rows=[
            row(
                "FLI1-ERG",
                "FLI1",
                COMPLEMENTARY,
                reconstitution_five_prime_span="0.0000-0.8193",
                reconstitution_three_prime_span="0.7229-1.0000",
                reconstitution_junction_gap_residues=-8.0,
                reconstitution_partner_breakpoint_mid_exon=True,
            ),
            row("X-ERG", "X", PARTNER_LACKS_DOMAIN),
            row("TMPRSS2-ERG", "TMPRSS2", None),
        ],
        summary={"domain_reconstitution_counts": {COMPLEMENTARY: 1, PARTNER_LACKS_DOMAIN: 1}},
    )
    text = "\n".join(_domain_reconstitution_lines(run))
    assert "2 fusion(s) cut ERG's key domain; in 1 the partner" in text
    assert (
        "| FLI1-ERG | S1 | unknown | 0.0000-0.8193 | 0.7229-1.0000 | -8.0 | ERG, FLI1 "
        "| partner keeps the rest of the domain |"
    ) in text
    assert "X-ERG" not in text

    run.rows = run.rows[1:]
    assert _domain_reconstitution_lines(run) == []
