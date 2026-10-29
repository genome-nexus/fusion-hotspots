from types import SimpleNamespace

import pytest

from cfh.cohort.outputs import (
    _domain_reconstitution_markdown,
    build_domain_reconstitution_rows,
    domain_reconstitution_overview,
    random_cut_completion_rate,
)
from cfh.mapping.domain_reconstitution import (
    COMPLEMENTARY,
    OVERLAPPING,
    PARTNER_LACKS_DOMAIN,
)


def _row(sample, patient, fusion, partner, role, status, **extra):
    return {
        "sample_id": sample,
        "patient_id": patient,
        "fusion_name": fusion,
        "partner_gene": partner,
        "target_role": role,
        "frame_status": None,
        "is_intronic_breakpoint": False,
        "domain_reconstitution": status,
        "reconstitution_pfam_accession": "PF00178" if status != PARTNER_LACKS_DOMAIN else None,
        **extra,
    }


def _outcome(gene, rows):
    return SimpleNamespace(gene_symbol=gene, run=SimpleNamespace(rows=rows))


def _result():
    fli1_erg = dict(
        reconstitution_five_prime_span="0.0000-0.8193",
        reconstitution_three_prime_span="0.7229-1.0000",
        reconstitution_junction_gap_residues=-8.0,
    )
    erg = _outcome(
        "ERG",
        [
            _row("S1", "P1", "FLI1-ERG", "FLI1", "three_prime", COMPLEMENTARY,
                 reconstitution_partner_breakpoint_mid_exon=True, **fli1_erg),
            _row("S2", "P2", "GABPA-ERG", "GABPA", "three_prime", OVERLAPPING,
                 reconstitution_partner_breakpoint_mid_exon=False),
            _row("S3", "P3", "TMPRSS2-ERG", "TMPRSS2", "three_prime", None),
        ],
    )  # fmt: skip
    fli1 = _outcome(
        "FLI1",
        [
            _row("S1", "P1", "FLI1-ERG", "ERG", "five_prime", COMPLEMENTARY,
                 reconstitution_partner_breakpoint_mid_exon=True, **fli1_erg),
        ],
    )  # fmt: skip
    alk = _outcome(
        "EML4",
        [_row("S4", "P4", "EML4-ALK", "ALK", "five_prime", PARTNER_LACKS_DOMAIN)],
    )
    failed = SimpleNamespace(gene_symbol="BROKEN", run=None)
    return SimpleNamespace(gene_outcomes=[erg, fli1, alk, failed])


def test_a_fusion_seen_from_both_genes_counts_once():
    rows = build_domain_reconstitution_rows(_result())

    assert [row["fusion_name"] for row in rows] == ["FLI1-ERG", "GABPA-ERG", "EML4-ALK"]
    fli1_erg = rows[0]
    assert fli1_erg["reported_by"] == "ERG, FLI1"
    assert (fli1_erg["five_prime_gene"], fli1_erg["three_prime_gene"]) == ("FLI1", "ERG")
    assert fli1_erg["five_prime_mid_exon"] is True and fli1_erg["three_prime_mid_exon"] is True
    gabpa_erg = rows[1]
    assert gabpa_erg["five_prime_mid_exon"] is False  # GABPA, the partner
    assert gabpa_erg["three_prime_mid_exon"] is True  # ERG, the target


def test_overview_reports_completion_frequency():
    rows = build_domain_reconstitution_rows(_result())
    overview = domain_reconstitution_overview(rows)

    assert overview["truncating_fusions"] == 3
    assert overview["truncating_patients"] == 3
    assert overview["partner_shares_family"] == 2
    assert overview["both_copies_cut"] == 2
    assert overview["complementary"] == 1
    assert overview["completion_rate_among_truncations"] == pytest.approx(1 / 3)
    assert overview["completion_rate_when_both_cut"] == 0.5
    assert overview["by_pfam_family"] == {
        "PF00178": {"partner_shares_family": 2, COMPLEMENTARY: 1, OVERLAPPING: 1}
    }
    assert overview["random_cut_completion_rate"] == pytest.approx(0.19)

    text = "\n".join(_domain_reconstitution_markdown(overview, rows))
    assert (
        "- Partner completes the domain: 1 (33% of truncating fusions; 50% of those "
        "cutting both copies, against 19% if both cuts fell at random)"
    ) in text
    assert "| FLI1-ERG | S1 | PF00178 | unknown | 0.0000-0.8193 | 0.7229-1.0000 | -8.0 |" in text
    assert "EML4-ALK" not in text


def test_empty_cohort_and_random_cut_reference():
    overview = domain_reconstitution_overview([])
    assert overview["completion_rate_among_truncations"] is None
    assert overview["completion_rate_when_both_cut"] is None
    assert "No fusion cut a scanned gene's key domain." in _domain_reconstitution_markdown(
        overview, []
    )
    assert random_cut_completion_rate(0.0) == 0.0
    assert random_cut_completion_rate(1.0) == 1.0
