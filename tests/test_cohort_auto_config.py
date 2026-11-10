"""Tests for auto-generated GeneConfig construction (mocked Genome Nexus
batch response) and the kinase/catalytic key-domain heuristic."""

from __future__ import annotations

from unittest.mock import MagicMock

import pytest
import requests

from cfh.cohort.auto_config import (
    PfamDescriptionSource,
    batch_fetch_canonical_transcripts,
    build_auto_gene_config,
    select_key_domain,
)
from cfh.mapping.genome_nexus_source import CanonicalTranscript, ExonRecord, PfamDomain


def _canonical(pfam_domains, *, protein_id="ENSP00000TEST", transcript_id="ENST00000TEST"):
    return CanonicalTranscript(
        transcript_id=transcript_id,
        refseq_mrna_id="NM_000000",
        protein_id=protein_id,
        protein_length=500,
        uniprot_id=None,
        pfam_domains=pfam_domains,
        exons=[],
        utrs=[],
    )


def _mock_description_source(descriptions: dict[str, str]) -> PfamDescriptionSource:
    source = PfamDescriptionSource.__new__(PfamDescriptionSource)
    source._cache = dict(descriptions)
    source.cache_dir = None
    return source


def test_select_key_domain_prefers_kinase_description():
    domains = [
        PfamDomain(pfam_id="PF00001", start_aa=10, end_aa=400),  # largest, non-kinase
        PfamDomain(pfam_id="PF07714", start_aa=458, end_aa=712),  # smaller, kinase
    ]
    description_source = _mock_description_source(
        {"PF00001": "Some large uninteresting domain", "PF07714": "Protein tyrosine kinase"}
    )

    chosen = select_key_domain(domains, description_source)

    assert chosen is not None
    assert chosen.accession == "PF07714"
    assert chosen.key == "kinase"
    assert "kinase" in chosen.name.lower()


def test_select_key_domain_matches_catalytic_keyword_too():
    domains = [PfamDomain(pfam_id="PF00069", start_aa=1, end_aa=100)]
    description_source = _mock_description_source({"PF00069": "Catalytic core domain"})

    chosen = select_key_domain(domains, description_source)

    assert chosen.key == "kinase"
    assert chosen.accession == "PF00069"


def test_select_key_domain_falls_back_to_largest_domain_when_no_kinase_match():
    domains = [
        PfamDomain(pfam_id="PF00001", start_aa=10, end_aa=400),  # span 390, largest
        PfamDomain(pfam_id="PF00002", start_aa=410, end_aa=420),  # span 10
    ]
    description_source = _mock_description_source(
        {"PF00001": "Some domain", "PF00002": "Another domain"}
    )

    chosen = select_key_domain(domains, description_source)

    assert chosen.accession == "PF00001"
    assert chosen.key == "auto_key_domain"


def test_select_key_domain_falls_back_to_largest_domain_without_description_source():
    domains = [
        PfamDomain(pfam_id="PF00001", start_aa=10, end_aa=50),
        PfamDomain(pfam_id="PF00002", start_aa=100, end_aa=500),
    ]

    chosen = select_key_domain(domains, description_source=None)

    assert chosen.accession == "PF00002"
    assert chosen.key == "auto_key_domain"


def test_select_key_domain_returns_none_for_no_domains():
    assert select_key_domain([], description_source=None) is None


def test_pfam_description_source_gracefully_degrades_on_network_failure():
    mock_session = MagicMock()
    mock_session.get.side_effect = requests.exceptions.ConnectionError("boom")
    source = PfamDescriptionSource(session=mock_session, max_retries=0)

    assert source.describe("PF07714") is None


def test_pfam_description_source_caches_across_calls():
    mock_session = MagicMock()
    response = MagicMock(status_code=200)
    response.json.return_value = {"metadata": {"name": {"name": "Protein kinase domain"}}}
    mock_session.get.return_value = response
    source = PfamDescriptionSource(session=mock_session)

    first = source.describe("PF00069")
    second = source.describe("PF00069")

    assert first == second == "Protein kinase domain"
    mock_session.get.assert_called_once()


def test_pfam_description_source_caches_on_disk(tmp_path):
    mock_session = MagicMock()
    response = MagicMock(status_code=200)
    response.json.return_value = {"metadata": {"name": {"name": "Protein kinase domain"}}}
    mock_session.get.return_value = response
    source = PfamDescriptionSource(session=mock_session, cache_dir=tmp_path)

    assert source.describe("PF00069") == "Protein kinase domain"
    assert (tmp_path / "PF00069.json").exists()

    second_source = PfamDescriptionSource(session=MagicMock(), cache_dir=tmp_path)
    assert second_source.describe("PF00069") == "Protein kinase domain"
    second_source.session.get.assert_not_called()


def test_build_auto_gene_config_populates_minimal_fields_only():
    canonical = _canonical(
        [
            PfamDomain(pfam_id="PF00001", start_aa=10, end_aa=40),
            PfamDomain(pfam_id="PF07714", start_aa=60, end_aa=90),
        ],
        protein_id="ENSP00000288602",
        transcript_id="ENST00000288602",
    )
    canonical.exons = [
        ExonRecord("e1", 1, 150, 1, 1),
        ExonRecord("e2", 151, 300, 2, 1),
    ]
    description_source = _mock_description_source(
        {"PF00001": "Regulatory domain", "PF07714": "Protein tyrosine kinase"}
    )

    config = build_auto_gene_config("braf", 673, canonical, description_source=description_source)

    assert config is not None
    assert config.gene_symbol == "BRAF"
    assert config.canonical_transcript_id == "NM_000000"
    assert config.protein_id == "ENSP00000288602"
    assert config.entrez_gene_id == 673
    assert len(config.key_domains) == 1
    assert config.key_domains[0].accession == "PF07714"
    assert [domain.name for domain in config.disruption_required_domains] == ["Regulatory domain"]
    assert config.expected_retained_exon_hint == "2"
    assert config.gene_pair is None


def test_build_auto_gene_config_handles_gene_with_no_pfam_domains():
    canonical = _canonical([], protein_id="ENSP00000TEST", transcript_id="ENST00000TEST")

    config = build_auto_gene_config("SOMEGENE", 42, canonical, description_source=None)

    assert config is not None
    assert config.key_domains == []


def test_build_auto_gene_config_returns_none_without_protein_id():
    canonical = _canonical([], protein_id=None, transcript_id="ENST00000TEST")
    assert build_auto_gene_config("SOMEGENE", 42, canonical) is None


@pytest.mark.network
def test_build_auto_gene_config_derives_defaults_from_live_genome_nexus():
    batch = batch_fetch_canonical_transcripts(["RET"])
    canonical = batch.by_gene_symbol["RET"]

    config = build_auto_gene_config("RET", 5979, canonical)

    assert config is not None
    assert config.expected_retained_exon_hint is not None
    assert all(
        domain.accession != config.key_domains[0].accession
        for domain in config.disruption_required_domains
    )


def test_batch_fetch_canonical_transcripts_posts_gene_list_in_one_call():
    mock_session = MagicMock()
    response = MagicMock(status_code=200)
    response.json.return_value = [
        {
            "transcriptId": "ENST00000288602",
            "refseqMrnaId": "NM_004333",
            "proteinId": "ENSP00000288602",
            "proteinLength": 766,
            "hugoSymbols": ["BRAF"],
            "pfamDomains": [
                {"pfamDomainId": "PF07714", "pfamDomainStart": 458, "pfamDomainEnd": 712}
            ],
            "exons": [],
            "utrs": [],
        },
        {
            "transcriptId": "ENST00000340058",
            "refseqMrnaId": "NM_020975",
            "proteinId": "ENSP00000340058",
            "proteinLength": 1114,
            "hugoSymbols": ["RET"],
            "pfamDomains": [],
            "exons": [],
            "utrs": [],
        },
    ]
    mock_session.post.return_value = response

    result = batch_fetch_canonical_transcripts(["BRAF", "RET", "NOTFOUND"], session=mock_session)

    assert mock_session.post.call_count == 1
    called_url = mock_session.post.call_args.args[0]
    _, kwargs = mock_session.post.call_args
    assert called_url.endswith("/ensembl/canonical-transcript/hgnc")
    assert kwargs["json"] == ["BRAF", "RET", "NOTFOUND"]

    assert set(result.by_gene_symbol) == {"BRAF", "RET"}
    assert result.by_gene_symbol["BRAF"].protein_id == "ENSP00000288602"
    assert result.unresolved_gene_symbols == ["NOTFOUND"]


def test_batch_fetch_canonical_transcripts_chunks_beyond_batch_size():
    mock_session = MagicMock()
    response = MagicMock(status_code=200)
    response.json.return_value = []
    mock_session.post.return_value = response

    batch_fetch_canonical_transcripts(
        [f"GENE{i}" for i in range(5)], session=mock_session, batch_size=2
    )

    assert mock_session.post.call_count == 3  # chunks of 2, 2, 1


def test_batch_fetch_canonical_transcripts_uses_disk_cache(tmp_path):
    mock_session = MagicMock()
    response = MagicMock(status_code=200)
    response.json.return_value = [
        {
            "transcriptId": "ENST00000288602",
            "refseqMrnaId": "NM_004333",
            "proteinId": "ENSP00000288602",
            "proteinLength": 766,
            "hugoSymbols": ["BRAF"],
            "pfamDomains": [],
            "exons": [],
            "utrs": [],
        }
    ]
    mock_session.post.return_value = response

    first = batch_fetch_canonical_transcripts(["BRAF"], session=mock_session, cache_dir=tmp_path)
    assert "BRAF" in first.by_gene_symbol
    assert mock_session.post.call_count == 1

    second_session = MagicMock()
    second = batch_fetch_canonical_transcripts(["BRAF"], session=second_session, cache_dir=tmp_path)
    assert "BRAF" in second.by_gene_symbol
    second_session.post.assert_not_called()


def test_batch_fetch_canonical_transcripts_retries_transient_request_errors():
    mock_session = MagicMock()
    response = MagicMock(status_code=200)
    response.json.return_value = []
    mock_session.post.side_effect = [requests.ReadTimeout("slow"), response]

    result = batch_fetch_canonical_transcripts(["GENE1"], session=mock_session, backoff_seconds=0)

    assert mock_session.post.call_count == 2
    assert result.unresolved_gene_symbols == ["GENE1"]


def test_batch_fetch_canonical_transcripts_reraises_after_max_retries():
    mock_session = MagicMock()
    mock_session.post.side_effect = requests.ConnectionError("down")

    with pytest.raises(requests.ConnectionError):
        batch_fetch_canonical_transcripts(
            ["GENE1"], session=mock_session, max_retries=1, backoff_seconds=0
        )
    assert mock_session.post.call_count == 2
