"""HGNC-verified single-gene fallback for batch canonical-transcript misses."""

from __future__ import annotations

from unittest.mock import MagicMock

import requests

from cfh.cohort.transcript_fallback import fetch_hgnc_identity, resolve_with_fallback
from cfh.mapping.genome_nexus_source import GenomeNexusGeneNotFound


def _hgnc_session(doc: dict | None):
    session = MagicMock()
    session.get.return_value.json.return_value = {"response": {"docs": [doc] if doc else []}}
    return session


def _gn_client(payload=None, error=None):
    client = MagicMock()
    if error is not None:
        client.fetch_canonical_transcript.side_effect = error
    else:
        client.fetch_canonical_transcript.return_value = payload
    return client


def _payload(gene_id: str, protein_id: str | None = "ENSP1", hugo=("ALIAS",)):
    return {
        "transcriptId": "ENST1",
        "geneId": gene_id,
        "proteinId": protein_id,
        "proteinLength": 100,
        "pfamDomains": [],
        "exons": [],
        "utrs": [],
        "hugoSymbols": list(hugo),
    }


TCF3 = {
    "symbol": "TCF3",
    "ensembl_gene_id": "ENSG00000071564",
    "locus_group": "protein-coding gene",
}


def test_mislabelled_payload_is_accepted_when_gene_id_matches():
    # Observed live: TCF3's Genome Nexus payload lists hugoSymbols ["TCF7L1"].
    resolution = resolve_with_fallback(
        "TCF3",
        6929,
        genome_nexus_client=_gn_client(_payload("ENSG00000071564", hugo=("TCF7L1",))),
        session=_hgnc_session(TCF3),
    )
    assert resolution.status == "resolved"
    assert resolution.canonical.protein_id == "ENSP1"


def test_payload_for_a_different_gene_is_rejected():
    resolution = resolve_with_fallback(
        "TCF3",
        6929,
        genome_nexus_client=_gn_client(_payload("ENSG_OTHER")),
        session=_hgnc_session(TCF3),
    )
    assert resolution.status == "unresolved"
    assert "ENSG_OTHER" in resolution.reason


def test_non_coding_and_missing_protein_are_distinguished():
    non_coding = resolve_with_fallback(
        "LINC00114",
        None,
        genome_nexus_client=_gn_client(_payload("ENSG_X")),
        session=_hgnc_session(
            {"symbol": "LINC00114", "ensembl_gene_id": "ENSG_X", "locus_group": "non-coding RNA"}
        ),
    )
    assert non_coding.status == "non_coding"
    no_protein = resolve_with_fallback(
        "TCF3",
        6929,
        genome_nexus_client=_gn_client(_payload("ENSG00000071564", protein_id=None)),
        session=_hgnc_session(TCF3),
    )
    assert no_protein.status == "unresolved"
    assert "no protein" in no_protein.reason


def test_data_source_failures_degrade_to_unresolved():
    failing = MagicMock()
    failing.get.side_effect = requests.ConnectionError("down")
    assert (
        resolve_with_fallback(
            "TCF3", 6929, genome_nexus_client=_gn_client(), session=failing
        ).status
        == "unresolved"
    )
    assert (
        resolve_with_fallback(
            "TCF3", 6929, genome_nexus_client=_gn_client(), session=_hgnc_session(None)
        ).reason
        == "HGNC has no record for this gene"
    )
    missing = resolve_with_fallback(
        "TCF3",
        6929,
        genome_nexus_client=_gn_client(error=GenomeNexusGeneNotFound("x")),
        session=_hgnc_session(TCF3),
    )
    assert missing.status == "unresolved"


def test_hgnc_lookup_prefers_entrez_and_caches(tmp_path):
    session = _hgnc_session(TCF3)
    first = fetch_hgnc_identity("TCF3", 6929, session=session, cache_dir=tmp_path)
    second = fetch_hgnc_identity("TCF3", 6929, session=session, cache_dir=tmp_path)
    assert first == second
    assert session.get.call_count == 1
    assert session.get.call_args.args[0].endswith("/fetch/entrez_id/6929")
    fetch_hgnc_identity("TCF3", None, session=session, cache_dir=tmp_path)
    assert session.get.call_args.args[0].endswith("/fetch/symbol/TCF3")
