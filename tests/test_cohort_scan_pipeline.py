"""Offline/mocked multi-gene end-to-end cohort-scan test.

Runs :func:`cfh.cohort.scan.run_cohort_scan` against a single injected mock
``requests.Session`` (no real network) covering: cohort-wide recurrence
gating, curated-vs-auto-generated ``GeneConfig`` resolution, the REAL
8-algorithm orchestrator (:mod:`cfh.orchestrator.run`, not bypassed or
stubbed), and cross-gene Benjamini-Hochberg FDR correction across more than
two genes (BRAF, RET, and two auto-configured genes -- one with a kinase
Pfam domain, one with none at all).
"""

from __future__ import annotations

import json
from unittest.mock import MagicMock

import pytest
import requests

from cfh.cohort import scan as scan_module
from cfh.cohort.outputs import build_summary_rows, write_cohort_scan_outputs
from cfh.cohort.scan import genes_needing_full_report, run_cohort_scan
from cfh.ingestion import cbioportal_api
from cfh.mapping.genome_nexus_source import GenomeNexusClient
from cfh.studies.registry import StudyConfig

_STUDY_ID = "test_cohort_study"
_EXON_GENOMIC_START = 1_000

_GENE_SPECS = {
    # gene_symbol: (entrez_gene_id, distinct_patient_count, protein_id, pfam_domains)
    # Complete synthetic annotations keep the offline test out of UniProt fallback.
    # Curated YAML still determines which domains are tested.
    "BRAF": (
        673,
        20,
        "P15056_CURATED",
        [
            {"pfamDomainId": "PF07714", "pfamDomainStart": 100, "pfamDomainEnd": 300},
            {"pfamDomainId": "PF02196", "pfamDomainStart": 10, "pfamDomainEnd": 30},
            {"pfamDomainId": "PF00130", "pfamDomainStart": 40, "pfamDomainEnd": 50},
        ],
    ),
    "RET": (
        5979,
        19,
        "P07949_CURATED",
        [
            {"pfamDomainId": "PF07714", "pfamDomainStart": 100, "pfamDomainEnd": 300},
            {"pfamDomainId": "PF00028", "pfamDomainStart": 10, "pfamDomainEnd": 30},
        ],
    ),
    "FAKE1": (
        9001,
        10,
        "ENSP00009001",
        [{"pfamDomainId": "PF07714", "pfamDomainStart": 40, "pfamDomainEnd": 90}],
    ),
    "FAKE2": (9002, 6, "ENSP00009002", []),  # no Pfam domains at all
    "SINGLETON": (9003, 1, "ENSP00009003", []),  # filtered out by the recurrence gate
}


def _canonical_payload(gene_symbol: str, protein_id: str, pfam_domains: list[dict]) -> dict:
    return {
        "transcriptId": f"ENST_{gene_symbol}",
        "transcriptIdVersion": "1",
        "geneId": f"ENSG_{gene_symbol}",
        "refseqMrnaId": f"NM_{gene_symbol}",
        "proteinId": protein_id,
        "proteinLength": 500,
        "pfamDomains": pfam_domains,
        "exons": [
            {
                "exonId": f"ENSE_{gene_symbol}_1",
                "exonStart": _EXON_GENOMIC_START,
                "exonEnd": _EXON_GENOMIC_START + 3 * 1000 - 1,
                "rank": 1,
                "strand": 1,
            }
        ],
        "utrs": [],
        "uniprotId": None,
        "hugoSymbols": [gene_symbol],
    }


def _breakpoint_genomic_for_protein_position(position: int) -> int:
    """Inverse of the plus-strand CDS-offset arithmetic (no UTRs, single
    exon starting at ``_EXON_GENOMIC_START``): choose a genomic breakpoint
    that maps to exactly ``position``."""
    offset = 3 * (position - 1)
    return _EXON_GENOMIC_START + offset


def _sv_call(
    sample_id: str,
    partner: str,
    gene_symbol: str,
    *,
    position: int,
    in_frame: bool,
    tumor_variant_count: int = 20,
) -> dict:
    frame_text = "in frame" if in_frame else "out of frame"
    return {
        "sampleId": sample_id,
        "patientId": sample_id,
        "site1HugoSymbol": partner,
        "site2HugoSymbol": gene_symbol,
        "connectionType": "3to3",
        "site2Position": _breakpoint_genomic_for_protein_position(position),
        "site2EffectOnFrame": "NA",
        "tumorPairedEndReadCount": 10,
        "tumorSplitReadCount": 5,
        "tumorVariantCount": tumor_variant_count,
        "eventInfo": f"Protein Fusion: {frame_text}  {{{partner}:{gene_symbol}}}",
    }


def _sv_calls_for_gene(gene_symbol: str) -> list[dict]:
    """3 in-frame domain-retained-position events + 2 in-frame
    domain-lost-position events + 2 out-of-frame events, against two
    distinct partner genes so frequency/composite_score have something to
    rank."""
    calls = []
    for i in range(3):
        calls.append(
            _sv_call(
                f"{gene_symbol}-IN-{i}",
                "PARTNERA",
                gene_symbol,
                position=60,
                in_frame=True,
                tumor_variant_count=30 + i,
            )
        )
    for i in range(2):
        calls.append(
            _sv_call(
                f"{gene_symbol}-INLOST-{i}",
                "PARTNERB",
                gene_symbol,
                position=400,
                in_frame=True,
                tumor_variant_count=5 + i,
            )
        )
    for i in range(2):
        calls.append(
            _sv_call(
                f"{gene_symbol}-OUT-{i}",
                "PARTNERA",
                gene_symbol,
                position=400,
                in_frame=False,
                tumor_variant_count=5 + i,
            )
        )
    return calls


def _recurrence_records() -> list[dict]:
    return [
        {
            "hugoGeneSymbol": symbol,
            "entrezGeneId": entrez_id,
            "numberOfAlteredCases": patients,
            "totalCount": patients,
        }
        for symbol, (entrez_id, patients, _protein_id, _domains) in _GENE_SPECS.items()
    ]


@pytest.fixture
def mock_session() -> MagicMock:
    session = MagicMock()

    def _post(url, json=None, **kwargs):
        response = MagicMock(status_code=200)
        if url.endswith("/structuralvariant-genes/fetch"):
            response.json.return_value = _recurrence_records()
        elif url.endswith("/mutations/fetch"):
            response.json.return_value = []
        elif url.endswith("/clinical-data/fetch"):
            response.json.return_value = [
                {"sampleId": sample, "clinicalAttributeId": attribute, "value": value}
                for sample in json["ids"]
                for attribute, value in [("CANCER_TYPE", "Glioma"), ("ONCOTREE_CODE", "PA")]
            ]
        elif url.endswith("/structural-variant/fetch"):
            entrez_ids = json["entrezGeneIds"]
            gene_symbol = next(
                symbol
                for symbol, (entrez_id, *_rest) in _GENE_SPECS.items()
                if entrez_id in entrez_ids
            )
            response.json.return_value = _sv_calls_for_gene(gene_symbol)
        elif url.endswith("/ensembl/canonical-transcript/hgnc"):
            requested = json
            payloads = []
            for symbol in requested:
                _entrez_id, _patients, protein_id, domains = _GENE_SPECS[symbol]
                payloads.append(_canonical_payload(symbol, protein_id, domains or []))
            response.json.return_value = payloads
        else:
            raise AssertionError(f"unexpected POST url in cohort-scan pipeline test: {url}")
        return response

    def _get(url, params=None, **kwargs):
        response = MagicMock(status_code=200)
        if url.endswith("/sample-ids"):
            response.json.return_value = [
                call["sampleId"] for call in _sv_calls_for_gene("BRAF")
            ] + ["NEG-1", "NEG-2"]
            return response
        for symbol, (_entrez_id, _patients, protein_id, domains) in _GENE_SPECS.items():
            if url.endswith(f"/ensembl/canonical-transcript/hgnc/{symbol}"):
                response.json.return_value = _canonical_payload(symbol, protein_id, domains or [])
                return response
        if "/interpro/wwwapi/entry/pfam/" in url:
            accession = url.rsplit("/", 1)[-1]
            response.json.return_value = {
                "metadata": {"name": {"name": f"Description of {accession}"}}
            }
            return response
        raise AssertionError(f"unexpected GET url in cohort-scan pipeline test: {url}")

    session.post.side_effect = _post
    session.get.side_effect = _get
    return session


def test_cohort_scan_end_to_end_offline(mock_session, tmp_path):
    result = run_cohort_scan(
        _STUDY_ID,
        min_distinct_patients=5,
        n_permutations=200,
        adaptive=True,
        n_permutations_small=20,
        cache_dir=tmp_path / "cache",
        session=mock_session,
    )

    # Recurrence gating: the full pre-gate universe is reported, and the
    # singleton gene is genuinely filtered out, not silently disappeared.
    assert result.total_genes_before_gating == len(_GENE_SPECS)
    assert result.genes_after_gating == len(_GENE_SPECS) - 1
    scanned_symbols = {outcome.gene_symbol for outcome in result.gene_outcomes}
    assert scanned_symbols == {"BRAF", "RET", "FAKE1", "FAKE2"}
    assert "SINGLETON" not in scanned_symbols

    # Curated configs win for BRAF/RET; FAKE1/FAKE2 are auto-generated.
    outcomes_by_gene = {outcome.gene_symbol: outcome for outcome in result.gene_outcomes}
    assert outcomes_by_gene["BRAF"].config_source == "curated"
    assert outcomes_by_gene["RET"].config_source == "curated"
    assert outcomes_by_gene["FAKE1"].config_source == "auto"
    assert outcomes_by_gene["FAKE2"].config_source == "auto"
    assert result.curated_gene_count == 2
    assert result.auto_config_gene_count == 2

    # No gene crashes the whole scan; every gene produced a real orchestrator run.
    for outcome in result.gene_outcomes:
        assert outcome.status == "ok", outcome.error
        assert outcome.run is not None
        assert all(row["tumor_type"] == "Glioma" for row in outcome.run.rows)
        assert all(row["oncotree_code"] == "PA" for row in outcome.run.rows)
        algorithm_names = {r.Algorithm for r in outcome.run.results}
        assert "composite_score" in algorithm_names
        assert "confidence_stats" in algorithm_names
        assert "frequency" in algorithm_names

    braf_confidence = next(
        r for r in outcomes_by_gene["BRAF"].run.results if r.Algorithm == "confidence_stats"
    )
    assert braf_confidence.Parameters["group_field"] == "Domain_retention_flags"
    assert braf_confidence.Parameters["numeric_field"] == "Tumor_variant_count"
    assert "group_field' is required" not in " ".join(braf_confidence.Warnings)
    assert braf_confidence.Summary["mle"]["groups"]["retained"]["n"] > 0

    # FAKE2 (no Pfam domains) must gracefully no-op domain_retention rather
    # than crash the gene -- the same opt-in/no-op pattern proven elsewhere.
    fake2_domain_retention = next(
        r for r in outcomes_by_gene["FAKE2"].run.results if r.Algorithm == "domain_retention"
    )
    assert fake2_domain_retention.Summary["fisher_p_value"] is None
    assert fake2_domain_retention.Warnings

    # Cross-gene FDR correction spans every scanned gene that produced an
    # applicable p-value -- more than just two genes -- via the real
    # cfh.stats.multiple_testing.benjamini_hochberg, not a stub. FAKE2 has
    # no Pfam domain at all, so every domain-dependent algorithm gracefully
    # no-ops/fails for it and it legitimately contributes zero hypotheses.
    fdr_genes = {row["gene"] for row in result.fdr_rows}
    assert fdr_genes == {"BRAF", "RET", "FAKE1"}
    assert len(fdr_genes) > 2
    assert all(0.0 <= row["bh_adjusted_q"] <= 1.0 for row in result.fdr_rows)

    # Adaptive permutations actually ran (small budget requested).
    braf_domain_retention = next(
        r for r in outcomes_by_gene["BRAF"].run.results if r.Algorithm == "domain_retention"
    )
    assert braf_domain_retention.Summary["adaptive_permutations"]["enabled"] is True

    # Summary-building and output writing must not crash either.
    rows = build_summary_rows(result)
    assert {row["gene_symbol"] for row in rows} == {"BRAF", "RET", "FAKE1", "FAKE2"}
    braf_row = next(row for row in rows if row["gene_symbol"] == "BRAF")
    assert [domain["name"] for domain in braf_row["key_domains"]] == ["Protein kinase domain"]
    assert [domain["name"] for domain in braf_row["disruption_required_domains"]] == [
        "RAS-binding domain",
        "Cysteine-rich domain",
    ]
    # Sorted by significance: no row with a real q-value sorts after one with none.
    q_values = [row["min_fdr_adjusted_q_value"] for row in rows]
    real_positions = [i for i, q in enumerate(q_values) if q is not None]
    none_positions = [i for i, q in enumerate(q_values) if q is None]
    assert all(r < n for r in real_positions for n in none_positions)

    paths = write_cohort_scan_outputs(result, tmp_path / "runs", pdf=False)
    assert paths["summary_tsv"].exists()
    assert paths["summary_json"].exists()
    assert paths["summary_markdown"].exists()
    summary_payload = json.loads(paths["summary_json"].read_text())
    assert summary_payload["genes_after_gating"] == 4
    assert len(summary_payload["genes"]) == 4

    # Full per-gene reports only for BRAF/RET (+ any FDR-significant gene);
    # never for every scanned gene.
    full_report_genes = genes_needing_full_report(result)
    assert {"BRAF", "RET"} <= set(full_report_genes)
    assert set(full_report_genes) <= scanned_symbols
    for gene_symbol in full_report_genes:
        payload = json.loads(paths["gene_reports"][gene_symbol]["json"].read_text())
        assert payload["events"]
        assert all(row["tumor_type"] == "Glioma" for row in payload["events"])
        assert all(row["oncotree_code"] == "PA" for row in payload["events"])

    # Cohort scan writes the complete per-gene Markdown report, not only its
    # JSON payload. Assert report content from the generated artifact so a
    # future wiring regression cannot silently leave gene_reports/<gene>/report.md
    # empty or stale.
    braf_report = paths["gene_reports"]["BRAF"]["markdown"].read_text()
    assert braf_report.startswith("# BRAF real-data fusion benchmark: test_cohort_study")
    assert "- Protein-fusion records found: 7" in braf_report
    assert (
        "![Domain retention diagram](visualizations/domain_retention_outliers.svg)" in braf_report
    )


def test_cohort_scan_never_crashes_on_one_malformed_gene(mock_session, tmp_path):
    """A gene whose per-gene SV fetch raises must be recorded as a failed
    outcome, not abort the whole scan."""
    original_side_effect = mock_session.post.side_effect

    def _flaky_post(url, json=None, **kwargs):
        if url.endswith("/structural-variant/fetch") and 9001 in json["entrezGeneIds"]:
            raise ConnectionError("simulated transient failure for FAKE1")
        return original_side_effect(url, json=json, **kwargs)

    mock_session.post.side_effect = _flaky_post

    result = run_cohort_scan(
        _STUDY_ID,
        min_distinct_patients=5,
        n_permutations=50,
        adaptive=True,
        n_permutations_small=10,
        cache_dir=tmp_path / "cache",
        session=mock_session,
    )

    outcomes_by_gene = {outcome.gene_symbol: outcome for outcome in result.gene_outcomes}
    assert outcomes_by_gene["FAKE1"].status == "failed"
    assert "ConnectionError" in outcomes_by_gene["FAKE1"].error
    # The other three genes still completed successfully.
    assert outcomes_by_gene["BRAF"].status == "ok"
    assert outcomes_by_gene["RET"].status == "ok"
    assert outcomes_by_gene["FAKE2"].status == "ok"


def test_cohort_scan_gracefully_skips_gene_genome_nexus_cannot_resolve(mock_session):
    """A gated gene whose canonical transcript Genome Nexus has no mapping
    for at all must be recorded as unresolved, not crash the scan."""

    def _post_missing_fake1(url, json=None, **kwargs):
        response = MagicMock(status_code=200)
        if url.endswith("/structuralvariant-genes/fetch"):
            response.json.return_value = _recurrence_records()
        elif url.endswith("/mutations/fetch"):
            response.json.return_value = []
        elif url.endswith("/clinical-data/fetch"):
            response.json.return_value = [
                {"sampleId": sample, "clinicalAttributeId": attribute, "value": value}
                for sample in json["ids"]
                for attribute, value in [("CANCER_TYPE", "Glioma"), ("ONCOTREE_CODE", "PA")]
            ]
        elif url.endswith("/structural-variant/fetch"):
            entrez_ids = json["entrezGeneIds"]
            gene_symbol = next(
                symbol for symbol, (eid, *_r) in _GENE_SPECS.items() if eid in entrez_ids
            )
            response.json.return_value = _sv_calls_for_gene(gene_symbol)
        elif url.endswith("/ensembl/canonical-transcript/hgnc"):
            payloads = []
            for symbol in json:
                if symbol == "FAKE1":
                    continue  # Genome Nexus has no mapping for this gene
                _eid, _patients, protein_id, domains = _GENE_SPECS[symbol]
                payloads.append(_canonical_payload(symbol, protein_id, domains or []))
            response.json.return_value = payloads
        else:
            raise AssertionError(f"unexpected POST url: {url}")
        return response

    mock_session.post.side_effect = _post_missing_fake1

    result = run_cohort_scan(
        _STUDY_ID,
        min_distinct_patients=5,
        n_permutations=50,
        adaptive=True,
        n_permutations_small=10,
        session=mock_session,
    )

    outcomes_by_gene = {outcome.gene_symbol: outcome for outcome in result.gene_outcomes}
    assert outcomes_by_gene["FAKE1"].config_source == "unresolved"
    assert outcomes_by_gene["FAKE1"].status == "failed"
    assert result.unresolved_gene_count == 1
    assert any("FAKE1" in warning for warning in result.warnings)
    # Everything else still ran fine.
    assert outcomes_by_gene["BRAF"].status == "ok"
    assert outcomes_by_gene["FAKE2"].status == "ok"


@pytest.mark.parametrize("gene_symbol", ["BRAF", "RET"])
def test_cohort_scan_never_overrides_curated_config(mock_session, gene_symbol):
    """Curated genes must always resolve to the checked-in YAML config, not
    an auto-generated one, even though Genome Nexus is queried for both."""
    result = run_cohort_scan(
        _STUDY_ID,
        min_distinct_patients=5,
        n_permutations=50,
        adaptive=True,
        n_permutations_small=10,
        session=mock_session,
    )
    outcome = next(o for o in result.gene_outcomes if o.gene_symbol == gene_symbol)
    assert outcome.config_source == "curated"


def test_genome_nexus_client_accepts_injected_session_for_testing():
    """Sanity check the test's own mocking assumption: GenomeNexusClient
    accepts a session so cohort-scan's single shared client can be mocked."""
    client = GenomeNexusClient(session=MagicMock())
    assert client.session is not None
    _fetch = cbioportal_api.fetch_structural_variants  # exercised via run_cohort_scan above
    assert callable(_fetch)


@pytest.mark.parametrize("failed_fetch", [None, "expression", "mutation"])
def test_cohort_scan_fetches_optional_evidence_and_reports_partial_failures(
    mock_session, monkeypatch, failed_fetch
):
    study_config = StudyConfig(
        study_ids=[_STUDY_ID],
        all_sample_list_template="{study_id}_eligible",
        mrna_expression_profile_template="{study_id}_expression",
    )
    monkeypatch.setattr(scan_module, "load_study_config", lambda _: study_config)
    original_post = mock_session.post.side_effect

    def evidence_post(url, json=None, **kwargs):
        if url.endswith("/molecular-data/fetch"):
            if failed_fetch == "expression":
                raise requests.ConnectionError("expression unavailable")
            assert json["sampleListId"] == f"{_STUDY_ID}_eligible"
            response = MagicMock(status_code=200)
            response.json.return_value = [
                {"sampleId": call["sampleId"], "value": index + 2.0}
                for index, call in enumerate(_sv_calls_for_gene("BRAF"))
            ] + [{"sampleId": "NEG-1", "value": 0.0}, {"sampleId": "NEG-2", "value": 1.0}]
            return response
        if url.endswith("/mutations/fetch"):
            if failed_fetch == "mutation":
                raise requests.ConnectionError("mutation unavailable")
            response = MagicMock(status_code=200)
            response.json.return_value = [{"sampleId": "NEG-1", "proteinChange": "V600E"}]
            return response
        return original_post(url, json=json, **kwargs)

    mock_session.post.side_effect = evidence_post
    result = run_cohort_scan(
        _STUDY_ID,
        max_genes=1,
        n_permutations=5,
        adaptive=False,
        algorithm_names=["frequency", "expression_association", "mutation_cooccurrence"],
        cbioportal_base_url="https://custom.example/api",
        session=mock_session,
    )
    outcome = result.gene_outcomes[0]
    assert outcome.status == "ok", outcome.error
    results = {r.Algorithm: r for r in outcome.run.results}
    expression = results["expression_association"]
    mutation = results["mutation_cooccurrence"]
    if failed_fetch == "expression":
        assert "fusion_positive_vs_negative" not in expression.Summary
        assert any("expression unavailable" in warning for warning in outcome.run.warnings)
    else:
        assert expression.Summary["fusion_positive_vs_negative"]["n_fusion_positive"] == 7
        assert expression.Summary["fusion_positive_vs_negative"]["n_fusion_negative"] == 2
        assert any(row["algorithm"] == "expression_association" for row in result.fdr_rows)
    if failed_fetch == "mutation":
        assert mutation.Summary["targets"] == []
        assert any("mutation unavailable" in warning for warning in outcome.run.warnings)
        assert not any(row["algorithm"] == "mutation_cooccurrence" for row in result.fdr_rows)
    else:
        assert mutation.Summary["targets"][0]["contingency_table"] == [[0, 7], [1, 1]]
        assert any(row["algorithm"] == "mutation_cooccurrence" for row in result.fdr_rows)
    evidence_urls = [
        call.args[0]
        for call in mock_session.mock_calls
        if call.args
        and isinstance(call.args[0], str)
        and any(
            part in call.args[0] for part in ["/mutations/", "/sample-lists/", "/molecular-data/"]
        )
    ]
    assert len(evidence_urls) == 3
    assert all(url.startswith("https://custom.example/api/") for url in evidence_urls)


def test_cohort_scan_does_not_fetch_unrequested_evidence(mock_session):
    result = run_cohort_scan(
        _STUDY_ID,
        max_genes=1,
        n_permutations=5,
        algorithm_names=["frequency"],
        session=mock_session,
    )
    assert result.gene_outcomes[0].status == "ok"
    assert not any(
        any(part in str(call) for part in ["/mutations/", "/sample-lists/", "/molecular-data/"])
        for call in mock_session.mock_calls
    )


def test_manhattan_svg_is_wired_into_all_three_cohort_scan_outputs(mock_session, tmp_path):
    """The genome-wide Manhattan/volcano SVG is written to
    ``cohort_scan/manhattan.svg``, embedded via Markdown image syntax in
    ``summary.md``, and embedded (not just linked) as a figure in
    ``summary.pdf`` -- for whatever genes actually scanned in this run, no
    gene names hardcoded here."""
    result = run_cohort_scan(
        _STUDY_ID,
        min_distinct_patients=5,
        n_permutations=50,
        adaptive=True,
        n_permutations_small=10,
        cache_dir=tmp_path / "cache",
        session=mock_session,
    )
    paths = write_cohort_scan_outputs(result, tmp_path / "runs", pdf=True)

    manhattan_svg_path = paths["manhattan_svg"]
    assert manhattan_svg_path.exists()
    assert manhattan_svg_path.name == "manhattan.svg"
    assert manhattan_svg_path.parent == paths["summary_markdown"].parent
    svg_text = manhattan_svg_path.read_text()
    assert svg_text.startswith("<svg")
    # BRAF and RET both produced FDR-tested p-values in this fixture (see
    # the module-level FDR assertions above), so both must appear as
    # plotted, gene-attributed points.
    assert 'data-gene="BRAF"' in svg_text
    assert 'data-gene="RET"' in svg_text

    markdown_text = paths["summary_markdown"].read_text()
    assert "![Genome-wide fusion-hotspot summary plot](manhattan.svg)" in markdown_text

    from pypdf import PdfReader

    pdf_text = "".join(
        page.extract_text() or "" for page in PdfReader(str(paths["summary_pdf"])).pages
    )
    assert "manhattan" in pdf_text
