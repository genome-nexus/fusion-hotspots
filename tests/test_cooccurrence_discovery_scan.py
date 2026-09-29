"""Wiring of co-occurrence discovery into the cohort scan and its outputs."""

from __future__ import annotations

from unittest.mock import MagicMock

from cfh.cohort.cooccurrence_data import load_discovery_inputs
from cfh.cohort.cooccurrence_discovery import SampleRecord
from cfh.cohort.outputs import select_discovery_rows


def test_written_rows_keep_significant_pairs_and_top_hits_per_gene():
    rows = [
        {"fusion_gene": "A", "cmh_p_value": p, "cmh_q_value": q}
        for p, q in [
            (1e-9, 1e-6),
            (0.2, 0.5),
            (0.3, 0.6),
            (0.4, 0.7),
            (0.5, 0.8),
            (0.6, 0.9),
            (0.7, 0.95),
        ]
    ] + [{"fusion_gene": "A", "cmh_p_value": None, "cmh_q_value": None}]
    selected = select_discovery_rows(rows)
    assert [row["cmh_p_value"] for row in selected] == [1e-9, 0.2, 0.3, 0.4, 0.5]


def test_loader_builds_samples_panels_and_alterations_from_api():
    session = MagicMock()

    def _post(url, json=None, **kwargs):
        response = MagicMock()
        if url.endswith("/gene-panel-data/fetch"):
            response.json.return_value = [
                {"sampleId": "S1", "profiled": True, "genePanelId": "P"},
                {"sampleId": "S2", "profiled": False, "genePanelId": "P"},
            ]
        elif "/mutations/fetch" in url:
            response.json.return_value = [{"sampleId": "S1", "entrezGeneId": 3845}]
        elif "/discrete-copy-number/fetch" in url:
            response.json.return_value = [
                {"sampleId": "S1", "entrezGeneId": 673, "alteration": 2},
                {"sampleId": "S1", "entrezGeneId": 3845, "alteration": -1},  # shallow: ignored
            ]
        else:
            raise AssertionError(url)
        return response

    def _get(url, params=None, **kwargs):
        response = MagicMock()
        if url.endswith("/gene-panels/P"):
            response.json.return_value = {
                "genes": [
                    {"entrezGeneId": 3845, "hugoGeneSymbol": "KRAS"},
                    {"entrezGeneId": 673, "hugoGeneSymbol": "BRAF"},
                ]
            }
        elif url.endswith("/clinical-data"):
            response.json.return_value = [{"sampleId": "S1", "value": "LUAD"}]
        else:
            raise AssertionError(url)
        return response

    session.post.side_effect = _post
    session.get.side_effect = _get
    inputs = load_discovery_inputs(
        "study",
        sv_profile_id="sv",
        mutation_profile_id="mut",
        cna_profile_id="cna",
        sample_list_id="all",
        base_url="https://portal.example/api",
        session=session,
    )
    by_id = {record.sample_id: record for record in inputs.samples}
    assert by_id["S1"] == SampleRecord("S1", "LUAD", "P", "P", "P")
    assert by_id["S2"].sv_panel is None and by_id["S2"].stratum == "UNKNOWN"
    assert inputs.alterations == {("KRAS", "mutation"): {"S1"}, ("BRAF", "amplification"): {"S1"}}
    assert inputs.panel_genes == {"P": {"KRAS", "BRAF"}}
    assert any("no ONCOTREE_CODE" in w for w in inputs.warnings)
