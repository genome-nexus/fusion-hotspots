import inspect
from unittest.mock import MagicMock

import pytest

from cfh.genes.registry import load_gene_config
from cfh.ingestion import cbioportal_api


def test_fetch_structural_variants_posts_expected_body_without_real_network():
    mock_session = MagicMock()
    mock_response = MagicMock()
    mock_response.json.return_value = [{"sampleId": "SAMPLE-001"}]
    mock_response.raise_for_status.return_value = None
    mock_session.post.return_value = mock_response

    gene_config = load_gene_config("braf")
    result = cbioportal_api.fetch_structural_variants(
        [gene_config.entrez_gene_id],
        [cbioportal_api.DEFAULT_SV_MOLECULAR_PROFILE_ID],
        session=mock_session,
    )

    assert result == [{"sampleId": "SAMPLE-001"}]
    mock_session.post.assert_called_once()
    _, kwargs = mock_session.post.call_args
    called_url = mock_session.post.call_args.args[0]
    assert called_url == f"{cbioportal_api.DEFAULT_BASE_URL}/structural-variant/fetch"
    assert kwargs["json"] == {
        "entrezGeneIds": [673],
        "molecularProfileIds": [cbioportal_api.DEFAULT_SV_MOLECULAR_PROFILE_ID],
    }


def test_fetch_structural_variants_supports_arbitrary_genes_without_real_network():
    """The client itself must stay gene-agnostic -- any Entrez id list works."""
    mock_session = MagicMock()
    mock_session.post.return_value.json.return_value = []

    cbioportal_api.fetch_structural_variants(
        [7157, 673, 5290],  # TP53, BRAF, PIK3CA -- arbitrary, not special-cased
        ["some_other_study_structural_variants"],
        session=mock_session,
    )

    _, kwargs = mock_session.post.call_args
    assert kwargs["json"]["entrezGeneIds"] == [7157, 673, 5290]
    assert kwargs["json"]["molecularProfileIds"] == ["some_other_study_structural_variants"]


def test_molecular_profile_ids_has_no_msk_specific_default():
    """Regression: molecular_profile_ids must be required, not silently
    defaulted to the MSK-IMPACT profile -- a caller supplying only gene IDs
    must not end up silently querying MSK-IMPACT data.
    """
    signature = inspect.signature(cbioportal_api.fetch_structural_variants)
    molecular_profile_ids_param = signature.parameters["molecular_profile_ids"]
    assert molecular_profile_ids_param.default is inspect.Parameter.empty

    with pytest.raises(TypeError):
        cbioportal_api.fetch_structural_variants([673])


def test_fetch_structural_variants_retries_transient_service_failure():
    mock_session = MagicMock()
    unavailable = MagicMock(status_code=503)
    recovered = MagicMock(status_code=200)
    recovered.json.return_value = [{"sampleId": "RECOVERED"}]
    mock_session.post.side_effect = [unavailable, recovered]

    result = cbioportal_api.fetch_structural_variants(
        [5979],
        ["study_structural_variants"],
        session=mock_session,
        max_retries=1,
        backoff_seconds=0,
    )

    assert result == [{"sampleId": "RECOVERED"}]
    assert mock_session.post.call_count == 2
    recovered.raise_for_status.assert_called_once_with()


def test_fetch_structural_variants_retries_more_than_three_times_by_default(monkeypatch):
    """Regression: a genome-wide cohort scan makes hundreds of sequential
    calls to this endpoint and can run alongside other callers hitting the
    same public API, so the default retry budget must tolerate more than a
    couple of transient 503s in a row, not give up after the first few.
    """
    monkeypatch.setattr(cbioportal_api.time, "sleep", lambda _seconds: None)
    mock_session = MagicMock()
    unavailable = MagicMock(status_code=503)
    recovered = MagicMock(status_code=200)
    recovered.json.return_value = [{"sampleId": "RECOVERED"}]
    mock_session.post.side_effect = [unavailable, unavailable, unavailable, unavailable, recovered]

    result = cbioportal_api.fetch_structural_variants(
        [5979],
        ["study_structural_variants"],
        session=mock_session,
    )

    assert result == [{"sampleId": "RECOVERED"}]
    assert mock_session.post.call_count == 5


def test_retry_sleep_seconds_is_capped_and_jittered_above_the_uncapped_value():
    """The exponential backoff must never exceed the cap, but a large enough
    attempt number must actually hit that cap (not silently keep growing) --
    and the jitter must only ever add delay, never subtract it."""
    uncapped = cbioportal_api._retry_sleep_seconds(1.0, attempt=1)
    assert 2.0 <= uncapped <= 2.0 * 1.25

    max_backoff = cbioportal_api._MAX_BACKOFF_SECONDS
    capped = cbioportal_api._retry_sleep_seconds(1.0, attempt=10)
    assert max_backoff <= capped <= max_backoff * 1.25


def test_structural_variant_api_rows_are_adapted_to_production_normalizer_schema():
    rows = cbioportal_api.structural_variants_to_dataframe(
        [
            {
                "sampleId": "SAMPLE-001",
                "site1HugoSymbol": "KIAA1549",
                "site2HugoSymbol": "BRAF",
                "site2Position": 140493152,
                "site2EffectOnFrame": "NA",
                "connectionType": "3to3",
                "tumorVariantCount": 35,
                "eventInfo": "Protein Fusion: in frame  {KIAA1549:BRAF}",
                "patientId": "PATIENT-001",
            }
        ]
    )

    assert rows.loc[0, "Sample_Id"] == "SAMPLE-001"
    assert rows.loc[0, "Site2_Effect_On_Frame"] == "NA"
    assert rows.loc[0, "Event_Info"] == "Protein Fusion: in frame  {KIAA1549:BRAF}"
    assert rows.loc[0, "Source_row_number"] == 1
    assert rows.loc[0, "Tumor_Variant_Count"] == 35
    assert rows.loc[0, "Extra_fields"]["patientId"] == "PATIENT-001"


def test_fetch_molecular_data_posts_expected_body_with_sample_list_id():
    mock_session = MagicMock()
    mock_response = MagicMock()
    mock_response.json.return_value = [{"sampleId": "SAMPLE-001", "value": 1.23}]
    mock_response.raise_for_status.return_value = None
    mock_session.post.return_value = mock_response

    result = cbioportal_api.fetch_molecular_data(
        [673],
        "thca_tcga_pan_can_atlas_2018_rna_seq_v2_mrna_median_Zscores",
        sample_list_id="thca_tcga_pan_can_atlas_2018_all",
        session=mock_session,
    )

    assert result == [{"sampleId": "SAMPLE-001", "value": 1.23}]
    called_url = mock_session.post.call_args.args[0]
    assert called_url == (
        f"{cbioportal_api.DEFAULT_BASE_URL}/molecular-profiles/"
        "thca_tcga_pan_can_atlas_2018_rna_seq_v2_mrna_median_Zscores/molecular-data/fetch"
    )
    _, kwargs = mock_session.post.call_args
    assert kwargs["json"] == {
        "entrezGeneIds": [673],
        "sampleListId": "thca_tcga_pan_can_atlas_2018_all",
    }


def test_fetch_molecular_data_supports_explicit_sample_ids_instead_of_a_sample_list():
    mock_session = MagicMock()
    mock_session.post.return_value.json.return_value = []

    cbioportal_api.fetch_molecular_data(
        [5979],
        "some_profile",
        sample_ids=["S1", "S2"],
        session=mock_session,
    )

    _, kwargs = mock_session.post.call_args
    assert kwargs["json"] == {"entrezGeneIds": [5979], "sampleIds": ["S1", "S2"]}


def test_fetch_molecular_data_requires_exactly_one_of_sample_list_id_or_sample_ids():
    with pytest.raises(ValueError, match="exactly one of"):
        cbioportal_api.fetch_molecular_data([673], "profile")
    with pytest.raises(ValueError, match="exactly one of"):
        cbioportal_api.fetch_molecular_data([673], "profile", sample_list_id="a", sample_ids=["b"])


def test_molecular_data_to_expression_by_sample_adapts_real_shape():
    records = [
        {"sampleId": "TCGA-1", "value": 0.5},
        {"sampleId": "TCGA-2", "value": -1.25},
    ]
    assert cbioportal_api.molecular_data_to_expression_by_sample(records) == {
        "TCGA-1": 0.5,
        "TCGA-2": -1.25,
    }


@pytest.mark.parametrize(
    "malformed_record",
    [
        {"sampleId": "TCGA-1"},  # missing value
        {"value": 1.0},  # missing sampleId
        {"sampleId": "TCGA-1", "value": None},
        {"sampleId": "TCGA-1", "value": "not-a-number"},
        {"sampleId": "TCGA-1", "value": float("nan")},
        "not-a-dict",
    ],
)
def test_molecular_data_to_expression_by_sample_skips_malformed_rows_without_raising(
    malformed_record,
):
    records = [malformed_record, {"sampleId": "TCGA-GOOD", "value": 2.0}]
    result = cbioportal_api.molecular_data_to_expression_by_sample(records)
    assert result == {"TCGA-GOOD": 2.0}


@pytest.mark.network
def test_fetch_structural_variants_real_network_call():
    """Excluded from default `pytest -m "not network"` runs."""
    gene_config = load_gene_config("braf")
    result = cbioportal_api.fetch_structural_variants(
        [gene_config.entrez_gene_id],
        [cbioportal_api.DEFAULT_SV_MOLECULAR_PROFILE_ID],
    )
    assert isinstance(result, list)
