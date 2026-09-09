"""Client for the cBioPortal structural-variant REST API.

Gene selection (which Entrez ids to fetch) is always caller-supplied --
typically from a ``GeneConfig``'s ``entrez_gene_id`` -- so this module
stays generic across genes. The real-network call in
:func:`fetch_structural_variants` is only ever exercised by tests marked
``@pytest.mark.network`` (excluded from the default ``pytest`` run);
everything else in this module is plain, mockable request-building logic.
"""

from __future__ import annotations

import random
import time
import warnings
from typing import Any, Iterable

import pandas as pd
import requests

from cfh.ingestion.sv_parser import OUTPUT_COLUMNS

DEFAULT_BASE_URL = "https://www.cbioportal.org/api"
DEFAULT_STUDY_ID = "msk_impact_50k_2026"
DEFAULT_SV_MOLECULAR_PROFILE_ID = "msk_impact_50k_2026_structural_variants"
DEFAULT_MUTATION_MOLECULAR_PROFILE_ID = "msk_impact_50k_2026_mutations"
DEFAULT_CNA_MOLECULAR_PROFILE_ID = "msk_impact_50k_2026_gistic"
DEFAULT_SAMPLE_LIST_ID = "msk_impact_50k_2026_all"
_RETRYABLE_STATUS_CODES = {429, 500, 502, 503, 504}
_MAX_BACKOFF_SECONDS = 30.0


def _retry_sleep_seconds(backoff_seconds: float, attempt: int) -> float:
    """Exponential backoff, capped, with up to 25% random jitter.

    A genome-wide cohort scan makes hundreds of sequential calls to this
    same endpoint, and may run concurrently with other callers hitting the
    same public API -- an uncapped exponential wait can grow unreasonably
    long, and every caller retrying on the exact same schedule synchronizes
    their retries into new bursts of contention instead of spreading them
    out. Both are transient-infrastructure-load mitigations, independent of
    any particular gene or query.
    """
    capped = min(backoff_seconds * (2**attempt), _MAX_BACKOFF_SECONDS)
    return capped * (1 + random.uniform(0, 0.25))


_API_TO_NORMALIZED_COLUMNS = {
    "sampleId": "Sample_Id",
    "site1HugoSymbol": "Site1_Hugo_Symbol",
    "site1Chromosome": "Site1_Chromosome",
    "site1Position": "Site1_Position",
    "site2HugoSymbol": "Site2_Hugo_Symbol",
    "site2Chromosome": "Site2_Chromosome",
    "site2Position": "Site2_Position",
    "site2EffectOnFrame": "Site2_Effect_On_Frame",
    "tumorSplitReadCount": "Tumor_Split_Read_Count",
    "tumorPairedEndReadCount": "Tumor_Paired_End_Read_Count",
    "tumorVariantCount": "Tumor_Variant_Count",
    "svStatus": "SV_Status",
    "ncbiBuild": "NCBI_Build",
    "connectionType": "Connection_Type",
    "breakpointType": "Breakpoint_Type",
    "annotation": "Annotation",
    "eventInfo": "Event_Info",
}


def fetch_structural_variants(
    entrez_gene_ids: Iterable[int],
    molecular_profile_ids: Iterable[str],
    *,
    base_url: str = DEFAULT_BASE_URL,
    session: "requests.Session | None" = None,
    timeout: float = 30,
    max_retries: int = 6,
    backoff_seconds: float = 1.0,
) -> list[dict]:
    """POST to ``/structural-variant/fetch`` and return the parsed JSON body.

    ``molecular_profile_ids`` is required and has no default: which cohort's
    SV profile to query is always caller-supplied (e.g. from ingestion
    config), never silently defaulted to a specific study like MSK-IMPACT.
    ``DEFAULT_SV_MOLECULAR_PROFILE_ID`` remains available for callers that
    do want the MSK-IMPACT 50k profile, but it's opt-in, not automatic.

    A genome-wide cohort scan calls this hundreds of times in one run, so
    the retry/backoff defaults here are more patient than a single one-off
    call needs (see :func:`_retry_sleep_seconds`) -- this endpoint is where
    transient upstream 503s under sustained per-gene load are actually
    observed in practice.
    """
    session = session or requests.Session()
    url = f"{base_url.rstrip('/')}/structural-variant/fetch"
    body: dict[str, Any] = {
        "entrezGeneIds": list(entrez_gene_ids),
        "molecularProfileIds": list(molecular_profile_ids),
    }
    attempt = 0
    while True:
        response = session.post(url, json=body, timeout=timeout)
        if response.status_code not in _RETRYABLE_STATUS_CODES or attempt >= max_retries:
            break
        time.sleep(_retry_sleep_seconds(backoff_seconds, attempt))
        attempt += 1
    response.raise_for_status()
    return response.json()


def fetch_mutations(
    entrez_gene_ids: Iterable[int],
    molecular_profile_ids: Iterable[str],
    *,
    base_url: str = DEFAULT_BASE_URL,
    session: "requests.Session | None" = None,
    timeout: float = 30,
    max_retries: int = 6,
    backoff_seconds: float = 1.0,
) -> list[dict]:
    """POST to ``/mutations/fetch`` and return the parsed JSON body.

    Mirrors :func:`fetch_structural_variants`'s request shape and retry
    behavior -- ``molecularProfileIds`` + ``entrezGeneIds`` -- for the
    point-mutation evidence layer. ``molecular_profile_ids`` has no default
    for the same reason as the structural-variant fetch: which cohort's
    mutation profile to query is always caller-supplied, never silently
    defaulted to a specific study.
    """
    session = session or requests.Session()
    url = f"{base_url.rstrip('/')}/mutations/fetch"
    body: dict[str, Any] = {
        "entrezGeneIds": list(entrez_gene_ids),
        "molecularProfileIds": list(molecular_profile_ids),
    }
    attempt = 0
    while True:
        response = session.post(url, json=body, timeout=timeout)
        if response.status_code not in _RETRYABLE_STATUS_CODES or attempt >= max_retries:
            break
        time.sleep(_retry_sleep_seconds(backoff_seconds, attempt))
        attempt += 1
    response.raise_for_status()
    return response.json()


def fetch_discrete_copy_number(
    entrez_gene_ids: Iterable[int],
    molecular_profile_id: str,
    sample_list_id: str,
    *,
    event_type: str = "ALL",
    base_url: str = DEFAULT_BASE_URL,
    session: "requests.Session | None" = None,
    timeout: float = 30,
    max_retries: int = 6,
    backoff_seconds: float = 1.0,
) -> list[dict]:
    """POST to ``/molecular-profiles/{molecularProfileId}/discrete-copy-number/fetch``.

    Unlike mutations/structural-variants, this cBioPortal endpoint is
    per-molecular-profile (not multi-profile) and requires a
    ``sampleListId`` -- both ``molecular_profile_id`` and ``sample_list_id``
    are always caller-supplied, never defaulted to a specific study or
    gene. ``event_type`` follows cBioPortal's ``DiscreteCopyNumberEventType``
    enum (``AMP``, ``HOMDEL``, ``GAIN``, ``HETLOSS``, ``DIPLOID``,
    ``HOMDEL_AND_AMP``, ``ALL``); ``"ALL"`` returns every alteration state
    so the caller -- not this client -- decides which states count as a
    "hit" (see ``cfh.normalization.alteration_normalizer``).
    """
    session = session or requests.Session()
    url = (
        f"{base_url.rstrip('/')}/molecular-profiles/{molecular_profile_id}"
        "/discrete-copy-number/fetch"
    )
    body: dict[str, Any] = {
        "sampleListId": sample_list_id,
        "entrezGeneIds": list(entrez_gene_ids),
    }
    attempt = 0
    while True:
        response = session.post(
            url,
            json=body,
            params={"discreteCopyNumberEventType": event_type},
            timeout=timeout,
        )
        if response.status_code not in _RETRYABLE_STATUS_CODES or attempt >= max_retries:
            break
        time.sleep(_retry_sleep_seconds(backoff_seconds, attempt))
        attempt += 1
    response.raise_for_status()
    return response.json()


def fetch_sample_list_ids(
    sample_list_id: str,
    *,
    base_url: str = DEFAULT_BASE_URL,
    session: "requests.Session | None" = None,
    timeout: float = 30,
    max_retries: int = 6,
    backoff_seconds: float = 1.0,
) -> list[str]:
    """GET ``/sample-lists/{sampleListId}/sample-ids``: every sample ID in a
    named cBioPortal sample list (e.g. a study's "_all" list). Used as the
    cohort-wide 2x2 background universe for the mutation/CNA co-occurrence
    test -- gene-agnostic and study-agnostic, ``sample_list_id`` is always
    caller-supplied.
    """
    session = session or requests.Session()
    url = f"{base_url.rstrip('/')}/sample-lists/{sample_list_id}/sample-ids"
    attempt = 0
    while True:
        response = session.get(url, timeout=timeout)
        if response.status_code not in _RETRYABLE_STATUS_CODES or attempt >= max_retries:
            break
        time.sleep(_retry_sleep_seconds(backoff_seconds, attempt))
        attempt += 1
    response.raise_for_status()
    return response.json()


def fetch_structural_variant_genes(
    study_ids: Iterable[str],
    *,
    base_url: str = DEFAULT_BASE_URL,
    session: "requests.Session | None" = None,
    timeout: float = 30,
    max_retries: int = 6,
    backoff_seconds: float = 1.0,
) -> list[dict]:
    """POST to ``/structuralvariant-genes/fetch`` for cohort-wide SV gene recurrence.

    Returns one ``AlterationCountByGene``-shaped dict per gene that has at
    least one structural-variant record anywhere in the requested
    studies -- every such gene in one call, each carrying its own
    ``numberOfAlteredCases`` (distinct-patient count) and ``totalCount``
    (raw SV record count). This is the whole-cohort recurrence signal used
    to gate a genome-wide scan down to recurrently-altered genes, as
    opposed to :func:`fetch_structural_variants`, which fetches per-event
    SV records for an already-selected set of Entrez gene ids.
    """
    session = session or requests.Session()
    url = f"{base_url.rstrip('/')}/structuralvariant-genes/fetch"
    body = {"studyIds": list(study_ids)}
    attempt = 0
    while True:
        response = session.post(url, json=body, timeout=timeout)
        if response.status_code not in _RETRYABLE_STATUS_CODES or attempt >= max_retries:
            break
        time.sleep(_retry_sleep_seconds(backoff_seconds, attempt))
        attempt += 1
    response.raise_for_status()
    return response.json()


_API_OUTPUT_COLUMNS = [*dict.fromkeys([*OUTPUT_COLUMNS, *_API_TO_NORMALIZED_COLUMNS.values()])]
"""``OUTPUT_COLUMNS`` plus every API field the live fetch path can supply.

``tumorVariantCount`` (unlike the split/paired-end read counts) has no
column in the offline ``data_sv.txt`` schema that ``sv_parser.OUTPUT_COLUMNS``
describes, so it must be added here rather than in ``sv_parser`` -- adding it
there would make every offline SV fixture missing that column (all of them)
warn about it.
"""


def structural_variants_to_dataframe(calls: Iterable[dict]) -> pd.DataFrame:
    """Adapt cBioPortal camelCase API objects to the production SV schema."""
    records = []
    for row_number, call in enumerate(calls, start=1):
        record = {
            destination: call.get(source)
            for source, destination in _API_TO_NORMALIZED_COLUMNS.items()
        }
        record["Extra_fields"] = {
            key: value for key, value in call.items() if key not in _API_TO_NORMALIZED_COLUMNS
        }
        record["Source_row_number"] = row_number
        record["Parse_warnings"] = None
        records.append(record)
    return pd.DataFrame.from_records(records, columns=_API_OUTPUT_COLUMNS)


def fetch_sample_tumor_types(
    study_id: str,
    sample_ids: Iterable[str],
    *,
    base_url: str = DEFAULT_BASE_URL,
    session: requests.Session | None = None,
    timeout: float = 30,
) -> pd.DataFrame:
    """Fetch optional sample annotations without changing patient/grouping metadata."""
    columns = ["Sample_id", "Tumor_type", "Oncotree_code"]
    ids = sorted(set(sample_ids))
    if not ids:
        return pd.DataFrame(columns=columns)
    session = session or requests.Session()
    try:
        response = session.post(
            f"{base_url.rstrip('/')}/studies/{study_id}/clinical-data/fetch",
            params={"clinicalDataType": "SAMPLE"},
            json={"ids": ids, "attributeIds": ["CANCER_TYPE", "ONCOTREE_CODE"]},
            timeout=timeout,
        )
        response.raise_for_status()
        records = response.json()
    except (requests.RequestException, ValueError) as exc:
        warnings.warn(f"Sample tumor annotations unavailable for {study_id}: {exc}", stacklevel=2)
        return pd.DataFrame(columns=columns)
    fields = {"CANCER_TYPE": "Tumor_type", "ONCOTREE_CODE": "Oncotree_code"}
    rows: dict[str, dict] = {}
    for record in records:
        sample_id = record.get("sampleId")
        field = fields.get(record.get("clinicalAttributeId"))
        if sample_id in ids and field:
            row = rows.setdefault(sample_id, dict.fromkeys(columns))
            row["Sample_id"] = sample_id
            row[field] = record.get("value") or None
    return pd.DataFrame(rows.values(), columns=columns)
