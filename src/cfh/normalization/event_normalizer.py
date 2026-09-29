"""Normalize raw cBioPortal SV rows into typed :class:`FusionEvent` objects.

Orientation (5'/3') and reading frame are only ever reported when the
source data says so explicitly; when the connection type is missing or
unrecognized, both are left ambiguous rather than guessed.
"""

from __future__ import annotations

import re
from typing import Any

import pandas as pd

from cfh.model.fusion_event import FusionEvent

DEFAULT_COHORT = "msk_impact_50k_2026"
"""Example cohort id for the MSK-IMPACT 50k ingestion path.

Not applied automatically inside :func:`normalize` -- cohort provenance is
caller-supplied so the same normalizer works for other cBioPortal cohorts
(e.g. a future TCGA PanCancer Atlas ingestion) without mislabeling them.
"""

_FRAME_STATUS_MAP = {"in-frame": "in-frame", "out-of-frame": "out-of-frame"}
_KNOWN_CONNECTION_TYPES = {"5to3", "3to5"}


def _text_blob(row: dict[str, Any]) -> str:
    return f"{row.get('Annotation') or ''} {row.get('Event_Info') or ''}".lower()


def _detect_frame_status(row: dict[str, Any]) -> str:
    effect = (row.get("Site2_Effect_On_Frame") or "").strip().lower()
    if effect in _FRAME_STATUS_MAP:
        return _FRAME_STATUS_MAP[effect]
    text = _text_blob(row)
    if "out-of-frame" in text or "out of frame" in text:
        return "out-of-frame"
    if "in-frame" in text or "in frame" in text:
        return "in-frame"
    return "unknown"


def _is_antisense(row: dict[str, Any]) -> bool:
    return "antisense" in _text_blob(row)


def _fusion_order_from_event_info(row: dict[str, Any]) -> tuple[str | None, str | None]:
    """Read the explicit 5'→3' gene order from cBioPortal fusion text.

    Real MSK-IMPACT API records use all four connection types because those
    values describe breakpoint-end directions, not necessarily transcript
    order.  Protein-fusion records separately encode transcript order as
    ``{FIVE_PRIME:THREE_PRIME}`` or ``(FIVE_PRIME-THREE_PRIME)``.
    """
    event_info = str(row.get("Event_Info") or "")
    if "fusion" not in event_info.lower():
        return None, None
    match = re.search(r"\{\s*([^:{}\s]+)\s*:\s*([^{}:\s]+)\s*\}", event_info)
    if match is None:
        match = re.search(r"\(\s*([^()\s-]+)\s*-\s*([^()\s-]+)\s*\)", event_info)
    if match is None:
        # TCGA PanCancer Atlas SV profiles use ``GENE1-GENE2 Fusion`` while
        # MSK profiles wrap the ordered genes in braces or parentheses.
        match = re.search(
            r"(?:^|\s)([A-Za-z0-9.]+)\s*-\s*([A-Za-z0-9.]+)\s+fusion\b",
            event_info,
            re.IGNORECASE,
        )
    if match is None:
        return None, None
    return match.group(1), match.group(2)


def _orientation(row: dict[str, Any]) -> tuple[str | None, str | None]:
    explicit_five_prime, explicit_three_prime = _fusion_order_from_event_info(row)
    if explicit_five_prime is not None:
        return explicit_five_prime, explicit_three_prime
    connection = (row.get("Connection_Type") or "").strip().lower()
    if connection == "5to3":
        return row.get("Site1_Hugo_Symbol"), row.get("Site2_Hugo_Symbol")
    if connection == "3to5":
        return row.get("Site2_Hugo_Symbol"), row.get("Site1_Hugo_Symbol")
    return None, None


def _classify_intragenic(row: dict[str, Any]) -> dict[str, Any]:
    """Same gene on both breakpoint sites: an intragenic rearrangement.

    This is a distinct code path from partner-gene fusions: intragenic
    events are never treated as a two-gene protein fusion.
    """
    text = _text_blob(row)
    if "inversion" in text:
        event_class = "inversion"
    elif "deletion" in text:
        event_class = "deletion"
    else:
        event_class = "unknown"

    gene = row.get("Site1_Hugo_Symbol")
    connection = (row.get("Connection_Type") or "").strip().lower()
    orientation_known = connection in _KNOWN_CONNECTION_TYPES
    frame_status = _detect_frame_status(row) if orientation_known else "unknown"

    return {
        "Event_class": event_class,
        "Fusion_name": f"{gene} intragenic {event_class}" if gene else None,
        "Is_protein_fusion": False,
        "Five_prime_gene": gene if orientation_known else None,
        "Three_prime_gene": gene if orientation_known else None,
        "Frame_status": frame_status,
    }


def _classify_intergenic(row: dict[str, Any]) -> dict[str, Any]:
    """Two different genes: a partner-gene fusion, inversion, or translocation."""
    text = _text_blob(row)
    chrom1, chrom2 = row.get("Site1_Chromosome"), row.get("Site2_Chromosome")

    if "inversion" in text:
        event_class = "inversion"
    elif "translocation" in text:
        event_class = "translocation"
    elif chrom1 and chrom2 and str(chrom1) != str(chrom2):
        event_class = "translocation"
    else:
        event_class = "fusion"

    five_prime, three_prime = _orientation(row)
    orientation_known = five_prime is not None or three_prime is not None
    frame_status = _detect_frame_status(row) if orientation_known else "unknown"

    gene1, gene2 = row.get("Site1_Hugo_Symbol"), row.get("Site2_Hugo_Symbol")
    fusion_name = f"{gene1}-{gene2}" if gene1 and gene2 else None

    return {
        "Event_class": event_class,
        "Fusion_name": fusion_name,
        "Is_protein_fusion": True,
        "Five_prime_gene": five_prime,
        "Three_prime_gene": three_prime,
        "Frame_status": frame_status,
    }


def _confidence_class(row: dict[str, Any]) -> str:
    breakpoint_type = (row.get("Breakpoint_Type") or "").strip().lower()
    if breakpoint_type == "precise":
        return "high"
    if breakpoint_type == "imprecise":
        return "low"
    return "unknown"


def _as_int(value: Any) -> int | None:
    """Coerce a numeric-valued scalar to ``int``.

    ``sv_parser``'s output DataFrame upcasts an all-numeric column that also
    contains missing values (``None``) to ``float64``, so a valid read-support
    count like ``10`` can arrive here as ``10.0`` (or ``numpy.float64(10.0)``).
    That must still become integer ``10``, not be dropped as "not an int".
    Truly missing (``None``/``NaN``) or non-numeric values still become
    ``None``.
    """
    if value is None:
        return None
    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return None
    if not numeric.is_integer():
        return None
    return int(numeric)


def _as_chromosome(value: Any) -> str | None:
    """Coerce a raw chromosome/build scalar to a plain string, or ``None``.

    Mirrors :func:`_as_int`'s missing-value handling for a non-numeric
    field: the row's own ``None``-normalization already catches most cases,
    but a stray ``NaN`` (e.g. from a DataFrame column upcast) must still
    become ``None`` rather than the literal string ``"nan"``.
    """
    if value is None:
        return None
    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass
    return str(value)


def _clinical_lookup(clinical_df: pd.DataFrame | None) -> dict[str, dict[str, Any]]:
    if clinical_df is None or "Sample_id" not in clinical_df.columns:
        return {}
    lookup: dict[str, dict[str, Any]] = {}
    for _, row in clinical_df.iterrows():
        sample_id = row.get("Sample_id")
        if sample_id:
            lookup[sample_id] = row.to_dict()
    return lookup


def normalize(
    raw_rows_df: pd.DataFrame,
    clinical_df: pd.DataFrame | None,
    cohort: str,
    *,
    sequencing_panel_id: str | None = None,
) -> list[FusionEvent]:
    """Normalize raw SV rows (see ``sv_parser``) into ``FusionEvent`` objects.

    ``cohort`` identifies the source cBioPortal study (e.g. ``"msk_impact_50k_2026"``)
    and is always caller-supplied -- never inferred or hardcoded -- so the same
    normalizer works unmodified for other cohorts sharing this SV schema.
    ``sequencing_panel_id`` is a fallback used only when the clinical data
    doesn't already carry a per-sample panel id.

    Output length always equals input row count: normalization never
    drops or merges rows (deduplication, if any, happens downstream).
    """
    clinical_lookup = _clinical_lookup(clinical_df)
    events: list[FusionEvent] = []

    for _, raw_row in raw_rows_df.iterrows():
        row = raw_row.where(pd.notnull(raw_row), None).to_dict()

        sample_id = row.get("Sample_Id")
        clinical_row = clinical_lookup.get(sample_id, {}) if sample_id else {}
        patient_id = clinical_row.get("Patient_id")
        if pd.isna(patient_id) or not patient_id:
            patient_id = row.get("Patient_Id")

        gene1, gene2 = row.get("Site1_Hugo_Symbol"), row.get("Site2_Hugo_Symbol")
        if gene1 and gene2 and gene1 == gene2:
            classification = _classify_intragenic(row)
        else:
            classification = _classify_intergenic(row)

        source_row_number = row.get("Source_row_number")
        event_id = f"EVT-{sample_id or 'UNKNOWN'}-{source_row_number}"

        paired_end_read_support = _as_int(row.get("Tumor_Paired_End_Read_Count"))
        split_read_support = _as_int(row.get("Tumor_Split_Read_Count"))
        tumor_variant_count = _as_int(row.get("Tumor_Variant_Count"))
        # cBioPortal uses negative counts (normally -1) as an
        # unavailable-value sentinel, not as biological measurements.
        paired_end_read_support = (
            paired_end_read_support
            if paired_end_read_support is not None and paired_end_read_support >= 0
            else None
        )
        split_read_support = (
            split_read_support
            if split_read_support is not None and split_read_support >= 0
            else None
        )
        tumor_variant_count = (
            tumor_variant_count
            if tumor_variant_count is not None and tumor_variant_count >= 0
            else None
        )
        # A partial sum when only one of the two components is available
        # (e.g. paired-end present but split absent) still gets summed here
        # -- Total_read_support is "whatever read-support evidence exists,"
        # not a guarantee both components contributed.
        read_support_values = [
            value for value in (paired_end_read_support, split_read_support) if value is not None
        ]

        events.append(
            FusionEvent(
                Event_id=event_id,
                Cohort=cohort,
                Sequencing_panel_id=clinical_row.get("Sequencing_panel_id") or sequencing_panel_id,
                Sample_id=sample_id,
                Patient_id=patient_id,
                Site1_gene=gene1,
                Site2_gene=gene2,
                Site1_chromosome=_as_chromosome(row.get("Site1_Chromosome")),
                Site1_position=_as_int(row.get("Site1_Position")),
                Site2_chromosome=_as_chromosome(row.get("Site2_Chromosome")),
                Site2_position=_as_int(row.get("Site2_Position")),
                Reference_build=_as_chromosome(row.get("NCBI_Build")),
                Five_prime_gene=classification["Five_prime_gene"],
                Three_prime_gene=classification["Three_prime_gene"],
                Fusion_name=classification["Fusion_name"],
                Event_class=classification["Event_class"],
                Connection_type=row.get("Connection_Type"),
                Frame_status=classification["Frame_status"],
                Is_protein_fusion=classification["Is_protein_fusion"],
                Is_antisense=_is_antisense(row),
                Confidence_class=_confidence_class(row),
                Paired_end_read_support=paired_end_read_support,
                Split_read_support=split_read_support,
                Total_read_support=(sum(read_support_values) if read_support_values else None),
                Tumor_variant_count=tumor_variant_count,
                Site1_description=(
                    f"{gene1}:{row.get('Site1_Chromosome')}:{row.get('Site1_Position')}"
                    if gene1
                    else None
                ),
                Site2_description=(
                    f"{gene2}:{row.get('Site2_Chromosome')}:{row.get('Site2_Position')}"
                    if gene2
                    else None
                ),
                Annotation=row.get("Annotation"),
                Event_info=row.get("Event_Info"),
                Source_row_number=source_row_number,
                Tumor_type=clinical_row.get("Tumor_type"),
                Oncotree_code=clinical_row.get("Oncotree_code"),
            )
        )

    return events
