"""Normalize cBioPortal mutation/discrete-copy-number API rows into
:class:`~cfh.model.alteration_event.AlterationEvent` records.

This lives alongside (and never touches) the structural-variant/
``FusionEvent`` normalization pipeline in
:mod:`cfh.normalization.event_normalizer`. It follows the same
never-crash-on-a-malformed-row convention already used throughout
ingestion/normalization: a row missing a required field is skipped and
recorded as a warning, never raised.
"""

from __future__ import annotations

from typing import Any

from cfh.model.alteration_event import AlterationEvent

CNA_ALTERATION_TYPE_BY_CODE: dict[int, str] = {
    2: "cna_amp",
    1: "cna_gain",
    0: "cna_diploid",
    -1: "cna_hetloss",
    -2: "cna_del",
}
"""cBioPortal's discrete copy-number ``alteration`` integer code, mapped to
this repo's ``Alteration_type`` vocabulary (matches the
``DiscreteCopyNumberEventType`` enum cBioPortal's REST API itself uses)."""


def normalize_mutations(
    calls: list[dict[str, Any]], gene_symbol: str, cohort: str
) -> tuple[list[AlterationEvent], list[str]]:
    """Adapt cBioPortal ``/mutations/fetch`` rows (see
    :func:`cfh.ingestion.cbioportal_api.fetch_mutations`) queried for a
    single gene.

    ``gene_symbol`` is caller-supplied -- the gene the Entrez id in the
    fetch was queried for -- because the API's default (SUMMARY-projection)
    mutation rows do not reliably carry the gene's Hugo symbol; this mirrors
    the "caller already knows which gene it queried for" pattern already
    used for structural-variant fusion mapping.
    """
    events: list[AlterationEvent] = []
    warnings: list[str] = []
    for row_number, call in enumerate(calls, start=1):
        try:
            sample_id = call.get("sampleId")
            if not sample_id:
                raise ValueError("missing sampleId")
            events.append(
                AlterationEvent(
                    Sample_id=str(sample_id),
                    Patient_id=call.get("patientId"),
                    Cohort=cohort,
                    Gene=gene_symbol.upper(),
                    Alteration_type="point_mutation",
                    Protein_change=call.get("proteinChange"),
                    Mutation_type=call.get("mutationType"),
                    Source_row_number=row_number,
                )
            )
        except Exception as exc:  # noqa: BLE001 - a malformed row must never crash a fetch
            warnings.append(
                f"Skipped malformed mutation row {row_number}: {type(exc).__name__}: {exc}"
            )
    return events, warnings


def normalize_discrete_copy_number(
    calls: list[dict[str, Any]], gene_symbol: str, cohort: str
) -> tuple[list[AlterationEvent], list[str]]:
    """Adapt cBioPortal ``discrete-copy-number/fetch`` rows (see
    :func:`cfh.ingestion.cbioportal_api.fetch_discrete_copy_number`) queried
    for a single gene. ``gene_symbol`` is caller-supplied for the same
    reason as :func:`normalize_mutations`.
    """
    events: list[AlterationEvent] = []
    warnings: list[str] = []
    for row_number, call in enumerate(calls, start=1):
        try:
            sample_id = call.get("sampleId")
            if not sample_id:
                raise ValueError("missing sampleId")
            alteration_code = call.get("alteration")
            if alteration_code is None:
                raise ValueError("missing alteration code")
            alteration_type = CNA_ALTERATION_TYPE_BY_CODE.get(int(alteration_code))
            if alteration_type is None:
                raise ValueError(f"unrecognized alteration code {alteration_code!r}")
            events.append(
                AlterationEvent(
                    Sample_id=str(sample_id),
                    Patient_id=call.get("patientId"),
                    Cohort=cohort,
                    Gene=gene_symbol.upper(),
                    Alteration_type=alteration_type,
                    Source_row_number=row_number,
                )
            )
        except Exception as exc:  # noqa: BLE001 - a malformed row must never crash a fetch
            warnings.append(
                f"Skipped malformed discrete copy-number row {row_number}: "
                f"{type(exc).__name__}: {exc}"
            )
    return events, warnings
