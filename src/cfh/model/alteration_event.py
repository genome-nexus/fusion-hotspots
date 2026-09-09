from typing import Optional

from pydantic import BaseModel, ConfigDict


class AlterationEvent(BaseModel):
    """A single normalized per-sample point-mutation or copy-number alteration call.

    Distinct from :class:`~cfh.model.fusion_event.FusionEvent` (structural-
    variant/fusion records): this represents cBioPortal
    ``mutations``/``discrete-copy-number`` rows, used by the mutation/CNA
    co-occurrence evidence layer to test a gene's fusion status against
    another alteration across a cohort. Lives alongside ``FusionEvent``
    rather than extending or otherwise touching its schema.
    """

    model_config = ConfigDict(extra="forbid")

    Sample_id: str
    Patient_id: Optional[str] = None
    Cohort: Optional[str] = None
    Gene: str
    Alteration_type: str
    """One of ``point_mutation``, ``cna_amp``, ``cna_gain``, ``cna_diploid``,
    ``cna_hetloss``, ``cna_del`` -- see
    :mod:`cfh.normalization.alteration_normalizer`."""
    Protein_change: Optional[str] = None
    Mutation_type: Optional[str] = None
    Source_row_number: Optional[int] = None
