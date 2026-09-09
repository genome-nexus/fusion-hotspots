from typing import Optional

from pydantic import BaseModel, ConfigDict


class FusionEvent(BaseModel):
    """A single normalized structural-variant / fusion event for one sample."""

    model_config = ConfigDict(extra="forbid")

    Event_id: str
    Cohort: str
    Sequencing_panel_id: Optional[str] = None
    Sample_id: Optional[str] = None
    Patient_id: Optional[str] = None
    Site1_gene: Optional[str] = None
    Site2_gene: Optional[str] = None
    Site1_chromosome: Optional[str] = None
    Site1_position: Optional[int] = None
    Site2_chromosome: Optional[str] = None
    Site2_position: Optional[int] = None
    Reference_build: Optional[str] = None
    """Genome assembly the raw ``Site1_position``/``Site2_position`` values
    are reported against (e.g. ``"GRCh37"``), copied verbatim from the
    source SV record's ``NCBI_Build``/``ncbiBuild`` field. Never inferred or
    normalized to a canonical spelling here -- a caller that pools positions
    across cohorts must compare this field before assuming two positions
    share a coordinate system, since different cBioPortal studies are not
    guaranteed to report the same build."""
    Five_prime_gene: Optional[str] = None
    Three_prime_gene: Optional[str] = None
    Fusion_name: Optional[str] = None
    Event_class: Optional[str] = None
    Connection_type: Optional[str] = None
    Frame_status: Optional[str] = None
    Is_protein_fusion: Optional[bool] = None
    Is_antisense: Optional[bool] = None
    Confidence_class: Optional[str] = None
    Paired_end_read_support: Optional[int] = None
    Split_read_support: Optional[int] = None
    Total_read_support: Optional[int] = None
    Tumor_variant_count: Optional[int] = None
    Site1_description: Optional[str] = None
    Site2_description: Optional[str] = None
    Annotation: Optional[str] = None
    Event_info: Optional[str] = None
    Source_row_number: Optional[int] = None
    Tumor_type: Optional[str] = None
    Oncotree_code: Optional[str] = None
