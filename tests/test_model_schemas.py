from cfh.model.algorithm_result import AlgorithmResult
from cfh.model.fusion_event import FusionEvent
from cfh.model.fusion_feature import FusionFeature

FUSION_EVENT_FIELDS = {
    "Event_id",
    "Cohort",
    "Sequencing_panel_id",
    "Sample_id",
    "Patient_id",
    "Site1_gene",
    "Site2_gene",
    "Site1_chromosome",
    "Site1_position",
    "Site2_chromosome",
    "Site2_position",
    "Reference_build",
    "Five_prime_gene",
    "Three_prime_gene",
    "Fusion_name",
    "Event_class",
    "Connection_type",
    "Frame_status",
    "Is_protein_fusion",
    "Is_antisense",
    "Confidence_class",
    "Paired_end_read_support",
    "Split_read_support",
    "Total_read_support",
    "Tumor_variant_count",
    "Site1_description",
    "Site2_description",
    "Annotation",
    "Event_info",
    "Source_row_number",
    "Tumor_type",
    "Oncotree_code",
}

FUSION_FEATURE_FIELDS = {
    "Event_id",
    "Gene",
    "Role",
    "Transcript_id",
    "Breakpoint_exon",
    "Breakpoint_intron",
    "Retained_exons",
    "Lost_exons",
    "Junction_position_aa",
    "Retained_domains",
    "Lost_domains",
    "Disrupted_domains",
    "Domain_retention_flags",
    "Domain_retention_details",
}

ALGORITHM_RESULT_FIELDS = {
    "Algorithm",
    "Algorithm_version",
    "Parameters",
    "Summary",
    "Tables",
    "Warnings",
    "Created_at",
    "Input_fingerprint",
}


def test_fusion_event_field_set_matches_prd():
    assert set(FusionEvent.model_fields) == FUSION_EVENT_FIELDS


def test_fusion_feature_field_set_matches_prd():
    assert set(FusionFeature.model_fields) == FUSION_FEATURE_FIELDS


def test_algorithm_result_field_set_matches_prd():
    assert set(AlgorithmResult.model_fields) == ALGORITHM_RESULT_FIELDS
