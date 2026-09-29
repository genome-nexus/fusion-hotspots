import pandas as pd

from cfh.ingestion import clinical_parser, sv_parser
from cfh.normalization.event_normalizer import DEFAULT_COHORT, normalize


def _events(sv_fixture_file, clinical_sample_fixture, cohort=DEFAULT_COHORT):
    raw = sv_parser.parse_sv_file(sv_fixture_file)
    clinical = clinical_parser.parse_clinical_sample(clinical_sample_fixture)
    return normalize(raw, clinical, cohort)


def _by_sample(events, sample_id):
    return next(e for e in events if e.Sample_id == sample_id)


def test_output_length_matches_input_row_count(sv_fixture_file, clinical_sample_fixture):
    raw = sv_parser.parse_sv_file(sv_fixture_file)
    events = _events(sv_fixture_file, clinical_sample_fixture)
    assert len(events) == len(raw)


def test_intragenic_deletion_is_classified_as_deletion_via_distinct_path(
    sv_fixture_file, clinical_sample_fixture
):
    events = _events(sv_fixture_file, clinical_sample_fixture)
    deletion_event = _by_sample(events, "SAMPLE-002")
    fusion_event = _by_sample(events, "SAMPLE-001")

    assert deletion_event.Event_class == "deletion"
    # Distinct code path, not just relabeling: intragenic events are never
    # flagged as a two-gene protein fusion.
    assert deletion_event.Is_protein_fusion is False
    assert fusion_event.Is_protein_fusion is True


def test_inversion_and_translocation_enum_values(sv_fixture_file, clinical_sample_fixture):
    events = _events(sv_fixture_file, clinical_sample_fixture)
    assert _by_sample(events, "SAMPLE-003").Event_class == "inversion"
    assert _by_sample(events, "SAMPLE-004").Event_class == "translocation"


def test_antisense_flag(sv_fixture_file, clinical_sample_fixture):
    events = _events(sv_fixture_file, clinical_sample_fixture)
    antisense_event = _by_sample(events, "SAMPLE-005")
    assert antisense_event.Is_antisense is True
    assert antisense_event.Frame_status == "out-of-frame"


def test_mid_exon_fusion_parses_and_classifies_as_fusion(sv_fixture_file, clinical_sample_fixture):
    events = _events(sv_fixture_file, clinical_sample_fixture)
    event = _by_sample(events, "SAMPLE-006")
    assert event.Event_class == "fusion"
    assert event.Frame_status == "in-frame"


def test_missing_annotation_row_has_unknown_frame_status_but_known_orientation(
    sv_fixture_file, clinical_sample_fixture
):
    events = _events(sv_fixture_file, clinical_sample_fixture)
    event = _by_sample(events, "SAMPLE-007")
    assert event.Frame_status == "unknown"
    assert event.Five_prime_gene == "BRAF"


def test_imprecise_numeric_row_still_normalizes(sv_fixture_file, clinical_sample_fixture):
    events = _events(sv_fixture_file, clinical_sample_fixture)
    event = _by_sample(events, "SAMPLE-008")
    assert event.Split_read_support is None  # unparsable "IMPPRECISE" -> not coerced
    assert event.Confidence_class == "low"


def test_negative_read_support_sentinels_are_missing(sv_fixture_file, clinical_sample_fixture):
    raw = sv_parser.parse_sv_file(sv_fixture_file).iloc[:1].copy()
    raw.loc[raw.index[0], "Tumor_Paired_End_Read_Count"] = -1
    raw.loc[raw.index[0], "Tumor_Split_Read_Count"] = -1
    clinical = clinical_parser.parse_clinical_sample(clinical_sample_fixture)

    event = normalize(raw, clinical, DEFAULT_COHORT)[0]
    assert event.Paired_end_read_support is None
    assert event.Split_read_support is None
    assert event.Total_read_support is None


def test_tumor_variant_count_is_populated_from_the_live_api_column(
    sv_fixture_file, clinical_sample_fixture
):
    # Tumor_Variant_Count has no column in the offline data_sv.txt schema
    # (sv_parser.EXPECTED_COLUMNS) -- it only arrives via the cBioPortal
    # structural-variant API path (see cbioportal_api._API_TO_NORMALIZED_COLUMNS),
    # so it's added directly to the raw frame here to exercise that column.
    raw = sv_parser.parse_sv_file(sv_fixture_file).iloc[:1].copy()
    raw["Tumor_Variant_Count"] = 35
    clinical = clinical_parser.parse_clinical_sample(clinical_sample_fixture)

    event = normalize(raw, clinical, DEFAULT_COHORT)[0]
    assert event.Tumor_variant_count == 35


def test_negative_tumor_variant_count_sentinel_is_missing(sv_fixture_file, clinical_sample_fixture):
    raw = sv_parser.parse_sv_file(sv_fixture_file).iloc[:1].copy()
    raw["Tumor_Variant_Count"] = -1
    clinical = clinical_parser.parse_clinical_sample(clinical_sample_fixture)

    event = normalize(raw, clinical, DEFAULT_COHORT)[0]
    assert event.Tumor_variant_count is None


def test_ambiguous_orientation_is_never_guessed(sv_fixture_file, clinical_sample_fixture):
    events = _events(sv_fixture_file, clinical_sample_fixture)
    event = _by_sample(events, "SAMPLE-009")
    assert event.Five_prime_gene is None
    assert event.Three_prime_gene is None
    assert event.Frame_status == "unknown"


def test_shared_patient_different_samples_produce_two_distinct_events(
    sv_fixture_file, clinical_sample_fixture
):
    events = _events(sv_fixture_file, clinical_sample_fixture)
    e1 = _by_sample(events, "SAMPLE-001")
    e2 = _by_sample(events, "SAMPLE-002")
    assert e1.Patient_id == e2.Patient_id == "PATIENT-001"
    assert e1.Sample_id != e2.Sample_id
    assert e1.Event_id != e2.Event_id
    assert e1 is not e2


def test_missing_sample_id_row_is_kept_not_dropped(sv_fixture_file, clinical_sample_fixture):
    events = _events(sv_fixture_file, clinical_sample_fixture)
    missing = [e for e in events if e.Sample_id is None]
    assert len(missing) == 1


def test_float_valued_read_support_column_is_coerced_to_int_not_dropped(
    sv_fixture_file, clinical_sample_fixture
):
    # Tumor_Paired_End_Read_Count has one missing value (SAMPLE-007) among
    # otherwise-valid integers, so sv_parser's output column upcasts to
    # float64 and every valid value arrives here as e.g. 10.0, not 10.
    raw = sv_parser.parse_sv_file(sv_fixture_file)
    assert raw["Tumor_Paired_End_Read_Count"].dtype.kind == "f"

    events = _events(sv_fixture_file, clinical_sample_fixture)
    event = _by_sample(events, "SAMPLE-001")
    assert event.Paired_end_read_support == 10
    assert isinstance(event.Paired_end_read_support, int)
    assert event.Total_read_support == (event.Paired_end_read_support + event.Split_read_support)


def test_cohort_is_caller_supplied_not_hardcoded(sv_fixture_file, clinical_sample_fixture):
    msk_events = _events(sv_fixture_file, clinical_sample_fixture, cohort="msk_impact_50k_2026")
    tcga_events = _events(
        sv_fixture_file, clinical_sample_fixture, cohort="tcga_pan_can_atlas_2018"
    )

    assert {e.Cohort for e in msk_events} == {"msk_impact_50k_2026"}
    assert {e.Cohort for e in tcga_events} == {"tcga_pan_can_atlas_2018"}


def test_oncotree_code_differs_across_samples_from_clinical_join(
    sv_fixture_file, clinical_sample_fixture
):
    events = _events(sv_fixture_file, clinical_sample_fixture)
    glioma_event = _by_sample(events, "SAMPLE-001")
    melanoma_event = _by_sample(events, "SAMPLE-002")

    assert glioma_event.Oncotree_code == "PA"
    assert glioma_event.Tumor_type == "Glioma"
    assert melanoma_event.Oncotree_code == "SKCM"
    assert melanoma_event.Tumor_type == "Melanoma"
    assert glioma_event.Oncotree_code != melanoma_event.Oncotree_code


def test_sample_missing_from_clinical_data_gets_none_oncotree_without_raising():
    raw = pd.DataFrame.from_records(
        [
            {
                "Sample_Id": "SAMPLE-NOT-IN-CLINICAL",
                "Site1_Hugo_Symbol": "BRAF",
                "Site2_Hugo_Symbol": "KIAA1549",
                "Site1_Chromosome": "7",
                "Site2_Chromosome": "7",
                "Site1_Position": 140487000,
                "Site2_Position": 138500000,
                "Site2_Effect_On_Frame": "in-frame",
                "Connection_Type": "5to3",
                "Breakpoint_Type": "PRECISE",
                "Annotation": "KIAA1549-BRAF fusion, exon16:exon9, in-frame",
                "Event_Info": "KIAA1549-BRAF fusion",
                "Tumor_Split_Read_Count": 10,
                "Tumor_Paired_End_Read_Count": 5,
                "Source_row_number": 1,
                "Parse_warnings": None,
            }
        ]
    )
    clinical = pd.DataFrame.from_records(
        [{"Sample_id": "SAMPLE-OTHER", "Patient_id": "PATIENT-OTHER", "Oncotree_code": "SKCM"}]
    )

    events = normalize(raw, clinical, DEFAULT_COHORT)

    assert len(events) == 1
    assert events[0].Oncotree_code is None
    assert events[0].Tumor_type is None
    assert events[0].Patient_id is None


def test_real_api_na_frame_falls_back_to_event_info_and_explicit_fusion_order():
    raw = pd.DataFrame.from_records(
        [
            {
                "Sample_Id": "SAMPLE-REAL",
                "Site1_Hugo_Symbol": "KIAA1549",
                "Site2_Hugo_Symbol": "BRAF",
                "Site2_Effect_On_Frame": "NA",
                "Connection_Type": "3to3",
                "Event_Info": "Protein Fusion: in frame  {KIAA1549:BRAF}",
                "Source_row_number": 1,
            }
        ]
    )

    event = normalize(raw, None, "real-study")[0]

    assert event.Frame_status == "in-frame"
    assert event.Five_prime_gene == "KIAA1549"
    assert event.Three_prime_gene == "BRAF"


def test_parenthesized_event_info_supplies_explicit_fusion_order():
    raw = pd.DataFrame.from_records(
        [
            {
                "Sample_Id": "SAMPLE-PAREN",
                "Site1_Hugo_Symbol": "KIAA1549",
                "Site2_Hugo_Symbol": "BRAF",
                "Site2_Effect_On_Frame": "NA",
                "Connection_Type": "5to5",
                "Event_Info": "Protein fusion: in frame (KIAA1549-BRAF)",
                "Source_row_number": 1,
            }
        ]
    )

    event = normalize(raw, None, "real-study")[0]

    assert event.Frame_status == "in-frame"
    assert event.Five_prime_gene == "KIAA1549"
    assert event.Three_prime_gene == "BRAF"


def test_tcga_pan_cancer_event_info_supplies_explicit_fusion_order():
    raw = pd.DataFrame.from_records(
        [
            {
                "Sample_Id": "TCGA-DE-A0Y2-01",
                "Site1_Hugo_Symbol": "MACF1",
                "Site2_Hugo_Symbol": "BRAF",
                "Site1_Chromosome": "1",
                "Site2_Chromosome": "7",
                "Site2_Effect_On_Frame": "in-frame",
                "Connection_Type": None,
                "Event_Info": "MACF1-BRAF Fusion",
                "Source_row_number": 1,
            }
        ]
    )

    event = normalize(raw, None, "thca_tcga_pan_can_atlas_2018")[0]

    assert event.Is_protein_fusion is True
    assert event.Frame_status == "in-frame"
    assert event.Five_prime_gene == "MACF1"
    assert event.Three_prime_gene == "BRAF"


def test_raw_patient_identity_fills_missing_clinical_identity():
    from cfh.ingestion.cbioportal_api import structural_variants_to_dataframe

    raw = structural_variants_to_dataframe(
        [
            {
                "sampleId": "S",
                "patientId": "RAW-PATIENT",
                "site1HugoSymbol": "PARTNER",
                "site2HugoSymbol": "BRAF",
                "eventInfo": "Protein Fusion: in frame {PARTNER:BRAF}",
            }
        ]
    )
    assert normalize(raw, None, "study")[0].Patient_id == "RAW-PATIENT"
    clinical = pd.DataFrame([{"Sample_id": "S", "Patient_id": float("nan")}], dtype=object)
    assert normalize(raw, clinical, "study")[0].Patient_id == "RAW-PATIENT"
    clinical.loc[0, "Patient_id"] = "CLINICAL-PATIENT"
    assert normalize(raw, clinical, "study")[0].Patient_id == "CLINICAL-PATIENT"
