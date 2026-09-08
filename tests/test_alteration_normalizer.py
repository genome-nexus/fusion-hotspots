from cfh.normalization.alteration_normalizer import (
    normalize_discrete_copy_number,
    normalize_mutations,
)


def test_normalize_mutations_maps_real_cbioportal_fields():
    calls = [
        {
            "sampleId": "P-0000217-T01-IM3",
            "patientId": "P-0000217",
            "entrezGeneId": 673,
            "proteinChange": "V600E",
            "mutationType": "Missense_Mutation",
        }
    ]

    events, warnings = normalize_mutations(calls, "BRAF", "msk_impact_50k_2026")

    assert warnings == []
    assert len(events) == 1
    event = events[0]
    assert event.Sample_id == "P-0000217-T01-IM3"
    assert event.Patient_id == "P-0000217"
    assert event.Gene == "BRAF"
    assert event.Alteration_type == "point_mutation"
    assert event.Protein_change == "V600E"
    assert event.Mutation_type == "Missense_Mutation"
    assert event.Source_row_number == 1


def test_normalize_mutations_skips_malformed_rows_without_crashing():
    calls = [
        {"sampleId": "GOOD-1", "proteinChange": "V600E"},
        {"patientId": "PATIENT-2"},  # missing sampleId entirely
        None,  # not even a dict
        {"sampleId": "GOOD-2", "proteinChange": "K601E"},
    ]

    events, warnings = normalize_mutations(calls, "BRAF", "cohort")

    assert [event.Sample_id for event in events] == ["GOOD-1", "GOOD-2"]
    assert len(warnings) == 2
    assert "row 2" in warnings[0]
    assert "row 3" in warnings[1]


def test_normalize_discrete_copy_number_maps_alteration_codes():
    calls = [
        {"sampleId": "S1", "patientId": "P1", "alteration": 2},
        {"sampleId": "S2", "alteration": -2},
        {"sampleId": "S3", "alteration": 0},
        {"sampleId": "S4", "alteration": 1},
        {"sampleId": "S5", "alteration": -1},
    ]

    events, warnings = normalize_discrete_copy_number(calls, "BRAF", "cohort")

    assert warnings == []
    assert [event.Alteration_type for event in events] == [
        "cna_amp",
        "cna_del",
        "cna_diploid",
        "cna_gain",
        "cna_hetloss",
    ]
    assert all(event.Gene == "BRAF" for event in events)


def test_normalize_discrete_copy_number_skips_malformed_rows_without_crashing():
    calls = [
        {"sampleId": "GOOD-1", "alteration": 2},
        {"sampleId": "MISSING-ALTERATION"},
        {"alteration": 2},  # missing sampleId
        {"sampleId": "BAD-CODE", "alteration": 999},
        {"sampleId": "BAD-TYPE", "alteration": "not-an-int"},
    ]

    events, warnings = normalize_discrete_copy_number(calls, "BRAF", "cohort")

    assert [event.Sample_id for event in events] == ["GOOD-1"]
    assert len(warnings) == 4
