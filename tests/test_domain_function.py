"""Sequence-level functional-domain checks (cfh.domain_function)."""

from __future__ import annotations

import gzip
from collections import Counter
from pathlib import Path

import pytest

from cfh.domain_function.afdb import ParentStructure, ss_class
from cfh.domain_function.assess import assess_record, hmm_report
from cfh.domain_function.metrics import (
    NOT_ASSESSED,
    align_parents,
    junction_structure_report,
    key_residue_report,
    map_structure,
    place_chimera,
    schema_report,
)
from cfh.domain_function.pfam_hmm import (
    HmmUnavailable,
    best_overlapping,
    gathering_thresholds,
    load_hmm,
    residues_in_inserts,
    scan,
)
from cfh.domain_function.record import (
    ChimeraRecord,
    ParentDomain,
    add_row,
    record_from_dict,
    record_to_dict,
)
from cfh.domain_function.run import shortlist, tier_counts
from cfh.domain_function.tiers import (
    LIKELY_BROKEN,
    LIKELY_INTACT,
    RULES,
    UNCERTAIN,
    assign_tier,
    disagreements,
)
from cfh.domain_function.uniprot import KeyFeature, UniProtEntry, is_key_feature, parse_entry

DATA = Path(__file__).parent / "data" / "domain_function"
ETS_HMM = DATA / "PF00178.hmm.gz"


def _fli1() -> str:
    lines = (DATA / "Q01543_FLI1.fasta").read_text().splitlines()
    return "".join(line for line in lines if not line.startswith(">"))


FLI1 = _fli1()
FLI1_ETS = (281, 361)  # UniProt DNA-binding ETS region


# --------------------------------------------------------------------------- HMMER
@pytest.fixture(scope="module")
def ets_hmm():
    return load_hmm(ETS_HMM)


def test_known_answer_parent_domain_scores_above_ga(ets_hmm):
    _, domain_ga = gathering_thresholds(ets_hmm)
    window = FLI1[FLI1_ETS[0] - 31 : FLI1_ETS[1] + 30]
    result = scan(ets_hmm, [("fli1", window)])["fli1"]
    hit, n_hits = best_overlapping(result, 31, 31 + 80)
    assert hit is not None and n_hits == 1
    assert hit.score > domain_ga + 50
    assert hit.match_state_coverage > 0.95
    assert hit.insert_positions == () or len(hit.insert_positions) < 5


def test_known_answer_large_internal_deletion_drops(ets_hmm):
    _, domain_ga = gathering_thresholds(ets_hmm)
    start, end = FLI1_ETS
    intact = FLI1[start - 31 : end + 30]
    # Delete 45 residues from the middle of the 81-residue ETS domain.
    deleted = FLI1[start - 31 : start + 15] + FLI1[start + 60 : end + 30]
    results = scan(ets_hmm, [("intact", intact), ("deleted", deleted)])
    full, _ = best_overlapping(results["intact"], 1, len(intact))
    cut, _ = best_overlapping(results["deleted"], 1, len(deleted))
    assert full is not None and full.score >= domain_ga
    assert cut is None or cut.score < domain_ga
    assert cut is None or cut.match_state_coverage < full.match_state_coverage


def test_unrelated_sequence_has_no_passing_hit(ets_hmm):
    _, domain_ga = gathering_thresholds(ets_hmm)
    result = scan(ets_hmm, [("poly", "GS" * 60)])["poly"]
    hit, _ = best_overlapping(result, 1, 120)
    assert hit is None or hit.score < domain_ga


def test_hmm_without_gathering_threshold_fails_closed(tmp_path):
    text = gzip.decompress(ETS_HMM.read_bytes()).decode()
    stripped = "\n".join(line for line in text.splitlines() if not line.startswith("GA "))
    path = tmp_path / "noga.hmm"
    path.write_text(stripped + "\n")
    with pytest.raises(HmmUnavailable):
        load_hmm(path)


def test_insert_state_count_is_none_without_hit():
    assert residues_in_inserts(None, [1, 2, 3]) is None


def test_hmm_report_without_hmm_is_not_assessed():
    record = _record(_PARENT, _PARTNER, five_last=30, three_first=31)
    assert hmm_report(record, None, [30, 31]) == {"hmm_status": NOT_ASSESSED}


# --------------------------------------------------------------------- synthetic records
_PARENT = "MSTEKLLRAGHPWVDEQNTYFKAILGSRDCPEMHWNKAVLTGEFSRQDYPLIAKTEHMNGWSA"
_PARTNER = "MSTDKLLRAGHPWIDEQNTYFRAILGSKDCPEMHWNKSVLTGEFARQDYPLLAKSEHMNGWTA"


def _record(
    five: str, three: str, *, five_last: int, three_first: int, hybrid: str = ""
) -> ChimeraRecord:
    protein = five[:five_last] + hybrid + three[three_first - 1 :]
    junction = five_last + len(hybrid)
    return ChimeraRecord(
        record_id="chim_test",
        pfam_family="PF00000",
        junction_class="test",
        class_reason="",
        five=ParentDomain("GENEA", "TXA", five, 5, 60, five_last),
        three=ParentDomain("GENEB", "TXB", three, 5, 60, three_first),
        hybrid_residue=hybrid,
        chimeric_protein=protein,
        domain_start=5,
        domain_end=junction + (60 - three_first + 1),
        junction_index=junction,
        pilot_indel_residues=None,
        pilot_indel_sequence="",
        pilot_alignment_parent="",
        pilot_alignment_identity=None,
    )


def test_origin_maps_each_chimera_residue_to_its_parent():
    record = _record(_PARENT, _PARTNER, five_last=30, three_first=31, hybrid="W")
    assert record.origin(30) == ("five", 30)
    assert record.origin(31) == ("hybrid", None)
    assert record.origin(32) == ("three", 31)
    assert record.chimeric_protein[31] == _PARTNER[30]


def test_clean_junction_has_no_inserted_or_deleted_columns():
    record = _record(_PARENT, _PARTNER, five_last=30, three_first=31)
    placed = place_chimera(record, align_parents(record))
    assert placed.inserted == [] and placed.deleted_columns == []


def test_tandem_overlap_is_counted_as_inserted_residues():
    record = _record(_PARENT, _PARTNER, five_last=30, three_first=26)
    placed = place_chimera(record, align_parents(record))
    assert len(placed.inserted) == 5
    assert [p for _, p in placed.inserted_sources] == [26, 27, 28, 29, 30]
    assert placed.deleted_columns == []


def test_skipped_columns_are_counted_as_deleted():
    record = _record(_PARENT, _PARTNER, five_last=30, three_first=38)
    placed = place_chimera(record, align_parents(record))
    assert placed.inserted == []
    assert len(placed.deleted_columns) == 7
    assert placed.deleted_residues["five"] == list(range(31, 38))


def _entry(gene: str, sequence: str, features: tuple[KeyFeature, ...]) -> UniProtEntry:
    return UniProtEntry(gene, f"P_{gene}", sequence, features)


def test_key_residues_not_assessed_without_uniprot():
    record = _record(_PARENT, _PARTNER, five_last=30, three_first=31)
    parents = align_parents(record)
    report = key_residue_report(
        record, parents, place_chimera(record, parents), {"five": None, "three": None}
    )
    assert report["key_features_assessed"] is False
    assert report["key_lost"] == NOT_ASSESSED
    assert report["junction_near_key"] == NOT_ASSESSED


def test_key_residue_lost_retained_and_swapped():
    record = _record(_PARENT, _PARTNER, five_last=30, three_first=38)
    parents = align_parents(record)
    placed = place_chimera(record, parents)
    five = _entry(
        "GENEA",
        _PARENT,
        (
            KeyFeature("Active site", 33, 33, ""),  # in the deleted columns
            KeyFeature("Binding site", 10, 10, "Zn(2+)"),  # kept 5' residue
        ),
    )
    # 3' residue 41 (K in partner, R in parent): chimera keeps the partner's residue.
    three = _entry("GENEB", _PARTNER, (KeyFeature("Site", 54, 54, ""),))
    report = key_residue_report(record, parents, placed, {"five": five, "three": three})
    assert report["key_features_assessed"] is True
    assert report["key_lost"] == 1
    assert "GENEA:Active site:33" in report["key_point_lost"]
    assert report["key_retained"] == 2
    assert report["junction_to_key_5p"] == 3
    assert report["junction_near_key"] is True


def test_key_residue_mapping_fails_closed_on_isoform_mismatch():
    record = _record(_PARENT, _PARTNER, five_last=30, three_first=31)
    parents = align_parents(record)
    unrelated = _entry("GENEA", "W" * 63, (KeyFeature("Active site", 10, 10, ""),))
    report = key_residue_report(
        record, parents, place_chimera(record, parents), {"five": unrelated, "three": None}
    )
    assert report["key_features_5p"] == NOT_ASSESSED
    assert report["key_features_assessed"] is False


def test_uniprot_feature_filter():
    assert is_key_feature({"type": "Active site"})
    assert not is_key_feature({"type": "Region", "description": "Disordered"})
    assert not is_key_feature({"type": "Site", "description": "Breakpoint for translocation"})
    assert not is_key_feature({"type": "Chain"})
    raw = {
        "primaryAccession": "Q1",
        "sequence": {"value": "MA"},
        "features": [
            {
                "type": "Binding site",
                "location": {"start": {"value": 2}, "end": {"value": 2}},
                "ligand": {"name": "Zn(2+)"},
            },
            {
                "type": "Motif",
                "location": {"start": {"value": None}, "end": {"value": 2}},
            },
        ],
    }
    entry = parse_entry("G", raw)
    assert entry is not None and [f.label for f in entry.features] == ["Binding site:2(Zn(2+))"]


def _structure(sequence: str, ss: str | None, contacts=(), rsa=None) -> ParentStructure:
    n = len(sequence)
    return ParentStructure(
        accession="P_X",
        model_id="AF-P_X-F1",
        model_version="6",
        sequence=sequence,
        first=1,
        plddt=tuple([90.0] * n),
        ss=ss,
        contacts=tuple(contacts),
        rsa=tuple(rsa or [0.05] * n),
    )


def test_ss_class_never_assumes_loop():
    assert ss_class(None) == NOT_ASSESSED
    assert ss_class("H") == "helix" and ss_class("E") == "strand" and ss_class("-") == "loop"


def test_indel_secondary_structure_not_assessed_without_model():
    record = _record(_PARENT, _PARTNER, five_last=30, three_first=26)
    placed = place_chimera(record, align_parents(record))
    report = junction_structure_report(record, placed, {"five": None, "three": None})
    assert report["indel_in_ss"] == NOT_ASSESSED
    assert report["junction_ss_5p"] == NOT_ASSESSED


def test_indel_secondary_structure_without_dssp_annotation_is_not_assessed():
    record = _record(_PARENT, _PARTNER, five_last=30, three_first=26)
    placed = place_chimera(record, align_parents(record))
    five = _entry("GENEA", _PARENT, ())
    three = _entry("GENEB", _PARTNER, ())
    mapped = {
        "five": map_structure(record.five, five, _structure(_PARENT, None))[0],
        "three": map_structure(record.three, three, _structure(_PARTNER, None))[0],
    }
    report = junction_structure_report(record, placed, mapped)
    assert report["indel_in_ss"] == NOT_ASSESSED


def test_indel_in_helix_is_reported():
    record = _record(_PARENT, _PARTNER, five_last=30, three_first=26)
    placed = place_chimera(record, align_parents(record))
    helix = "-" * 20 + "H" * 20 + "-" * (len(_PARTNER) - 40)
    five = _entry("GENEA", _PARENT, ())
    three = _entry("GENEB", _PARTNER, ())
    mapped = {
        "five": map_structure(record.five, five, _structure(_PARENT, helix))[0],
        "three": map_structure(record.three, three, _structure(_PARTNER, helix))[0],
    }
    report = junction_structure_report(record, placed, mapped)
    assert report["indel_ss"] == "helix"
    assert report["indel_in_ss"] is True


def test_clean_junction_reports_no_indel():
    record = _record(_PARENT, _PARTNER, five_last=30, three_first=31)
    placed = place_chimera(record, align_parents(record))
    report = junction_structure_report(record, placed, {"five": None, "three": None})
    assert report["indel_ss"] == "no_indel" and report["indel_in_ss"] is False


def test_map_structure_rejects_sequence_mismatch():
    record = _record(_PARENT, _PARTNER, five_last=30, three_first=31)
    entry = _entry("GENEA", _PARENT, ())
    mapped, status, detail = map_structure(record.five, entry, _structure(_PARTNER, "-" * 63))
    assert mapped is None and status == "afdb_sequence_mismatch" and "differs" in detail
    assert map_structure(record.five, None, None)[:2] == (None, "no_uniprot_entry")
    assert map_structure(record.five, entry, None)[:2] == (None, "no_afdb_model")
    assert map_structure(record.five, entry, None, "RuntimeError: 429")[:2] == (
        None,
        "afdb_fetch_failed",
    )


def test_schema_counts_broken_lost_and_inserted_contacts():
    # Parents differ at 14 (V/I), 22 (K/R), 38 (A/S) and 45 (S/A), among others.
    record = _record(_PARENT, _PARTNER, five_last=30, three_first=31)
    parents = align_parents(record)
    placed = place_chimera(record, parents)
    contacts = [(14, 38), (22, 45), (10, 12), (14, 22)]
    five = _entry("GENEA", _PARENT, ())
    three = _entry("GENEB", _PARTNER, ())
    mapped = {
        "five": map_structure(record.five, five, _structure(_PARENT, "-" * 63, contacts))[0],
        "three": None,
    }
    report = schema_report(record, parents, placed, mapped)
    # (14,38): chimera V(5') + S(3') -> parent pairs (V,A) and (I,S): broken.
    # (22,45): chimera K + A -> (K,S) / (R,A): broken. Others are within one block.
    assert report["schema_E_5p"] == 2
    assert report["schema_E_core_5p"] == 2
    assert report["schema_contacts_5p"] == 4
    assert report["schema_E_3p"] == NOT_ASSESSED

    deletion = _record(_PARENT, _PARTNER, five_last=30, three_first=38)
    parents = align_parents(deletion)
    placed = place_chimera(deletion, parents)
    mapped = {
        "five": map_structure(deletion.five, five, _structure(_PARENT, "-" * 63, [(14, 33)]))[0],
        "three": map_structure(deletion.three, three, _structure(_PARTNER, "-" * 63, [(27, 50)]))[
            0
        ],
    }
    report = schema_report(deletion, parents, placed, mapped)
    assert report["schema_lost_to_deletion_5p"] == 1 and report["schema_E_5p"] == 0

    insertion = _record(_PARENT, _PARTNER, five_last=30, three_first=26)
    parents = align_parents(insertion)
    placed = place_chimera(insertion, parents)
    mapped = {
        "five": None,
        "three": map_structure(insertion.three, three, _structure(_PARTNER, "-" * 63, [(27, 50)]))[
            0
        ],
    }
    report = schema_report(insertion, parents, placed, mapped)
    assert report["schema_touching_inserted_3p"] == 1


def test_low_plddt_contacts_are_excluded():
    record = _record(_PARENT, _PARTNER, five_last=30, three_first=31)
    parents = align_parents(record)
    placed = place_chimera(record, parents)
    structure = _structure(_PARENT, "-" * 63, [(14, 38)])
    low = ParentStructure(**{**structure.__dict__, "plddt": tuple([50.0] * 63)})
    entry = _entry("GENEA", _PARENT, ())
    mapped = {"five": map_structure(record.five, entry, low)[0], "three": None}
    report = schema_report(record, parents, placed, mapped)
    assert report["schema_contacts_5p"] == 0
    assert report["schema_plddt_coverage_5p"] == 0.0


# ------------------------------------------------------------------------------ tiers
_PASS = {
    "ga_domain": 20.0,
    "chimera_domain_bits": 80.0,
    "parent5_domain_bits": 85.0,
    "parent3_domain_bits": 84.0,
    "chimera_match_coverage": 0.95,
    "chimera_hit_spans_junction": True,
    "key_features_assessed": True,
    "key_point_lost": "",
    "key_point_altered": "",
    "indel_ss": "no_indel",
    "indel_in_ss": False,
    "junction_near_key": False,
    "key_annotation_status_5p": "annotated",
    "key_annotation_status_3p": "no_annotations_in_window",
    "key_unresolved": 0,
    "key_span_lost_or_altered": "",
    "key_motif_lost_or_altered": "",
    "structure_status_5p": "ok",
    "structure_status_3p": "ok",
}


def test_tier_intact_only_when_everything_assessed_and_passing():
    assert assign_tier(dict(_PASS))[0] == LIKELY_INTACT


@pytest.mark.parametrize(
    "missing",
    ["ga_domain", "parent5_domain_bits", "key_features_assessed", "indel_in_ss"],
)
def test_tier_missing_input_is_never_intact(missing):
    metrics = dict(_PASS)
    metrics[missing] = NOT_ASSESSED
    tier, rules, missing_inputs = assign_tier(metrics)
    assert tier != LIKELY_INTACT
    assert tier == NOT_ASSESSED and rules == ["N1_missing"] and missing_inputs


def test_tier_empty_metrics_is_not_assessed():
    assert assign_tier({})[0] == NOT_ASSESSED


def test_tier_below_ga_is_broken():
    tier, rules, _ = assign_tier({**_PASS, "chimera_domain_bits": 10.0})
    assert tier == LIKELY_BROKEN and "B1_below_GA" in rules
    tier, rules, _ = assign_tier({**_PASS, "chimera_domain_bits": None})
    assert tier == LIKELY_BROKEN


def test_tier_lost_active_site_is_broken_even_with_missing_hmm():
    tier, rules, _ = assign_tier({"key_point_lost": "A:Active site:33"})
    assert tier == LIKELY_BROKEN and rules == ["B3_point_feature_lost"]


def test_tier_parent_below_ga_is_uncertain():
    tier, rules, _ = assign_tier({**_PASS, "parent3_domain_bits": 15.0})
    assert tier == UNCERTAIN and "U1_parent_below_GA" in rules


@pytest.mark.parametrize(
    ("change", "rule"),
    [
        ({"chimera_domain_bits": 60.0}, "U2_bit_drop"),
        ({"chimera_match_coverage": 0.7}, "U3_partial_coverage"),
        ({"key_point_altered": "A:Site:3"}, "U4_point_feature_altered"),
        ({"indel_ss": "helix", "indel_in_ss": True}, "U5_indel_in_ss"),
        ({"indel_ss": "loop", "junction_near_key": True}, "U6_indel_near_key"),
        ({"chimera_hit_spans_junction": False}, "U7_junction_outside_hit"),
    ],
)
def test_tier_uncertain_rules(change, rule):
    tier, rules, _ = assign_tier({**_PASS, **change})
    assert tier == UNCERTAIN and rule in rules


def test_every_fired_rule_is_documented():
    for metrics in ({}, _PASS, {**_PASS, "chimera_domain_bits": 1.0}):
        _, rules, _ = assign_tier(metrics)
        assert all(rule in RULES for rule in rules)


def test_disagreements_are_reported():
    metrics = {**_PASS, "indel_in_ss": True, "junction_residues_in_inserts": 3}
    flags = disagreements(metrics, UNCERTAIN)
    assert "passes_GA_but_indel_in_helix_or_strand" in flags
    assert "passes_GA_with_junction_residues_in_insert_states" in flags
    below = {**_PASS, "chimera_domain_bits": 5.0, "schema_E_5p": 0, "schema_E_3p": 0}
    assert "below_GA_but_no_broken_contacts" in disagreements(below, LIKELY_BROKEN)


# ----------------------------------------------------------------- records and outputs
def test_record_roundtrip_and_row_aggregation():
    record = _record(_PARENT, _PARTNER, five_last=30, three_first=31)
    flags = {
        "paralog_sequence_flag": "pass",
        "genomic_similarity_flag": "pass",
        "breakpoint_repeat_flag": "pass",
        "sample_burden_flag": "pass",
        "read_evidence_flag": "pass",
        "normal_recurrence_flag": "pass",
    }
    add_row(record, {"cohort": "c1", "sample_id": "s1", "tumor_split_reads": "4", **flags})
    add_row(
        record,
        {
            "cohort": "c2",
            "sample_id": "s2",
            "tumor_split_reads": "",
            "read_evidence_flag": "flagged",
        },
    )
    from cfh.domain_function.record import finalize_support

    finalize_support(record)
    assert record.n_rows == 2 and record.cohorts == ["c1", "c2"]
    assert record.artefact_clean_rows == 1
    assert record.artefact_core_clean_rows == 1
    assert record.max_split_reads == 4
    assert record.artefact_flag_counts["paralog_sequence_flag"] == {"pass": 1, "not_assessed": 1}
    data = record_to_dict(record)
    assert data["chimeric_domain"] == record.chimeric_domain
    assert record_from_dict(data) == record


def test_assess_record_end_to_end_without_external_data(ets_hmm):
    record = _record(_PARENT, _PARTNER, five_last=30, three_first=31)
    metrics = assess_record(record, None, {"five": None, "three": None}, {})
    assert metrics["tier"] == NOT_ASSESSED
    assert metrics["tier_kind"] == "HEURISTIC"
    assert metrics["structure_5p"] == NOT_ASSESSED


def _shortlist_row(rid, fusion, tier, e, clean=1, core=None, reads=1, family="PF1"):
    return {
        "record_id": rid,
        "fusion": fusion,
        "pfam_family": family,
        "tier": tier,
        "schema_E_5p": e,
        "schema_E_3p": NOT_ASSESSED,
        "chimera_ga_margin": 10.0,
        "max_split_reads": reads,
        "n_samples": 1,
        "artefact_clean_rows": clean,
        "artefact_core_clean_rows": clean if core is None else core,
        "genomic_similarity_flag_counts": "pass:1" if clean else "flagged:1",
        "sample_burden_flag_counts": "pass:1",
    }


def test_shortlist_ranks_by_tier_then_e_and_keeps_named_cases():
    row = _shortlist_row
    rows = [
        row("a", "A::B", UNCERTAIN, 0),
        row("b", "C::D", LIKELY_INTACT, 5),
        row("c", "E::F", LIKELY_INTACT, 1),
        row("d", "FLI1::ERG", LIKELY_BROKEN, 9, clean=0, core=0),
        row("e", "G::H", LIKELY_INTACT, 0, clean=0, core=0),
    ]
    chosen = shortlist(rows, size=3)
    assert [r["record_id"] for r in chosen] == ["c", "b", "a", "d"]
    assert [r["shortlist_rank"] for r in chosen] == [1, 2, 3, 4]
    assert chosen[0]["shortlist_reason"] == "ranked"
    assert chosen[-1]["shortlist_reason"] == "named_case;flagged:genomic_similarity"


def test_shortlist_fills_with_core_clean_records_and_names_their_flags():
    row = _shortlist_row
    rows = [
        row("a", "A::B", UNCERTAIN, 0),
        row("f", "I::J", LIKELY_INTACT, 0, clean=0, core=2),
        row("g", "K::L", LIKELY_INTACT, 0, clean=0, core=0),
    ]
    chosen = shortlist(rows, size=5)
    assert [r["record_id"] for r in chosen] == ["a", "f"]
    assert chosen[1]["shortlist_reason"] == "fill;flagged:genomic_similarity"


def test_shortlist_uses_the_clean_record_of_a_group():
    row = _shortlist_row
    rows = [
        row("best_but_flagged", "A::B", LIKELY_INTACT, 0, clean=0, core=0),
        row("clean", "A::B", UNCERTAIN, 3),
    ]
    chosen = shortlist(rows, size=5)
    assert [r["record_id"] for r in chosen] == ["clean"]


def test_tier_counts_by_cohort_and_class():
    record = _record(_PARENT, _PARTNER, five_last=30, three_first=31)
    record.rows_per_cohort = {"c1": 3, "c2": 1}
    counts = tier_counts([{"record_id": "chim_test", "tier": UNCERTAIN}], [record])
    assert counts["records"]["c1|test"] == {UNCERTAIN: 1}
    assert counts["rows"]["c1|test"] == {UNCERTAIN: 3}


# ---------------------------------------------------------------------------- network
@pytest.mark.network
def test_fetch_real_pfam_hmm_and_release(tmp_path):
    from cfh.domain_function.pfam_hmm import fetch_hmm, fetch_pfam_release

    hmm = load_hmm(fetch_hmm("PF00178", tmp_path))
    assert hmm.accession.startswith(b"PF00178") if isinstance(hmm.accession, bytes) else True
    assert "pfam_version" in fetch_pfam_release()


@pytest.mark.network
def test_real_uniprot_and_afdb(tmp_path):
    from cfh.domain_function.afdb import AfdbClient
    from cfh.domain_function.uniprot import UniProtClient

    entry = UniProtClient(tmp_path / "u").entry("ERG")
    assert entry is not None and entry.accession == "P11308"
    assert any(f.type == "DNA binding" for f in entry.features)
    structure = AfdbClient(tmp_path / "a").structure("P11308")
    assert structure is not None and structure.ss is not None
    assert structure.sequence == entry.sequence
    assert ss_class(structure.ss_at(330)) in ("helix", "strand", "loop")
    assert not list((tmp_path / "a").glob("*.cif"))


def test_records_file_roundtrip_and_schema_check(tmp_path):
    from cfh.domain_function.record import read_records_json
    from cfh.domain_function.run import write_records

    record = _record(_PARENT, _PARTNER, five_last=30, three_first=31)
    path = tmp_path / "records.json.gz"
    write_records([record], path, {"code_sha": "x"})
    assert read_records_json(path) == [record]
    bad = tmp_path / "bad.json"
    bad.write_text('{"schema_version": "other", "records": []}')
    with pytest.raises(ValueError):
        read_records_json(bad)


def test_build_core_fails_closed_on_internal_stop_codon():
    from types import SimpleNamespace

    from cfh.domain_function.record import build_core
    from cfh.mapping.junction_translation import TranscriptExon, TranscriptModel

    def model(gene: str, cds: str) -> TranscriptModel:
        return TranscriptModel(
            gene,
            f"T{gene}",
            "GRCh38",
            1,
            (TranscriptExon(1, 100, 100 + len(cds) - 1),),
            100,
            100 + len(cds) - 1,
            cds,
        )

    models = {"A": model("A", "ATGTGAGCCTAA"), "B": model("B", "ATGGCCGCCTAA")}
    source = SimpleNamespace(
        resolve=lambda gene, annotated=None: gene,
        model=lambda gene: (models[gene], None),
    )
    row = {
        "five_prime_gene": "A",
        "three_prime_gene": "B",
        "five_prime_breakpoint": "105",
        "three_prime_breakpoint": "105",
        "genome_build": "GRCh38",
        "pfam_family": "PF1",
        "reconstitution_class": "clean_reconstitution",
    }
    record, problem = build_core(row, source)  # type: ignore[arg-type]
    assert record is None and problem is not None and "internal stop" in problem


def test_split_domain_uses_union_coverage_and_is_uncertain():
    metrics = {
        **_PASS,
        "chimera_match_coverage": 0.45,
        "chimera_match_coverage_union": 0.9,
        "chimera_hits_in_domain": 2,
        "parent5_hits_in_domain": 1,
        "parent3_hits_in_domain": 1,
    }
    tier, rules, _ = assign_tier(metrics)
    assert tier == UNCERTAIN and "U8_domain_split" in rules and "B2_low_coverage" not in rules
    tier, rules, _ = assign_tier({**metrics, "chimera_match_coverage_union": 0.4})
    assert tier == LIKELY_BROKEN and "B2_low_coverage" in rules


def test_union_match_coverage(ets_hmm):
    from cfh.domain_function.pfam_hmm import overlapping_hits, union_match_coverage

    assert union_match_coverage([]) is None
    window = FLI1[FLI1_ETS[0] - 31 : FLI1_ETS[1] + 30]
    result = scan(ets_hmm, [("fli1", window)])["fli1"]
    hits = overlapping_hits(result, 31, 111)
    assert len(hits) == 1 and union_match_coverage(hits) == hits[0].match_state_coverage
    assert overlapping_hits(None, 1, 10) == []


# ------------------------------------------- regression tests for the #155 review
def _make(
    five: str,
    three: str,
    *,
    five_last: int,
    three_first: int,
    five_domain: tuple[int, int] = (5, 60),
    three_domain: tuple[int, int] = (5, 60),
    hybrid: str = "",
    five_tx: str = "TXA",
    three_tx: str = "TXB",
) -> ChimeraRecord:
    protein = five[:five_last] + hybrid + three[three_first - 1 :]
    junction = five_last + len(hybrid)
    return ChimeraRecord(
        record_id="chim_test",
        pfam_family="PF00000",
        junction_class="test",
        class_reason="",
        five=ParentDomain("GENEA", five_tx, five, five_domain[0], five_domain[1], five_last),
        three=ParentDomain("GENEB", three_tx, three, three_domain[0], three_domain[1], three_first),
        hybrid_residue=hybrid,
        chimeric_protein=protein,
        domain_start=five_domain[0],
        domain_end=junction + (three_domain[1] - three_first + 1),
        junction_index=junction,
        pilot_indel_residues=None,
        pilot_indel_sequence="",
        pilot_alignment_parent="",
        pilot_alignment_identity=None,
    )


def _fli1_self_chimera() -> ChimeraRecord:
    """A FLI1/FLI1 junction inside the ETS domain: the chimera equals the parent."""
    record = _make(
        FLI1,
        FLI1,
        five_last=320,
        three_first=321,
        five_domain=FLI1_ETS,
        three_domain=FLI1_ETS,
    )
    record.five.gene = record.three.gene = "FLI1"
    return record


def _fli1_structure(sequence: str = FLI1) -> ParentStructure:
    n = len(sequence)
    return ParentStructure(
        accession="Q01543",
        model_id="AF-Q01543-F1",
        model_version="6",
        sequence=sequence,
        first=1,
        plddt=tuple([90.0] * n),
        ss="-" * n,
        contacts=(),
        rsa=tuple([0.5] * n),
    )


# Review item 1: missing or rejected AlphaFold DB data must never give likely_intact.
def test_review1_tier_requires_validated_structures():
    old_style = {k: v for k, v in _PASS.items() if not k.startswith("structure_status")}
    tier, _, missing = assign_tier({**old_style, "structure_5p": NOT_ASSESSED})
    assert tier == NOT_ASSESSED
    assert any(item.startswith("structure_5p") for item in missing)
    tier, _, missing = assign_tier({**_PASS, "structure_status_3p": "afdb_sequence_mismatch"})
    assert tier == NOT_ASSESSED and "structure_3p:afdb_sequence_mismatch" in missing


def test_review1_missing_and_mismatched_models_end_to_end(ets_hmm):
    record = _fli1_self_chimera()
    entry = UniProtEntry("FLI1", "Q01543", FLI1, ())
    entries = {"five": entry, "three": entry}
    rejected = {"five": None, "three": _fli1_structure(FLI1[:-1] + "W")}
    metrics = assess_record(record, ets_hmm, entries, rejected)
    assert metrics["chimera_passes_ga"] is True  # everything else would pass
    assert metrics["tier"] == NOT_ASSESSED
    assert "structure_5p:no_afdb_model" in metrics["tier_missing_inputs"]
    assert "structure_3p:afdb_sequence_mismatch" in metrics["tier_missing_inputs"]

    valid = {"five": _fli1_structure(), "three": _fli1_structure()}
    assert assess_record(record, ets_hmm, entries, valid)["tier"] == LIKELY_INTACT


# Review item 2: absent key residues must not be reported as retained.
_FIVE_BASE = _PARENT + _PARTNER[:40]
_THREE_BASE = _PARTNER + _PARENT[:40]
_FIVE_WITH_INSERT = _FIVE_BASE[:79] + "C" + _FIVE_BASE[79:]  # residue 80 has no 3' partner


def test_review2_parent_only_column_after_crossover_is_lost():
    record = _make(
        _FIVE_WITH_INSERT,
        _THREE_BASE,
        five_last=50,
        three_first=51,
        five_domain=(5, 95),
        three_domain=(5, 94),
    )
    parents = align_parents(record)
    assert parents.columns[parents.five_col[80]][1] is None  # gap in the 3' parent
    five = _entry("GENEA", _FIVE_WITH_INSERT, (KeyFeature("Active site", 80, 80, ""),))
    three = _entry("GENEB", _THREE_BASE, ())
    report = key_residue_report(
        record, parents, place_chimera(record, parents), {"five": five, "three": three}
    )
    assert "GENEA:Active site:80" not in report["key_retained_list"]
    assert report["key_lost"] == 1 and "GENEA:Active site:80" in report["key_point_lost"]
    assert assign_tier({**_PASS, **report})[0] == LIKELY_BROKEN


def test_review2_survival_is_decided_from_junction_coordinates():
    record = _make(_PARENT, _PARTNER, five_last=30, three_first=31)
    parents = align_parents(record)
    placed = place_chimera(record, parents)
    from cfh.domain_function.metrics import residue_status

    assert residue_status(record, parents, placed, "five", 14) == "retained"
    assert residue_status(record, parents, placed, "three", 38) == "retained"
    # 5' residue 38 (A) is past the crossover; the chimera has the 3' parent's S there.
    assert residue_status(record, parents, placed, "five", 38) == "swapped"
    # 5' residue 33 is past the crossover but both parents have the same residue.
    assert residue_status(record, parents, placed, "five", 33) == "conserved"


def test_review2_unaligned_uniprot_position_is_unresolved_not_retained():
    record = _make(_PARENT, _PARTNER, five_last=30, three_first=31)
    parents = align_parents(record)
    uniprot = _PARENT[:29] + "W" + _PARENT[29:]  # UniProt residue 30 has no transcript residue
    five = _entry("GENEA", uniprot, (KeyFeature("Active site", 30, 30, ""),))
    three = _entry("GENEB", _PARTNER, ())
    report = key_residue_report(
        record, parents, place_chimera(record, parents), {"five": five, "three": three}
    )
    # 334c050 dropped the unaligned position, leaving an "assessed" empty feature set.
    assert report["key_features_5p"] == 1
    assert report["key_unresolved"] == 1
    assert report["key_retained"] == 0
    tier, _, missing = assign_tier({**_PASS, **report})
    assert tier != LIKELY_INTACT and "key_residue_mapping" in missing


# Review item 3: unusable UniProt coordinates are "incomplete", not "nothing annotated".
_FUZZY_RAW = {
    "primaryAccession": "Q1",
    "sequence": {"value": _PARENT},
    "features": [
        {
            "type": "Active site",
            "location": {
                "start": {"value": None, "modifier": "UNKNOWN"},
                "end": {"value": 40, "modifier": "EXACT"},
            },
        },
        {
            "type": "Binding site",
            "location": {
                "start": {"value": 20, "modifier": "UNSURE"},
                "end": {"value": 20, "modifier": "UNSURE"},
            },
        },
    ],
}


def test_review3_fuzzy_features_are_kept_as_incomplete():
    entry = parse_entry("GENEA", _FUZZY_RAW)
    assert entry.features == ()
    assert len(entry.incomplete) == 2
    assert entry.incomplete[0].may_overlap(30, 50) and not entry.incomplete[0].may_overlap(41, 50)


def test_review3_incomplete_annotation_is_not_assessed():
    record = _make(_PARENT, _PARTNER, five_last=30, three_first=31)
    parents = align_parents(record)
    placed = place_chimera(record, parents)
    incomplete = parse_entry("GENEA", _FUZZY_RAW)
    empty = _entry("GENEB", _PARTNER, ())
    report = key_residue_report(record, parents, placed, {"five": incomplete, "three": empty})
    assert report["key_features_assessed"] is False
    assert report["key_annotation_status_5p"] == "incomplete_annotation"
    assert report["key_annotation_status_3p"] == "no_annotations_in_window"
    assert report["key_features_assessed"] is False
    tier, _, missing = assign_tier({**_PASS, **report})
    assert tier == NOT_ASSESSED and "key_annotation_5p:incomplete_annotation" in missing

    clean = key_residue_report(
        record, parents, placed, {"five": _entry("GENEA", _PARENT, ()), "three": empty}
    )
    assert clean["key_features_assessed"] is True
    assert clean["key_annotation_status_5p"] == "no_annotations_in_window"


# Review item 4: count each biological sample once; never sum one sample's rows.
def test_review4_samples_and_reads_are_counted_once(monkeypatch):
    import cfh.domain_function.record as record_module

    core = _record(_PARENT, _PARTNER, five_last=30, three_first=31)
    monkeypatch.setattr(record_module, "build_core", lambda row, source: (core, None))

    def row(cohort, sample, split, five_bp="1"):
        return {
            "cohort": cohort,
            "sample_id": sample,
            "tumor_split_reads": str(split),
            "reconstitution_class": "clean_reconstitution",
            "genome_build": "GRCh38",
            "five_prime_breakpoint": five_bp,
        }

    rows = [
        row("target_consolidated", "TARGET-30-PASEGA-01A", 2),
        row("target_starfusion", "TARGET-30-PASEGA-01A", 2),
        row("depmap_22q2", "ACH-001552", 48),
        row("depmap_22q2", "ACH-001552", 3, five_bp="2"),  # second event row, same sample
        row("depmap_22q2", "ACH-000555", 5),
    ]
    records, _ = record_module.build_records(rows, {"GRCh38": object()})  # type: ignore[dict-item]
    (record,) = records
    assert record.n_rows == 5
    assert record.n_samples == 3
    assert record.max_split_reads == 48
    assert record.total_split_reads == 2 + 48 + 5
    assert record.samples_per_cohort == {
        "depmap_22q2": 2,
        "target_consolidated": 1,
        "target_starfusion": 1,
    }
    assert {o["biological_sample"] for o in record.observations} == {
        "TARGET:TARGET-30-PASEGA-01A",
        "depmap_22q2:ACH-001552",
        "depmap_22q2:ACH-000555",
    }
    assert record_from_dict(record_to_dict(record)) == record


# Review item 5: the same sequence from different origins stays separate.
def test_review5_record_identity_includes_origin():
    from cfh.domain_function.record import record_key, sequence_key

    # Parents are identical at residues 29-37, so both crossovers give one sequence.
    first = _record(_PARENT, _PARTNER, five_last=30, three_first=31)
    second = _record(_PARENT, _PARTNER, five_last=31, three_first=32)
    assert first.chimeric_domain_flanked == second.chimeric_domain_flanked
    assert sequence_key(first) == sequence_key(second)
    assert record_key(first) != record_key(second)
    moved_domain = _make(_PARENT, _PARTNER, five_last=30, three_first=31, five_domain=(6, 60))
    assert record_key(moved_domain) != record_key(first)


def test_review5_build_records_keeps_each_origin(monkeypatch):
    import cfh.domain_function.record as record_module
    from cfh.domain_function.record import record_key

    cores = {}
    for five_last in (30, 31):
        core = _record(_PARENT, _PARTNER, five_last=five_last, three_first=five_last + 1)
        core.record_id = record_key(core)
        cores[str(five_last)] = core
    monkeypatch.setattr(
        record_module, "build_core", lambda row, source: (cores[row["five_prime_breakpoint"]], None)
    )
    rows = [
        {
            "cohort": "c",
            "sample_id": f"s{bp}",
            "five_prime_breakpoint": bp,
            "reconstitution_class": "clean_reconstitution",
            "genome_build": "GRCh38",
        }
        for bp in ("30", "31")
    ]
    records, _ = record_module.build_records(rows, {"GRCh38": object()})  # type: ignore[dict-item]
    assert len(records) == 2
    assert {r.five.junction_residue for r in records} == {30, 31}


# Review item 6: the records cache is reused only when its fingerprint matches.
def test_review6_records_cache_is_fingerprinted(tmp_path):
    import pickle

    from cfh.domain_function.run import load_or_build_records

    builds = []

    def build():
        builds.append(1)
        return [len(builds)], Counter({"p": len(builds)})

    cache = tmp_path / "records.pickle"
    cache.write_bytes(pickle.dumps(([0], Counter())))  # pre-fingerprint cache format
    fingerprint = {"screen_sha256": "a"}
    records, _, _, reused = load_or_build_records(cache, lambda: fingerprint, build)
    assert (records, reused) == ([1], False)
    records, _, _, reused = load_or_build_records(cache, lambda: fingerprint, build)
    assert (records, reused) == ([1], True)
    records, _, _, reused = load_or_build_records(cache, lambda: {"screen_sha256": "b"}, build)
    assert (records, reused) == ([2], False)


def test_review6_fingerprint_tracks_screen_and_transcript_cache(tmp_path):
    from cfh.domain_function.run import records_fingerprint

    screen = tmp_path / "screen.tsv.gz"
    screen.write_bytes(gzip.compress(b"a\n"))
    transcripts = tmp_path / "tx"
    transcripts.mkdir()
    (transcripts / "GENE.json").write_text("{}")
    base = records_fingerprint(screen, {"GRCh38": transcripts})
    assert records_fingerprint(screen, {"GRCh38": transcripts}) == base
    (transcripts / "OTHER.json").write_text("{}")
    assert records_fingerprint(screen, {"GRCh38": transcripts}) != base
    after_cache = records_fingerprint(screen, {"GRCh38": transcripts})
    screen.write_bytes(gzip.compress(b"b\n"))
    assert records_fingerprint(screen, {"GRCh38": transcripts}) != after_cache
    assert set(base["code"]) >= {"cfh.domain_function.record", "cfh.mapping.junction_translation"}


# Non-blocking review item: lost/altered span features block likely_intact.
def test_review_span_feature_lost_blocks_likely_intact():
    record = _make(_PARENT, _PARTNER, five_last=30, three_first=38)  # 7 deleted columns
    parents = align_parents(record)
    five = _entry("GENEA", _PARENT, (KeyFeature("Motif", 31, 36, "catalytic"),))
    three = _entry("GENEB", _PARTNER, ())
    report = key_residue_report(
        record, parents, place_chimera(record, parents), {"five": five, "three": three}
    )
    tier, rules, _ = assign_tier({**_PASS, **report})
    assert tier == UNCERTAIN and "U9_motif_lost_or_altered" in rules
    assert report["key_span_lost_or_altered"] == "GENEA:Motif:31-36(catalytic)"
    assert report["key_motif_lost_or_altered"] == "GENEA:Motif:31-36(catalytic)"


def test_review_long_region_lost_is_reported_not_tiered():
    record = _make(_PARENT, _PARTNER, five_last=30, three_first=38)  # 7 deleted columns
    parents = align_parents(record)
    five = _entry("GENEA", _PARENT, (KeyFeature("Region", 20, 50, "interaction"),))
    three = _entry("GENEB", _PARTNER, ())
    report = key_residue_report(
        record, parents, place_chimera(record, parents), {"five": five, "three": three}
    )
    assert report["key_span_lost_or_altered"] == "GENEA:Region:20-50(interaction)"
    assert report["key_motif_lost_or_altered"] == ""
    assert assign_tier({**_PASS, **report})[0] == LIKELY_INTACT


def _synthetic_model(gene: str, protein: str, start: int):
    from cfh.mapping.junction_translation import TranscriptDomain, TranscriptExon, TranscriptModel

    codon = {
        "A": "GCC", "C": "TGC", "D": "GAC", "E": "GAG", "F": "TTC", "G": "GGC", "H": "CAC",
        "I": "ATC", "K": "AAG", "L": "CTG", "M": "ATG", "N": "AAC", "P": "CCC", "Q": "CAG",
        "R": "CGC", "S": "AGC", "T": "ACC", "V": "GTG", "W": "TGG", "Y": "TAC",
    }  # fmt: skip
    cds = "".join(codon[aa] for aa in protein) + "TAA"
    end = start + len(cds) - 1
    return TranscriptModel(
        gene,
        f"ENST_{gene}",
        "GRCh38",
        1,
        (TranscriptExon(1, start, end),),
        start,
        end,
        cds,
        (TranscriptDomain("PF00001", 5, 60),),
    )


def test_review5_real_rebuild_keeps_same_sequence_from_two_crossovers():
    """Two junctions give one chimeric sequence (parents identical at 29-37)."""
    from types import SimpleNamespace

    import cfh.domain_function.record as record_module

    models = {"GENEA": _synthetic_model("GENEA", _PARENT, 1000)}
    models["GENEB"] = _synthetic_model("GENEB", _PARTNER, 5000)
    source = SimpleNamespace(
        resolve=lambda gene, annotated=None: gene, model=lambda gene: (models[gene], None)
    )

    def row(after_residue: int) -> dict[str, str]:
        # both_breakpoint_bases_retained: 5' keeps through the bp, 3' starts at it.
        return {
            "cohort": "c",
            "sample_id": f"s{after_residue}",
            "five_prime_gene": "GENEA",
            "three_prime_gene": "GENEB",
            "five_prime_breakpoint": str(1000 + 3 * after_residue - 1),
            "three_prime_breakpoint": str(5000 + 3 * after_residue),
            "genome_build": "GRCh38",
            "pfam_family": "PF00001",
            "reconstitution_class": "clean_reconstitution",
        }

    rows = [row(30), row(31)]
    cores = [record_module.build_core(r, source)[0] for r in rows]  # type: ignore[arg-type]
    assert all(core is not None for core in cores)
    assert cores[0].chimeric_domain_flanked == cores[1].chimeric_domain_flanked
    records, problems = record_module.build_records(rows, {"GRCh38": source})  # type: ignore[dict-item]
    assert not problems
    assert len(records) == 2
    assert sorted(r.five.junction_residue for r in records) == [30, 31]
    assert len({getattr(r, "sequence_id", None) for r in records}) == 1
