import json
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from cfh.cohort.domain_screen import (
    SCREEN_COLUMNS,
    is_screenable_fusion,
    screen_cohort,
    write_screen_outputs,
)
from cfh.ingestion.sv_parser import parse_sv_file
from cfh.mapping.junction_domain_classifier import (
    CLEAN_RECONSTITUTION,
    GAP_OR_DELETION,
    NO_TRANSCRIPT_MODEL,
    NOT_SPLIT_ON_BOTH,
    OUT_OF_FRAME,
    PREMATURE_STOP_AT_JUNCTION,
    RECONSTITUTION_CLASSES,
    SPLICE_SITE_UNRESOLVED,
    TANDEM_DUPLICATION_OR_INSERTION,
    align_global,
    classify_family,
)
from cfh.mapping.junction_event import (
    MODE_KEY_DOMAIN,
    NO_SHARED_FAMILY,
    assess_event,
    junction_call_columns,
    read_support_columns,
)
from cfh.mapping.junction_translation import (
    ACCEPTOR_SITE_NT,
    CONVENTION_THREE_PRIME_BASE_DROPPED,
    DOMAIN_INTACT,
    DOMAIN_LOST,
    DOMAIN_TRUNCATED,
    DONOR_SITE_NT,
    FRAME_IN,
    FRAME_OUT,
    FRAME_UNKNOWN,
    MECHANISM_MID_EXON,
    MECHANISM_NONE,
    MECHANISM_SPLICE,
    MECHANISM_SPLICE_SITE_UNRESOLVED,
    REGION_CODING_EXON,
    REGION_INTRON,
    TranscriptDomain,
    TranscriptExon,
    TranscriptModel,
    locate_breakpoint,
    normalize_genome_build,
    predict_junction,
    translate,
)
from cfh.mapping.transcript_models import (
    CANONICAL_ANNOTATION_UNMATCHED,
    CANONICAL_NO_ANNOTATION,
    SOURCE_ANNOTATED,
    EnsemblCdsClient,
    TranscriptModelSource,
    genome_build_of_genome_nexus,
    source_annotated_transcripts,
)

FIXTURES = Path(__file__).parent / "fixtures"
TRANSCRIPTS = json.loads((FIXTURES / "junction" / "transcripts.json").read_text())
MSK_ARCHIVE_FIXTURE = FIXTURES / "sv" / "msk_archive_sv.txt"
BUILD = "GRCh37"
ETS = "PF00178"
ANK = "PF12796"
ARF = "PF07392"

FLI1_ERG_ANNOTATION = (
    "FLI1 (NM_002017) - ERG (NM_182918) rearrangement: "
    "t(11;21)(q24.3;q22.2)(chr11:g.128680564::chr21:g.39755655)"
)
GABPA_ERG_ANNOTATION = (
    "GABPA (NM_002040) - ERG (NM_182918) Rearrangement :  c.1137-1426:GABPA_c.958:ERGinv"
)
CDKN2_ANNOTATION = (
    "CDKN2B (NM_004936) - CDKN2A (NM_058195) rearrangement: c.338:CDKN2B_c.303:CDKN2Adel"
)
SPLIT_SITES = {"tumor_split_reads": None}


class FakeGenomeNexus:
    """Serves the committed fixture, with an optional canonical override."""

    base_url = "https://example.invalid"
    session = None

    def __init__(self, canonical_override=None):
        self.canonical_ids = {**TRANSCRIPTS["canonical"], **(canonical_override or {})}
        self.isoform_calls = []

    def _payloads(self, gene):
        return TRANSCRIPTS["transcripts"][gene]

    def fetch_canonical_transcript(self, gene):
        wanted = self.canonical_ids[gene]
        for payload in self._payloads(gene):
            if payload["transcriptId"] == wanted:
                return payload
        # GABPA's real canonical isoform is not in the fixture: a stand-in with
        # that id, so "differs from canonical" is exercised on real ids.
        return {**self._payloads(gene)[0], "transcriptId": wanted, "refseqMrnaId": None}

    def fetch_transcripts_for_gene(self, gene):
        self.isoform_calls.append(gene)
        return self._payloads(gene)


class FakeCds:
    def __init__(self):
        self.calls = []

    def fetch_cds(self, transcript_id):
        self.calls.append(transcript_id)
        for payloads in TRANSCRIPTS["transcripts"].values():
            for payload in payloads:
                if payload["transcriptId"] == transcript_id:
                    return payload["cds"]
        return None


def _source(canonical_override=None):
    return TranscriptModelSource(
        BUILD, genome_nexus_client=FakeGenomeNexus(canonical_override), cds_client=FakeCds()
    )


def _resolved_model(gene, transcript_id, source=None):
    source = source or _source()
    for payload in TRANSCRIPTS["transcripts"][gene]:
        if payload["transcriptId"] == transcript_id:
            resolved = source.resolve(gene, [payload["refseqMrnaId"]])
            assert resolved.transcript_id == transcript_id
            model, problem = source.model(resolved)
            assert problem is None
            return model
    raise AssertionError(transcript_id)


def _fli1_erg(**kwargs):
    return predict_junction(
        "FLI1",
        128680564,
        _resolved_model("FLI1", "ENST00000527786"),
        "ERG",
        39755655,
        _resolved_model("ERG", "ENST00000288319"),
        BUILD,
        **kwargs,
    )


def test_fli1_erg_is_a_tandem_duplication_with_the_published_peptide():
    fli1 = _resolved_model("FLI1", "ENST00000527786")
    erg = _resolved_model("ERG", "ENST00000288319")
    prediction = predict_junction("FLI1", 128680564, fli1, "ERG", 39755655, erg, BUILD)

    assert prediction.frame_status == FRAME_IN
    assert prediction.mechanism == MECHANISM_MID_EXON
    assert "LSRALRYYYDKNI|YYYDKNIMTKV" in prediction.junction_peptide
    assert prediction.five_prime_location.region == REGION_CODING_EXON
    assert prediction.five_prime_cds_retained == 1040
    assert prediction.three_prime_cds_start == 1110

    call = classify_family(prediction, fli1, erg, ETS)
    assert call.classification == TANDEM_DUPLICATION_OR_INSERTION
    assert call.junction_indel.length == 7
    assert call.junction_indel.residues == "YYYDKNI"
    assert call.junction_indel.is_tandem
    assert call.alignment_identity == 1.0


def test_fli1_erg_alternative_convention_frameshifts():
    prediction = _fli1_erg(coordinate_convention=CONVENTION_THREE_PRIME_BASE_DROPPED)
    assert prediction.frame_status == FRAME_OUT
    assert prediction.coordinate_convention == CONVENTION_THREE_PRIME_BASE_DROPPED


def test_gabpa_erg_has_no_transcript_model_and_says_why():
    source = _source()
    rows = assess_event(
        five_prime_gene="GABPA",
        five_prime_breakpoint=27139889,
        three_prime_gene="ERG",
        three_prime_breakpoint=39755807,
        annotation=GABPA_ERG_ANNOTATION,
        genome_build=BUILD,
        source=source,
    )
    ets = next(row for row in rows if row["pfam_family"] == ETS)
    assert ets["reconstitution_class"] == NO_TRANSCRIPT_MODEL
    assert ets["junction_mechanism"] == MECHANISM_NONE
    assert ets["translated_frame_status"] == FRAME_UNKNOWN
    assert "no downstream splice acceptor" in ets["class_reason"]
    assert ets["junction_peptide"] is None
    # The intronic break is reported, not silently snapped.
    assert "2791 nt past exon 9" in ets["five_prime_break_location"]
    assert "exon 9's splice donor" in ets["junction_notes"]
    assert not any(
        row["reconstitution_class"]
        in (CLEAN_RECONSTITUTION, TANDEM_DUPLICATION_OR_INSERTION, GAP_OR_DELETION)
        for row in rows
    )


def test_cdkn2b_cdkn2a_uses_the_annotated_transcript_and_positions_differ():
    source = _source(canonical_override={"CDKN2A": "ENST00000304494"})
    annotated = source.resolve("CDKN2A", source_annotated_transcripts(CDKN2_ANNOTATION, "CDKN2A"))
    assert annotated.source == SOURCE_ANNOTATED
    assert annotated.refseq_id == "NM_058195"
    assert annotated.transcript_id == "ENST00000579755"  # CCDS-bearing of the two NM_058195 hits
    assert annotated.differs_from_canonical is True
    assert annotated.canonical_transcript_id == "ENST00000304494"
    assert annotated.note == (
        "NM_058195 maps to ENST00000361570, ENST00000579755; used ENST00000579755"
    )

    position = 21971098
    p14, _ = source.model(annotated)
    p16, _ = source.model(source.resolve("CDKN2A", ["NM_000077"]))
    assert locate_breakpoint(p14, position).cds_position == 303  # aa 101
    assert locate_breakpoint(p16, position).cds_position == 260  # aa 87

    cdkn2b = _resolved_model("CDKN2B", "ENST00000276925", source)
    on_p14 = predict_junction("CDKN2B", 22006065, cdkn2b, "CDKN2A", position, p14, BUILD)
    on_p16 = predict_junction("CDKN2B", 22006065, cdkn2b, "CDKN2A", position, p16, BUILD)
    p14_spans = {s.family: s for s in on_p14.domain_spans if s.partner == "three_prime"}
    p16_spans = {s.family: s for s in on_p16.domain_spans if s.partner == "three_prime"}
    assert set(p14_spans) == {ARF}
    assert p14_spans[ARF].status == DOMAIN_LOST
    assert p16_spans[ANK].status == DOMAIN_TRUNCATED
    assert p16_spans[ANK].retained_start_aa == 88
    assert on_p14.frame_status == FRAME_IN and on_p16.frame_status == FRAME_OUT


def test_cdkn2b_cdkn2a_event_is_evaluated_on_the_annotated_isoform():
    source = _source()
    rows = assess_event(
        five_prime_gene="CDKN2B",
        five_prime_breakpoint=22006065,
        three_prime_gene="CDKN2A",
        three_prime_breakpoint=21971098,
        annotation=CDKN2_ANNOTATION,
        genome_build=BUILD,
        source=source,
    )
    assert [row["reconstitution_class"] for row in rows] == [NO_SHARED_FAMILY]
    assert rows[0]["three_prime_transcript"] == "NM_058195"
    assert rows[0]["three_prime_transcript_source"] == SOURCE_ANNOTATED
    assert rows[0]["three_prime_transcript_differs_from_canonical"] is False


def test_annotation_transcript_parsing():
    assert source_annotated_transcripts(FLI1_ERG_ANNOTATION, "ERG") == ["NM_182918"]
    assert source_annotated_transcripts("CDKN2A ( NM_000077) x", "CDKN2A") == ["NM_000077"]
    assert source_annotated_transcripts("CDKN2A (NM_000077, NM_058195) x", "CDKN2A") == [
        "NM_000077",
        "NM_058195",
    ]
    assert source_annotated_transcripts("CDKN2Ap14ARF (NM_058195.3) x", "CDKN2A") == ["NM_058195"]
    assert source_annotated_transcripts("ERGIC1 (NM_001) x", "ERG") == []
    assert source_annotated_transcripts(None, "ERG") == []


# --- Synthetic transcripts -------------------------------------------------

_CODON = {
    "A": "GCC",
    "C": "TGC",
    "D": "GAC",
    "E": "GAG",
    "F": "TTC",
    "G": "GGC",
    "H": "CAC",
    "I": "ATC",
    "K": "AAG",
    "L": "CTG",
    "M": "ATG",
    "N": "AAC",
    "P": "CCC",
    "Q": "CAG",
    "R": "CGC",
    "S": "AGC",
    "T": "ACC",
    "V": "GTC",
    "W": "TGG",
    "Y": "TAC",
}
DOMAIN = "ACDEFGHIKLMNPQRSTVWYACDEFGHIKL"  # 30 residues
FAMILY = "PF99999"


def _cds(protein):
    return "".join(_CODON[residue] for residue in protein) + "TAA"


def _model(
    gene, protein, domain_start, domain_end, *, exons=None, strand=1, family=FAMILY, domains=None
):
    cds = _cds(protein)
    if exons is None:
        exons = (TranscriptExon(1, 1000, 1000 + len(cds) - 1),)
    low = min(e.start for e in exons)
    high = max(e.end for e in exons)
    return TranscriptModel(
        gene=gene,
        transcript_id=f"T_{gene}",
        genome_build=BUILD,
        strand=strand,
        exons=tuple(exons),
        cds_start=low,
        cds_end=high,
        cds_sequence=cds,
        domains=domains
        if domains is not None
        else (TranscriptDomain(family, domain_start, domain_end),),
    )


def _genomic_at_cds(model, cds_position):
    """Single-exon + strand helper: genomic position of a CDS base."""
    return model.exons[0].start + cds_position - 1


def _split_case(five_keep_residues, three_start_residue, three_offset_bases=0):
    """Both partners carry DOMAIN at residues 11-40 of a 60-residue protein."""
    flank = "M" + "S" * 9
    protein = flank + DOMAIN + "T" * 20
    five = _model("AAA", protein, 11, 40)
    three = _model("BBB", protein, 11, 40)
    five_break = _genomic_at_cds(five, five_keep_residues * 3)
    three_break = _genomic_at_cds(three, (three_start_residue - 1) * 3 + 1 + three_offset_bases)
    return five, five_break, three, three_break


def _predict(five, five_break, three, three_break, build=BUILD, **kwargs):
    return predict_junction("AAA", five_break, five, "BBB", three_break, three, build, **kwargs)


def test_clean_reconstitution_rebuilds_the_whole_domain():
    five, five_break, three, three_break = _split_case(25, 26)
    prediction = _predict(five, five_break, three, three_break)
    assert prediction.frame_status == FRAME_IN
    assert prediction.mechanism == MECHANISM_MID_EXON
    call = classify_family(prediction, five, three, FAMILY)
    assert call.classification == CLEAN_RECONSTITUTION
    assert call.chimeric_domain == DOMAIN
    assert call.junction_indel is None


def test_tandem_duplication_when_the_three_prime_side_restarts_early():
    five, five_break, three, three_break = _split_case(25, 19)
    call = classify_family(_predict(five, five_break, three, three_break), five, three, FAMILY)
    assert call.classification == TANDEM_DUPLICATION_OR_INSERTION
    assert call.junction_indel.length == 7
    assert call.junction_indel.is_tandem


def test_gap_when_residues_are_skipped():
    five, five_break, three, three_break = _split_case(20, 26)
    call = classify_family(_predict(five, five_break, three, three_break), five, three, FAMILY)
    assert call.classification == GAP_OR_DELETION
    assert call.junction_indel.length == -5


REPEAT_DOMAIN = "M" + "A" * 30 + "W"  # 32 residues: a low-complexity run


def _repeat_case(five_last_domain_residue, three_first_domain_residue):
    """Identical parents carrying REPEAT_DOMAIN at residues 11-42."""
    protein = "M" + "S" * 9 + REPEAT_DOMAIN + "T" * 20
    five = _model("AAA", protein, 11, 42)
    three = _model("BBB", protein, 11, 42)
    five_break = _genomic_at_cds(five, (10 + five_last_domain_residue) * 3)
    three_break = _genomic_at_cds(three, (10 + three_first_domain_residue - 1) * 3 + 1)
    return classify_family(_predict(five, five_break, three, three_break), five, three, FAMILY)


def test_duplication_inside_a_repeat_is_not_hidden_by_gap_placement():
    # Join after domain residue 21, restart at residue 12: ten extra alanines.
    # The aligner may put that gap anywhere along the run, far from the junction.
    call = _repeat_case(21, 12)
    assert len(call.chimeric_domain) == len(REPEAT_DOMAIN) + 10
    assert call.classification == TANDEM_DUPLICATION_OR_INSERTION
    assert call.junction_indel.length == 10
    assert call.junction_indel.residues == "A" * 10
    assert call.junction_indel.is_tandem


def test_deletion_inside_a_repeat_is_not_hidden_by_gap_placement():
    call = _repeat_case(21, 30)  # skips domain residues 22-29
    assert call.classification == GAP_OR_DELETION
    assert call.junction_indel.length == -8


def test_repeat_calls_do_not_depend_on_where_in_the_run_the_join_falls():
    for five_last in range(5, 25):
        assert _repeat_case(five_last, five_last + 1).classification == CLEAN_RECONSTITUTION
        dup = _repeat_case(five_last, five_last - 3)
        assert dup.classification == TANDEM_DUPLICATION_OR_INSERTION
        assert dup.junction_indel.length == 4


def test_clean_is_never_called_when_the_chimera_length_differs():
    for five_last, three_first in ((21, 12), (21, 30), (15, 14), (15, 18)):
        call = _repeat_case(five_last, three_first)
        assert len(call.chimeric_domain) != len(REPEAT_DOMAIN)
        assert call.classification != CLEAN_RECONSTITUTION


def test_paralog_indel_away_from_the_junction_does_not_count():
    # The 3' paralog has two extra residues near its C-terminus; the junction
    # itself is in register, so the call is clean.
    five_protein = "M" + "S" * 9 + DOMAIN + "T" * 20
    three_domain = DOMAIN[:27] + "GG" + DOMAIN[27:]
    three_protein = "M" + "S" * 9 + three_domain + "T" * 20
    five = _model("AAA", five_protein, 11, 40)
    three = _model("BBB", three_protein, 11, 42)
    five_break = _genomic_at_cds(five, 25 * 3)
    three_break = _genomic_at_cds(three, 25 * 3 + 1)
    call = classify_family(_predict(five, five_break, three, three_break), five, three, FAMILY)
    assert call.classification == CLEAN_RECONSTITUTION
    assert call.junction_indel is None


def test_out_of_frame_junction():
    five, five_break, three, three_break = _split_case(25, 26, three_offset_bases=1)
    prediction = _predict(five, five_break, three, three_break)
    assert prediction.frame_status == FRAME_OUT
    assert classify_family(prediction, five, three, FAMILY).classification == OUT_OF_FRAME


def test_domain_not_cut_on_one_side():
    five, five_break, three, three_break = _split_case(40, 26)  # 5' keeps it whole
    prediction = _predict(five, five_break, three, three_break)
    call = classify_family(prediction, five, three, FAMILY)
    assert call.classification == NOT_SPLIT_ON_BOTH
    assert "AAA: intact" in call.reason and "BBB: truncated" in call.reason
    spans = {span.partner: span.status for span in prediction.domain_spans}
    assert spans["five_prime"] == DOMAIN_INTACT


def test_junction_stop_codon_is_called_out():
    # The 5' side ends one base into a tyrosine codon (T of TAC), and the 3'
    # side resumes at the second base of a lysine codon (AAG), forming TAG.
    five_protein = "M" + "S" * 9 + DOMAIN[:9] + "Y" + DOMAIN[10:] + "T" * 20
    five = _model("AAA", five_protein, 11, 40)
    three_protein = "M" + "S" * 9 + DOMAIN[:9] + "K" + DOMAIN[10:] + "T" * 20
    three = _model("BBB", three_protein, 11, 40)
    five_break = _genomic_at_cds(five, 19 * 3 + 1)
    three_break = _genomic_at_cds(three, 19 * 3 + 2)
    prediction = _predict(five, five_break, three, three_break)
    assert prediction.frame_status == FRAME_IN
    assert prediction.hybrid_codon == "TAG"
    assert prediction.premature_stop_at_junction
    call = classify_family(prediction, five, three, FAMILY)
    assert call.classification == PREMATURE_STOP_AT_JUNCTION


def _two_exon(gene, protein, exon_one_end_cds, strand=1, intron=1000):
    """Two coding exons separated by ``intron`` nt; low coordinates first."""
    cds_length = (len(protein) + 1) * 3
    if strand == 1:
        first = TranscriptExon(1, 1000, 1000 + exon_one_end_cds - 1)
        second = TranscriptExon(
            2, first.end + intron + 1, first.end + intron + cds_length - exon_one_end_cds
        )
        return _model(gene, protein, 11, 40, exons=(first, second))
    second = TranscriptExon(2, 1000, 1000 + cds_length - exon_one_end_cds - 1)
    first = TranscriptExon(1, second.end + intron + 1, second.end + intron + exon_one_end_cds)
    return _model(gene, protein, 11, 40, exons=(first, second), strand=-1)


_PROTEIN = "M" + "S" * 9 + DOMAIN + "T" * 20


def test_two_intronic_breaks_are_splice_mediated_and_recorded():
    # Exon 1 ends at residue 25, exon 2 of the 3' gene starts at residue 26.
    five = _two_exon("AAA", _PROTEIN, 75)
    three = _two_exon("BBB", _PROTEIN, 75)
    prediction = _predict(five, 1000 + 75 + 300, three, 1000 + 75 + 400)
    assert prediction.mechanism == MECHANISM_SPLICE
    assert prediction.five_prime_location.region == REGION_INTRON
    assert prediction.five_prime_cds_retained == 75
    assert prediction.three_prime_cds_start == 76
    assert prediction.five_prime_location.distance_from_upstream_exon_nt == 301
    assert any("exon 1's splice donor" in note for note in prediction.notes)
    assert any("exon 2's splice acceptor" in note for note in prediction.notes)
    call = classify_family(prediction, five, three, FAMILY)
    assert call.classification == CLEAN_RECONSTITUTION


def test_exon_boundary_breaks_are_an_exonic_dna_join_without_snapping_notes():
    five = _two_exon("AAA", _PROTEIN, 75)
    three = _two_exon("BBB", _PROTEIN, 75)
    five_last_exon_one = 1000 + 75 - 1
    three_first_exon_two = three.exons[1].start
    prediction = _predict(five, five_last_exon_one, three, three_first_exon_two)
    assert prediction.mechanism == MECHANISM_MID_EXON
    assert prediction.five_prime_cds_retained == 75 and prediction.three_prime_cds_start == 76
    assert prediction.notes == ()


def _exact_exon_edge_join(strand):
    """5' exon 1 is 9 nt; the 3' break is on the first base of 3' exon 2 (c.10)."""
    five = _two_exon("AAA", _PROTEIN, 9, strand=strand)
    three = _two_exon("BBB", _PROTEIN, 9, strand=strand)
    five_last = five.exons[0].end if strand == 1 else five.exons[0].start
    three_first = three.exons[1].start if strand == 1 else three.exons[1].end
    assert locate_breakpoint(five, five_last).at_exon_edge == "last_base"
    assert locate_breakpoint(three, three_first).cds_position == 10
    return five, five_last, three, three_first


@pytest.mark.parametrize("strand", [1, -1])
def test_exact_exon_edge_join_applies_both_coordinate_conventions(strand):
    five, five_last, three, three_first = _exact_exon_edge_join(strand)
    kept = _predict(five, five_last, three, three_first)
    assert kept.mechanism == MECHANISM_MID_EXON
    assert kept.frame_status == FRAME_IN
    assert kept.five_prime_cds_retained == 9 and kept.three_prime_cds_start == 10

    dropped = _predict(
        five,
        five_last,
        three,
        three_first,
        coordinate_convention=CONVENTION_THREE_PRIME_BASE_DROPPED,
    )
    assert dropped.mechanism == MECHANISM_MID_EXON
    assert dropped.three_prime_cds_start == 11
    assert dropped.frame_status == FRAME_OUT
    assert dropped.three_prime_location.genomic_position == three_first + strand


@pytest.mark.parametrize("strand", [1, -1])
def test_dropping_the_base_before_an_acceptor_lands_on_the_exon_first_base(strand):
    five, five_last, three, three_first = _exact_exon_edge_join(strand)
    before_exon = three_first - strand  # intron base -1
    dropped = _predict(
        five,
        five_last,
        three,
        before_exon,
        coordinate_convention=CONVENTION_THREE_PRIME_BASE_DROPPED,
    )
    assert dropped.mechanism == MECHANISM_MID_EXON
    assert dropped.three_prime_cds_start == 10
    assert dropped.frame_status == FRAME_IN


def _intronic_position(model, exon_index, nt, side, strand):
    """``nt`` intronic bases after (``"donor"``) or before (``"acceptor"``) an exon."""
    exon = model.exons[exon_index]
    if side == "donor":
        return exon.end + nt if strand == 1 else exon.start - nt
    return exon.start - nt if strand == 1 else exon.end + nt


@pytest.mark.parametrize("strand", [1, -1])
def test_five_prime_break_inside_the_donor_site_is_unresolved(strand):
    five = _two_exon("AAA", _PROTEIN, 75, strand=strand)
    three = _two_exon("BBB", _PROTEIN, 75, strand=strand)
    three_break = _intronic_position(three, 1, 300, "acceptor", strand)
    for nt in range(1, DONOR_SITE_NT):
        prediction = _predict(
            five, _intronic_position(five, 0, nt, "donor", strand), three, three_break
        )
        assert prediction.mechanism == MECHANISM_SPLICE_SITE_UNRESOLVED
        assert not prediction.has_product
        assert prediction.frame_status == FRAME_UNKNOWN
        assert f"keeps only {nt} nt of the exon 1 splice donor" in prediction.no_product_reason
        call = classify_family(prediction, five, three, FAMILY)
        assert call.classification == SPLICE_SITE_UNRESOLVED
    intact = _predict(
        five, _intronic_position(five, 0, DONOR_SITE_NT, "donor", strand), three, three_break
    )
    assert intact.mechanism == MECHANISM_SPLICE and intact.frame_status == FRAME_IN


@pytest.mark.parametrize("strand", [1, -1])
def test_three_prime_break_inside_the_acceptor_site_is_unresolved(strand):
    five = _two_exon("AAA", _PROTEIN, 75, strand=strand)
    three = _two_exon("BBB", _PROTEIN, 75, strand=strand)
    five_break = _intronic_position(five, 0, 300, "donor", strand)
    for nt in range(1, ACCEPTOR_SITE_NT):
        prediction = _predict(
            five, five_break, three, _intronic_position(three, 1, nt, "acceptor", strand)
        )
        assert prediction.mechanism == MECHANISM_SPLICE_SITE_UNRESOLVED
        assert "exon 2 splice acceptor" in prediction.no_product_reason
        region = "essential AG" if nt < 3 else "polypyrimidine tract"
        assert region in prediction.no_product_reason
        call = classify_family(prediction, five, three, FAMILY)
        assert call.classification == SPLICE_SITE_UNRESOLVED
    intact = _predict(
        five, five_break, three, _intronic_position(three, 1, ACCEPTOR_SITE_NT, "acceptor", strand)
    )
    assert intact.mechanism == MECHANISM_SPLICE and intact.frame_status == FRAME_IN


def test_splice_site_check_covers_an_intronic_break_paired_with_a_mid_exon_one():
    five = _two_exon("AAA", _PROTEIN, 75)
    three = _two_exon("BBB", _PROTEIN, 60)
    prediction = _predict(five, _intronic_position(five, 0, 2, "donor", 1), three, 1000 + 30)
    assert prediction.mechanism == MECHANISM_SPLICE_SITE_UNRESOLVED


def test_intron_on_five_prime_with_mid_exon_three_prime_uses_next_acceptor():
    five = _two_exon("AAA", _PROTEIN, 75)
    protein = _PROTEIN
    three = _two_exon("BBB", protein, 60)
    # 3' break inside exon 1, whose acceptor is lost: exon 2 (c.61) supplies one.
    prediction = _predict(five, 1000 + 75 + 300, three, 1000 + 30)
    assert prediction.mechanism == MECHANISM_SPLICE
    assert prediction.three_prime_cds_start == 61
    assert any("drops the exon 1 piece" in note for note in prediction.notes)


def test_mid_exon_five_prime_with_intronic_three_prime_uses_previous_donor():
    five = _two_exon("AAA", _PROTEIN, 60)
    three = _two_exon("BBB", _PROTEIN, 75)
    prediction = _predict(five, five.exons[1].start + 5, three, 1000 + 75 + 300)
    assert prediction.mechanism == MECHANISM_SPLICE
    assert prediction.five_prime_cds_retained == 60
    assert any("drops the exon 2 piece" in note for note in prediction.notes)


def test_mid_exon_five_prime_in_first_exon_with_intronic_three_prime_has_no_product():
    five = _two_exon("AAA", _PROTEIN, 60)
    three = _two_exon("BBB", _PROTEIN, 75)
    prediction = _predict(five, 1000 + 10, three, 1000 + 75 + 300)
    assert not prediction.has_product
    assert prediction.mechanism == MECHANISM_NONE
    assert "no upstream splice donor" in prediction.no_product_reason


def test_minus_strand_positions_run_high_to_low():
    three = _two_exon("BBB", _PROTEIN, 75, strand=-1)
    location = locate_breakpoint(three, three.exons[0].end)  # genomic high end of exon 1
    assert location.exon_rank == 1 and location.cds_position == 1
    assert location.at_exon_edge == "first_base"
    last = locate_breakpoint(three, three.exons[0].start)
    assert last.cds_position == 75 and last.at_exon_edge == "last_base"
    second = locate_breakpoint(three, three.exons[1].end)
    assert second.exon_rank == 2 and second.cds_position == 76


def test_breaks_outside_the_coding_region_have_no_product():
    five = _model("AAA", _PROTEIN, 11, 40)
    three = _model("BBB", _PROTEIN, 11, 40)
    before = _predict(five, 10, three, 1100)
    assert "outside the transcript" in before.no_product_reason

    utr5 = TranscriptModel(
        "CCC",
        "T_CCC",
        BUILD,
        1,
        (TranscriptExon(1, 900, 1000 + len(_cds(_PROTEIN)) - 1),),
        1000,
        1000 + len(_cds(_PROTEIN)) - 1,
        _cds(_PROTEIN),
    )
    in_utr = predict_junction("CCC", 950, utr5, "BBB", 1100, three, BUILD)
    assert "before its coding region" in in_utr.no_product_reason


def test_genome_build_is_explicit_and_checked():
    five = _model("AAA", _PROTEIN, 11, 40)
    three = _model("BBB", _PROTEIN, 11, 40)
    assert "no genome build" in _predict(five, 1050, three, 1100, build=None).no_product_reason
    mismatch = _predict(five, 1050, three, 1100, build="GRCh38")
    assert "GRCh37" in mismatch.no_product_reason and "GRCh38" in mismatch.no_product_reason
    assert _predict(five, 1050, three, 1100, build="hg19").has_product
    assert normalize_genome_build("hg38") == "GRCh38"
    with pytest.raises(ValueError):
        TranscriptModelSource(None)


def test_missing_or_inconsistent_transcript_model_reasons():
    three = _model("BBB", _PROTEIN, 11, 40)
    missing = predict_junction("AAA", 1050, None, "BBB", 1100, three, BUILD)
    assert "no transcript model for 5' partner AAA" == missing.no_product_reason

    bad = _model("AAA", _PROTEIN, 11, 40)
    bad = TranscriptModel(**{**bad.__dict__, "cds_sequence": bad.cds_sequence[:-3]})
    assert "exon structure" in _predict(bad, 1050, three, 1100).no_product_reason
    assert classify_family(missing, None, three, FAMILY).classification == NO_TRANSCRIPT_MODEL


def test_unknown_convention_is_rejected():
    five = _model("AAA", _PROTEIN, 11, 40)
    with pytest.raises(ValueError):
        _predict(five, 1050, five, 1100, coordinate_convention="nope")


def test_translate_and_alignment():
    assert translate("ATGGCCTAA") == "MA*"
    aligned_a, aligned_b, score = align_global("ACDEFGHIK", "ACDEGHIK")
    assert aligned_a == "ACDEFGHIK"
    assert aligned_b == "ACDE-GHIK"
    assert score > 0


def test_classes_are_documented_constants():
    assert CLEAN_RECONSTITUTION in RECONSTITUTION_CLASSES
    assert NO_TRANSCRIPT_MODEL in RECONSTITUTION_CLASSES
    assert SPLICE_SITE_UNRESOLVED in RECONSTITUTION_CLASSES


# --- Transcript source --------------------------------------------------------


def test_canonical_fallback_and_unmatched_annotation_are_labelled():
    source = _source()
    assert source.resolve("ERG").source == CANONICAL_NO_ANNOTATION
    unmatched = source.resolve("ERG", ["NM_999999"])
    assert unmatched.source == CANONICAL_ANNOTATION_UNMATCHED
    assert "NM_999999 not found" in unmatched.note
    assert unmatched.differs_from_canonical is False


def test_gene_absent_from_genome_nexus_is_reported():
    from cfh.mapping.genome_nexus_source import GenomeNexusGeneNotFound

    class Absent(FakeGenomeNexus):
        def fetch_canonical_transcript(self, gene):
            raise GenomeNexusGeneNotFound(gene)

    source = TranscriptModelSource(BUILD, genome_nexus_client=Absent(), cds_client=FakeCds())
    resolved = source.resolve("NOPE")
    assert resolved.transcript is None
    assert source.model(resolved) == (None, "Genome Nexus has no transcript for NOPE")


def test_cds_mismatched_to_protein_length_is_not_used():
    class Short(FakeCds):
        def fetch_cds(self, transcript_id):
            return super().fetch_cds(transcript_id)[3:]

    source = TranscriptModelSource(BUILD, genome_nexus_client=FakeGenomeNexus(), cds_client=Short())
    model, problem = source.model(source.resolve("ERG", ["NM_182918"]))
    assert model is None and "Genome Nexus reports 479" in problem


def test_models_and_isoforms_are_fetched_once():
    client = FakeGenomeNexus()
    cds = FakeCds()
    source = TranscriptModelSource(BUILD, genome_nexus_client=client, cds_client=cds)
    for _ in range(2):
        source.model(source.resolve("ERG", ["NM_182918"]))
    assert client.isoform_calls == ["ERG"]
    assert cds.calls == ["ENST00000288319"]


def test_on_disk_caches_are_reused(tmp_path):
    client = FakeGenomeNexus()
    source = TranscriptModelSource(
        BUILD, genome_nexus_client=client, cds_client=FakeCds(), cache_dir=tmp_path
    )
    source.resolve("ERG", ["NM_182918"])
    assert (tmp_path / "genome_nexus_canonical_transcripts_mskcc" / BUILD / "ERG.json").exists()
    assert (tmp_path / "genome_nexus_transcripts" / BUILD / "ERG.json").exists()

    class Offline(FakeGenomeNexus):
        def fetch_canonical_transcript(self, gene):
            raise AssertionError("should come from disk")

        def fetch_transcripts_for_gene(self, gene):
            raise AssertionError("should come from disk")

    reused = TranscriptModelSource(
        BUILD, genome_nexus_client=Offline(), cds_client=FakeCds(), cache_dir=tmp_path
    )
    assert reused.resolve("ERG", ["NM_182918"]).refseq_id == "NM_182918"


def test_ensembl_cds_client_builds_per_build_and_caches(tmp_path):
    from unittest.mock import MagicMock

    session = MagicMock()
    session.get.return_value.status_code = 200
    session.get.return_value.json.return_value = {"seq": "ATGTAA"}
    client = EnsemblCdsClient("GRCh38", session=session, cache_dir=tmp_path)
    assert client.fetch_cds("ENST1") == "ATGTAA"
    assert client.fetch_cds("ENST1") == "ATGTAA"
    assert session.get.call_count == 1
    assert session.get.call_args.args[0].startswith("https://rest.ensembl.org/")

    session.get.return_value.status_code = 400
    assert EnsemblCdsClient("GRCh37", session=session).fetch_cds("ENST2") is None
    assert EnsemblCdsClient("T2T", session=session).fetch_cds("ENST3") is None


def _shifted_payloads(gene, shift):
    """The fixture's payloads for ``gene`` with every coordinate moved by ``shift``."""
    payloads = json.loads(json.dumps(TRANSCRIPTS["transcripts"][gene]))
    for payload in payloads:
        for exon in payload.get("exons") or []:
            exon["exonStart"] += shift
            exon["exonEnd"] += shift
        for utr in payload.get("utrs") or []:
            utr["start"] += shift
            utr["end"] += shift
    return payloads


class FakeGrch38GenomeNexus(FakeGenomeNexus):
    """A GRCh38 instance: same transcripts, coordinates moved by a fixed offset."""

    base_url = "https://grch38.genomenexus.org"
    SHIFT = 1_000_000

    def _payloads(self, gene):
        return _shifted_payloads(gene, self.SHIFT)


def test_grch37_then_grch38_in_one_cache_dir_do_not_cross_contaminate(tmp_path):
    grch37 = TranscriptModelSource(
        "GRCh37", genome_nexus_client=FakeGenomeNexus(), cds_client=FakeCds(), cache_dir=tmp_path
    )
    model37, _ = grch37.model(grch37.resolve("ERG", ["NM_182918"]))
    canonical37 = grch37.canonical("FLI1")

    grch38_client = FakeGrch38GenomeNexus()
    grch38 = TranscriptModelSource(
        "GRCh38", genome_nexus_client=grch38_client, cds_client=FakeCds(), cache_dir=tmp_path
    )
    model38, _ = grch38.model(grch38.resolve("ERG", ["NM_182918"]))
    canonical38 = grch38.canonical("FLI1")

    shift = FakeGrch38GenomeNexus.SHIFT
    assert grch38_client.isoform_calls == ["ERG"]  # fetched, not read from the GRCh37 cache
    assert model38.genome_build == "GRCh38" and model37.genome_build == "GRCh37"
    assert [e.start for e in model38.exons] == [e.start + shift for e in model37.exons]
    assert canonical38.exons[0].start == canonical37.exons[0].start + shift
    for kind in ("genome_nexus_canonical_transcripts_mskcc", "genome_nexus_transcripts"):
        for build in ("GRCh37", "GRCh38"):
            marker = json.loads((tmp_path / kind / build / "_provenance.json").read_text())
            assert marker["genome_build"] == build

    # A second GRCh38 run reuses its own cache, still with GRCh38 coordinates.
    class Offline(FakeGrch38GenomeNexus):
        def fetch_canonical_transcript(self, gene):
            raise AssertionError("should come from disk")

        def fetch_transcripts_for_gene(self, gene):
            raise AssertionError("should come from disk")

    reused = TranscriptModelSource(
        "GRCh38", genome_nexus_client=Offline(), cds_client=FakeCds(), cache_dir=tmp_path
    )
    again, _ = reused.model(reused.resolve("ERG", ["NM_182918"]))
    assert [e.start for e in again.exons] == [e.start for e in model38.exons]


def test_cache_written_for_another_build_is_rejected_on_load(tmp_path):
    TranscriptModelSource(
        "GRCh37", genome_nexus_client=FakeGenomeNexus(), cds_client=FakeCds(), cache_dir=tmp_path
    ).resolve("ERG", ["NM_182918"])
    # Simulate a cache mislabelled as GRCh38 that really holds GRCh37 data.
    (tmp_path / "genome_nexus_transcripts" / "GRCh37").rename(
        tmp_path / "genome_nexus_transcripts" / "GRCh38"
    )
    source = TranscriptModelSource(
        "GRCh38",
        genome_nexus_client=FakeGrch38GenomeNexus(),
        cds_client=FakeCds(),
        cache_dir=tmp_path,
    )
    with pytest.raises(ValueError, match="GRCh37"):
        source.resolve("ERG", ["NM_182918"])


def test_genome_nexus_instance_follows_the_requested_build():
    assert TranscriptModelSource("GRCh38").client.base_url == "https://grch38.genomenexus.org"
    assert TranscriptModelSource("hg19").client.base_url == "https://www.genomenexus.org"
    with pytest.raises(ValueError, match="serves GRCh37"):
        TranscriptModelSource("GRCh38", genome_nexus_client=FakeGenomeNexus())
    with pytest.raises(ValueError, match="no public Genome Nexus"):
        TranscriptModelSource("T2T")


def test_ensembl_cds_cache_rejects_another_build(tmp_path):
    from unittest.mock import MagicMock

    session = MagicMock()
    session.get.return_value.status_code = 200
    session.get.return_value.json.return_value = {"seq": "ATGTAA"}
    EnsemblCdsClient("GRCh37", session=session, cache_dir=tmp_path).fetch_cds("ENST1")
    with pytest.raises(ValueError, match="GRCh37"):
        EnsemblCdsClient("GRCh38", session=session, cache_dir=tmp_path).fetch_cds("ENST1")


def test_genome_build_of_instance():
    assert genome_build_of_genome_nexus("https://grch38.genomenexus.org") == "GRCh38"
    assert genome_build_of_genome_nexus("https://www.genomenexus.org") == "GRCh37"


# --- Cohort screen -----------------------------------------------------------

KEY_FAMILY = "PF07714"
OTHER_FAMILY = "PF00018"


class SyntheticGenomeNexus:
    """Two genes sharing a key family and a non-key family on one 80-aa protein."""

    base_url = "https://example.invalid"
    session = None

    def __init__(self):
        self.protein = "M" + "S" * 9 + DOMAIN + "G" * 40
        self.cds = _cds(self.protein)

    def _payload(self, gene, offset):
        start = 1000 + offset
        return {
            "transcriptId": f"ENST_{gene}",
            "refseqMrnaId": {"AAA": "NM_100001", "BBB": "NM_100002"}[gene],
            "proteinLength": len(self.protein),
            "proteinId": f"ENSP_{gene}",
            "pfamDomains": [
                {"pfamDomainId": KEY_FAMILY, "pfamDomainStart": 55, "pfamDomainEnd": 75},
                {"pfamDomainId": OTHER_FAMILY, "pfamDomainStart": 11, "pfamDomainEnd": 40},
            ],
            "exons": [
                {
                    "exonId": f"E_{gene}",
                    "exonStart": start,
                    "exonEnd": start + len(self.cds) - 1,
                    "rank": 1,
                    "strand": 1,
                }
            ],
            "utrs": [],
        }

    def fetch_canonical_transcript(self, gene):
        return self._payload(gene, 0 if gene == "AAA" else 5_000_000)

    def fetch_transcripts_for_gene(self, gene):
        return [self.fetch_canonical_transcript(gene)]


class SyntheticCds:
    def __init__(self, cds):
        self.cds = cds

    def fetch_cds(self, transcript_id):
        return self.cds


def _screen_frame():
    five_break = 1000 + 25 * 3 - 1
    three_break = 5_001_000 + 25 * 3
    row = {
        "Sample_Id": "S1",
        "Site1_Hugo_Symbol": "AAA",
        "Site1_Chromosome": "1",
        "Site1_Position": five_break,
        "Site2_Hugo_Symbol": "BBB",
        "Site2_Chromosome": "2",
        "Site2_Position": three_break,
        "Site2_Effect_On_Frame": "in-frame",
        "Tumor_Split_Read_Count": 7,
        "Tumor_Paired_End_Read_Count": 3,
        "Tumor_Variant_Count": 9,
        "Normal_Variant_Count": 0,
        "SV_Status": "SOMATIC",
        "NCBI_Build": BUILD,
        "Connection_Type": "3to5",
        "Breakpoint_Type": "IMPPRECISE",
        "Annotation": "AAA (NM_100001) - BBB (NM_100002) rearrangement",
        "Event_Info": "Protein Fusion: mid-exon {AAA:BBB}",
        "Extra_fields": None,
        "Source_row_number": 1,
        "Parse_warnings": None,
    }
    return pd.DataFrame([row])


def _synthetic_source():
    client = SyntheticGenomeNexus()
    source = TranscriptModelSource(
        BUILD, genome_nexus_client=client, cds_client=SyntheticCds(client.cds)
    )
    source.prefetch_canonical = lambda genes: None
    return source


def test_any_domain_mode_finds_a_non_key_domain_key_domain_mode_misses():
    source = _synthetic_source()
    any_rows, any_summary = screen_cohort(_screen_frame(), source, cohort="demo")
    by_family = {row["pfam_family"]: row for row in any_rows}
    assert set(by_family) == {KEY_FAMILY, OTHER_FAMILY}
    assert by_family[OTHER_FAMILY]["reconstitution_class"] == CLEAN_RECONSTITUTION
    assert by_family[KEY_FAMILY]["reconstitution_class"] == NOT_SPLIT_ON_BOTH
    assert any_summary["class_counts"] == {CLEAN_RECONSTITUTION: 1, NOT_SPLIT_ON_BOTH: 1}

    key_rows, key_summary = screen_cohort(
        _screen_frame(),
        _synthetic_source(),
        cohort="demo",
        mode=MODE_KEY_DOMAIN,
        key_family=lambda gene: KEY_FAMILY,
    )
    assert [row["pfam_family"] for row in key_rows] == [KEY_FAMILY]
    assert key_rows[0]["reconstitution_class"] == NOT_SPLIT_ON_BOTH
    assert key_summary["family_rows"] == 1


def test_key_domain_mode_requires_the_family_mapping():
    with pytest.raises(ValueError):
        screen_cohort(_screen_frame(), _synthetic_source(), cohort="demo", mode=MODE_KEY_DOMAIN)
    with pytest.raises(ValueError):
        screen_cohort(_screen_frame(), _synthetic_source(), cohort="demo", mode="nope")


def test_screen_reports_read_support_transcripts_and_columns(tmp_path):
    rows, summary = screen_cohort(_screen_frame(), _synthetic_source(), cohort="demo")
    row = next(item for item in rows if item["pfam_family"] == OTHER_FAMILY)
    assert row["tumor_split_reads"] == 7
    assert row["tumor_paired_end_reads"] == 3
    assert row["tumor_variant_count"] == 9
    assert row["normal_variant_count"] == 0
    assert row["breakpoint_precision"] == "IMPRECISE"
    assert row["five_prime_transcript"] == "NM_100001"
    assert row["five_prime_transcript_source"] == SOURCE_ANNOTATED
    assert row["five_prime_transcript_differs_from_canonical"] is False
    assert row["junction_mechanism"] == MECHANISM_MID_EXON

    paths = write_screen_outputs(rows, summary, tmp_path)
    with paths["domain_screen_tsv"].open() as handle:
        header = handle.readline().rstrip("\n").split("\t")
    assert header == list(SCREEN_COLUMNS)
    assert json.loads(paths["domain_screen_summary"].read_text())["family_rows"] == 2


def test_screen_skips_non_protein_fusions():
    frame = _screen_frame()
    frame.loc[0, "Event_Info"] = "Antisense Fusion"
    rows, summary = screen_cohort(frame, _synthetic_source(), cohort="demo")
    assert rows == [] and summary["screened_fusion_events"] == 0
    event = SimpleNamespace(
        Five_prime_gene="A",
        Three_prime_gene="A",
        Is_protein_fusion=True,
        Is_antisense=False,
        Event_info="Protein Fusion",
    )
    assert not is_screenable_fusion(event)


def test_screen_records_a_breakpoint_outside_the_locus():
    frame = _screen_frame()
    frame.loc[0, "Site2_Position"] = 99_000_000
    frame.loc[0, "Site2_Chromosome"] = "7"
    rows, _ = screen_cohort(frame, _synthetic_source(), cohort="demo")
    row = next(item for item in rows if item["pfam_family"] == OTHER_FAMILY)
    assert row["reconstitution_class"] == NO_TRANSCRIPT_MODEL
    assert "fusion has no site breakpoint" in row["class_reason"]


# --- Read support and archive parsing -----------------------------------------


def test_real_msk_header_read_counts_parse():
    df = parse_sv_file(MSK_ARCHIVE_FIXTURE)
    fli1 = df[df["Sample_Id"] == "P-0053737-T01-IM6"].iloc[0]
    assert fli1["Tumor_Paired_End_Read_Count"] == 6
    assert fli1["Tumor_Split_Read_Count"] == 0
    assert fli1["Tumor_Variant_Count"] == 8
    assert fli1["Normal_Variant_Count"] == 0
    gabpa = df[df["Sample_Id"] == "P-0023818-T01-IM6"].iloc[0]
    assert gabpa["Tumor_Split_Read_Count"] == 10
    assert gabpa["Tumor_Paired_End_Read_Count"] == 6
    for warnings in df["Parse_warnings"].dropna():
        assert "Tumor_Split_Read_Count" not in warnings
        assert "Tumor_Paired_End_Read_Count" not in warnings
    # Aliased and read-support columns are not reported as unknown extras.
    assert "Split_Read_Support" not in (fli1["Extra_fields"] or {})
    assert "Normal_Variant_Count" not in (fli1["Extra_fields"] or {})


def test_na_in_a_numeric_read_support_column_is_missing_not_an_error():
    df = parse_sv_file(MSK_ARCHIVE_FIXTURE)
    row = df[df["Normal_Variant_Count"].isna()].iloc[0]
    assert pd.isna(row["Normal_Variant_Count"])
    assert "could not parse Normal_Variant_Count" not in str(row["Parse_warnings"])


def test_read_support_sentinels_become_blank():
    support = read_support_columns(
        {
            "Tumor_Split_Read_Count": -1,
            "Tumor_Paired_End_Read_Count": "bad",
            "Breakpoint_Type": "precise",
        }
    )
    assert support == {
        "tumor_split_reads": None,
        "tumor_paired_end_reads": None,
        "tumor_variant_count": None,
        "normal_variant_count": None,
        "breakpoint_precision": "PRECISE",
    }


def test_event_columns_keep_one_column_set():
    filled = junction_call_columns({"reconstitution_class": "x", "pfam_family": "PF1"}, {})
    blank = junction_call_columns(None)
    assert filled.keys() == blank.keys()
    assert filled["reconstitution_class"] == "x"
    assert filled["reconstitution_pfam_family"] == "PF1"
    assert "reconstitution_tumor_split_reads" in blank
    assert set(blank.values()) == {None}


def test_per_gene_pipeline_call_uses_the_key_family_and_partner_locus():
    from cfh.real_benchmark import _junction_call

    row = {
        "Site1_Hugo_Symbol": "FLI1",
        "Site1_Chromosome": "21",
        "Site1_Position": 39755655,  # the archive swaps the two sites' genes
        "Site2_Hugo_Symbol": "ERG",
        "Site2_Chromosome": "11",
        "Site2_Position": 128680564,
        "Annotation": FLI1_ERG_ANNOTATION,
    }
    call = _junction_call(
        row,
        SimpleNamespace(Reference_build=BUILD),
        "ERG",
        "three_prime",
        39755655,
        "FLI1",
        ETS,
        _source(),
    )
    assert call["pfam_family"] == ETS
    assert call["reconstitution_class"] == TANDEM_DUPLICATION_OR_INSERTION
    assert "LSRALRYYYDKNI|YYYDKNIMTKV" in call["junction_peptide"]
    columns = junction_call_columns(call, {"Tumor_Paired_End_Read_Count": 6})
    assert columns["reconstitution_class"] == TANDEM_DUPLICATION_OR_INSERTION
    assert columns["reconstitution_tumor_paired_end_reads"] == 6
