import csv
import gzip
import json
import random

import pytest

from cfh.cohort.artefact_flags import (
    FLAG_COLUMNS,
    FLAGGED,
    NOT_APPLICABLE,
    NOT_ASSESSED,
    PASS,
    EnsemblSequenceSource,
    RepeatTracks,
    annotate_rows,
    best_single_sequence_match,
    breakpoint_similarity,
    chromosome_name,
    genomic_flags,
    junction_key,
    local_align,
    normal_flag,
    paralog_flag,
    reverse_complement,
    sample_call_counts,
    tcga_gtex_index,
)
from cfh.ingestion.open_fusion_tables import SCHEMA


def _random_dna(length, seed):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


def _mutate(sequence, every):
    swap = {"A": "C", "C": "G", "G": "T", "T": "A"}
    return "".join(swap[b] if i % every == every - 1 else b for i, b in enumerate(sequence))


# --- alignment ----------------------------------------------------------------------


def test_local_align_identical_mismatched_gapped_and_unrelated():
    a = _random_dna(200, 1)
    same = local_align(a, a)
    assert (same.identity, same.columns, same.query_span) == (1.0, 200, 200)

    mismatched = local_align(a, _mutate(a, 10))
    assert mismatched.columns >= 190 and mismatched.identity == pytest.approx(0.9, abs=0.01)

    gapped = local_align(a, a[:100] + "TTT" + a[100:])
    assert gapped.matches == 200 and gapped.columns == 203

    unrelated = local_align(_random_dna(200, 2), _random_dna(200, 3))
    assert unrelated.columns < 60


def test_local_align_embeds_query_in_longer_target_and_never_matches_n():
    query = _random_dna(100, 4)
    target = _random_dna(500, 5) + query + _random_dna(500, 6)
    hit = local_align(query, target)
    assert (hit.identity, hit.query_span) == (1.0, 100)
    assert local_align("NNNNNNNN", "NNNNNNNN").score == 0
    assert local_align("", "ACGT").columns == 0


def test_reverse_complement_and_chromosome_names():
    assert reverse_complement("AACGTn") == "nACGTT"
    assert chromosome_name("chr7") == "7" and chromosome_name("X") == "X"
    assert chromosome_name("chrM") == "MT" and chromosome_name("") is None


def test_single_sequence_explains_a_window_only_when_it_covers_it():
    five, three = _random_dna(400, 7), _random_dna(400, 8)
    chimera = five[:200] + three[200:]
    window = chimera[150:250]
    match = best_single_sequence_match(window, {"five": five, "three": three})
    assert match is not None and not match.explains
    assert match.coverage == pytest.approx(0.5, abs=0.05)

    paralog = five[:200] + _mutate(three[200:], 25)  # three's 3' half, 96% identical
    explained = best_single_sequence_match(
        window, {"five": five, "three": three, "paralog": paralog}
    )
    assert explained is not None and explained.explains and explained.target == "paralog"


def test_batched_target_scores_match_single_alignments_and_respect_spacers(monkeypatch):
    import cfh.cohort.artefact_flags as flags

    query = _random_dna(100, 20)
    targets = {f"t{i}": _random_dna(300 + 37 * i, 30 + i) for i in range(6)}
    targets["half"] = _random_dna(200, 40) + query[:50]
    targets["tail"] = query[50:] + _random_dna(200, 41)  # next to "half" when joined
    monkeypatch.setattr(flags, "CHUNK_BASES", 900)  # force several chunks
    scores = flags.target_scores(query, targets)
    assert scores == {name: local_align(query, seq).score for name, seq in targets.items()}
    assert scores["half"] == 50 and scores["tail"] == 50  # no alignment bridges the two
    assert best_single_sequence_match("AAAA", {"t": "CCCC"}) is None


def test_breakpoint_similarity_reports_orientation():
    a = _random_dna(301, 9)
    assert breakpoint_similarity(a, a, 150).high
    opposite = breakpoint_similarity(a, reverse_complement(a), 150)
    assert opposite.high and opposite.orientation == "opposite_strand"
    assert not breakpoint_similarity(a, _random_dna(301, 10), 150).high
    short = a[:90]
    assert not breakpoint_similarity(short, short, 150).high  # < 100 aligned bp


# --- repeat tracks -----------------------------------------------------------------


def _gz(path, rows):
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerows(rows)
    return path


def _segdup(chrom, start, end):
    row = ["0", chrom, str(start), str(end), "chr9:1", "0", "+", "chr9", "100", "200"]
    row += ["0"] * 16 + ["0.987"] + ["0"] * 3
    return row


def _rmsk(chrom, start, end, name="AluY", klass="SINE", family="Alu"):
    return ["0", "0", "0", "0", "0", chrom, str(start), str(end), "0", "+", name, klass, family]


def test_repeat_tracks_use_zero_based_half_open_intervals(tmp_path):
    tracks = RepeatTracks(
        "GRCh38",
        _gz(tmp_path / "s.gz", [_segdup("chr1", 100, 200)]),
        _gz(tmp_path / "r.gz", [_rmsk("chr2", 10, 20)]),
    )
    hits = tracks.overlaps([("chr1", 100), ("1", 101), ("chr1", 200), ("chr1", 201), ("2", 15)])
    assert ("1", 100) not in hits
    assert [h.track for h in hits[("1", 101)]] == ["segdup"]
    assert "fracMatch=0.987" in hits[("1", 200)][0].detail
    assert ("1", 201) not in hits
    assert hits[("2", 15)][0].name == "AluY" and hits[("2", 15)][0].detail == "SINE/Alu"


def test_repeat_tracks_report_each_record_once_for_a_repeated_position(tmp_path):
    tracks = RepeatTracks(
        "GRCh38",
        _gz(tmp_path / "s.gz", [_segdup("chr1", 100, 200)]),
        _gz(tmp_path / "r.gz", [_rmsk("chr1", 150, 160)]),
    )
    hits = tracks.overlaps([("chr1", 155), ("1", 155), ("chr1", 155)])
    assert [h.track for h in hits[("1", 155)]] == ["segdup", "rmsk"]


# --- fakes -------------------------------------------------------------------------


class FakeSequences:
    def __init__(self, genome_build="GRCh38", regions=None, cds=None, paralogs=None):
        self.genome_build = genome_build
        self.regions = regions or {}
        self.cds = cds or {}
        self.paralogs = paralogs or {}
        self.region_calls = []

    def region(self, chromosome, start, end):
        self.region_calls.append((chromosome, start, end))
        base = self.regions.get(chromosome)
        return base[start - 1 : end] if base else None

    def paralog_cds(self, gene):
        return self.paralogs.get(gene, {})

    def transcript_cds(self, transcript_id):
        return self.cds.get(transcript_id)


def _row(**overrides):
    row = {
        "cohort": "TARGET",
        "sample_id": "S1",
        "genome_build": "GRCh38",
        "reconstitution_class": "clean_reconstitution",
        "five_prime_gene": "AAA",
        "three_prime_gene": "BBB",
        "five_prime_chromosome": "chr1",
        "three_prime_chromosome": "chr2",
        "five_prime_breakpoint": 1000,
        "three_prime_breakpoint": 2000,
        "five_prime_ensembl_transcript": "ENST_A",
        "three_prime_ensembl_transcript": "ENST_B",
        "junction_context_nt": None,
        "tumor_split_reads": 5,
        "tumor_paired_end_reads": 2,
        "breakpoint_precision": None,
        "Open_source_table": "STAR-Fusion.v1.10.0.TARGETdb.tsv.gz",
        "Open_tissue_class": "tumor",
    }
    row.update(overrides)
    return row


# --- paralog flag -------------------------------------------------------------------


def test_paralog_flag_is_not_assessed_without_a_fused_cds():
    out = paralog_flag(_row(reconstitution_class="splice_site_unresolved"), FakeSequences())
    assert out["paralog_sequence_flag"] == NOT_ASSESSED
    assert "splice_site_unresolved" in out["paralog_sequence_reason"]


def test_paralog_flag_passes_a_true_chimera_and_flags_a_paralog_explained_one():
    five, three = _random_dna(400, 11), _random_dna(400, 12)
    window = (five[:200] + three[200:])[150:250]
    sequences = FakeSequences(cds={"ENST_A": five, "ENST_B": three})
    out = paralog_flag(_row(junction_context_nt=window), sequences)
    assert out["paralog_sequence_flag"] == PASS
    assert "best single sequence (2 tested)" in out["paralog_sequence_reason"]

    sequences.paralogs = {"BBB": {"BBB:paralog:ENSG_P": five[:200] + three[200:]}}
    flagged = paralog_flag(_row(junction_context_nt=window), sequences)
    assert flagged["paralog_sequence_flag"] == FLAGGED
    assert flagged["paralog_best_target"] == "BBB:paralog:ENSG_P"
    assert flagged["paralog_best_identity"] == 1.0


def test_paralog_flag_is_not_assessed_when_no_sequence_is_available():
    out = paralog_flag(_row(junction_context_nt="ACGT" * 25), FakeSequences())
    assert out["paralog_sequence_flag"] == NOT_ASSESSED


# --- genomic flags -----------------------------------------------------------------


def _genome(similar):
    one = _random_dna(4000, 13)
    two = _random_dna(4000, 14)
    if similar:
        two = two[:1849] + one[849:1150] + two[2150:]  # bases 850-1150 copied to 1850-2150
    return {"1": one, "2": two}


def test_genomic_flags_report_similarity_and_both_windows():
    out = genomic_flags(_row(), FakeSequences(regions=_genome(True)), {}, 1)
    assert out["genomic_similarity_flag"] == FLAGGED
    assert out["genomic_identity_150"] == 1.0 and out["genomic_aligned_length_150"] >= 290
    assert out["genomic_aligned_length_500"] >= 290
    assert out["breakpoint_repeat_flag"] == PASS
    assert out["sample_burden_flag"] == PASS and out["read_evidence_flag"] == PASS

    clean = genomic_flags(_row(), FakeSequences(regions=_genome(False)), {}, 1)
    assert clean["genomic_similarity_flag"] == PASS
    assert clean["genomic_aligned_length_500"] < 100


def test_genomic_flags_never_pass_silently_without_repeat_tracks_or_sequence():
    out = genomic_flags(
        _row(), FakeSequences(regions=_genome(False)), None, 1, "no tracks for this build"
    )
    assert out["breakpoint_repeat_flag"] == NOT_ASSESSED
    assert out["breakpoint_repeat_reason"] == "no tracks for this build"
    assert out["genomic_similarity_flag"] == NOT_ASSESSED

    missing = genomic_flags(_row(three_prime_chromosome=None), FakeSequences(), {}, 1)
    assert missing["genomic_similarity_flag"] == NOT_ASSESSED
    assert "missing" in missing["genomic_similarity_reason"]


def test_genomic_flags_sample_burden_repeats_and_read_evidence():
    from cfh.cohort.artefact_flags import TrackHit

    hits = {("1", 1000): [TrackHit("segdup", "chr9:1", "other=chr9:101-200 fracMatch=0.990")]}
    out = genomic_flags(
        _row(tumor_split_reads=0, breakpoint_precision="IMPRECISE"),
        FakeSequences(regions=_genome(False)),
        hits,
        3,
    )
    assert out["genomic_similarity_flag"] == FLAGGED
    assert (
        out["breakpoint_repeat_flag"] == FLAGGED
        and "5' in segdup" in out["breakpoint_repeat_reason"]
    )
    assert out["sample_burden_flag"] == FLAGGED and out["sample_shared_family_calls"] == 3
    assert out["read_evidence_flag"] == FLAGGED
    assert "0 split/junction reads (2 pairs)" in out["read_evidence_reason"]
    assert "imprecise" in out["genomic_similarity_reason"]


def test_sample_call_counts_count_distinct_shared_family_events():
    rows = [
        _row(),
        _row(pfam_family="PF2"),  # same event, second family
        _row(three_prime_gene="CCC"),
        _row(reconstitution_class="no_shared_family", three_prime_gene="DDD"),
        _row(sample_id="S2"),
    ]
    assert sample_call_counts(rows) == {("TARGET", "S1"): 2, ("TARGET", "S2"): 1}


# --- normal recurrence ---------------------------------------------------------------


def _tcga_table(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(SCHEMA)
        for values in rows:
            row = {column: "" for column in SCHEMA}
            row.update(values)
            row.setdefault("annotation_flags", "{}")
            writer.writerow([row[column] for column in SCHEMA])
    return path


def _tcga_row(sample, tissue_class, gene_5p="AAA", gene_3p="BBB", pos_5p="1000", label=""):
    return {
        "sample_id": sample,
        "tissue_class": tissue_class,
        "gene_5p": gene_5p,
        "gene_3p": gene_3p,
        "chrom_5p": "chr1",
        "pos_5p": pos_5p,
        "chrom_3p": "chr2",
        "pos_3p": "2000",
        "annotation_flags": json.dumps({"tissue_subtype": label}),
    }


def test_normal_recurrence_counts_exact_normal_junctions_only(tmp_path):
    table = _tcga_table(
        tmp_path / "tcga.tsv",
        [
            _tcga_row("GTEX-1-A", "normal", label="Whole Blood"),
            _tcga_row("GTEX-2-A", "cell_line", label="Cells - EBV-transformed lymphocytes"),
            _tcga_row("GTEX-3-A", "normal", gene_3p="BBB-AS1"),
            _tcga_row("GTEX-4-A", "normal", pos_5p="1001"),
            _tcga_row("LUSC-1-TP", "tumor"),
            _tcga_row("S1", "normal"),  # the candidate's own sample
        ],
    )
    key = ("1", 1000, "2", 2000)
    exact, pairs = tcga_gtex_index(table, {key}, {("AAA", "BBB")})
    out = normal_flag(_row(), key, None, exact, pairs)
    assert out["normal_recurrence_flag"] == FLAGGED
    assert out["normal_exact_junction_samples"] == 1
    assert out["cell_line_exact_junction_samples"] == 1
    assert "Whole Blood x1" in out["normal_recurrence_reason"]
    assert "other gene symbols AAA::BBB-AS1 (not counted)" in out["normal_recurrence_reason"]
    assert out["normal_gene_pair_samples"] == 3

    only_cells = tcga_gtex_index(
        _tcga_table(tmp_path / "cells.tsv", [_tcga_row("GTEX-2-A", "cell_line")]), {key}, set()
    )
    assert normal_flag(_row(), key, None, *only_cells)["normal_recurrence_flag"] == PASS


def test_normal_recurrence_flags_a_call_from_a_normal_sample_and_reports_missing_keys():
    own = _row(Open_source_table="TCGA_n_GTEx.STAR-Fusion.v1.7.fusion_preds.tsv.gz")
    own["Open_tissue_class"] = "normal"
    assert normal_flag(own, ("1", 1, "2", 2), None, {}, {})["normal_recurrence_flag"] == FLAGGED
    missing = normal_flag(_row(), None, "no liftover", {}, {})
    assert missing["normal_recurrence_flag"] == NOT_ASSESSED
    assert "no liftover" in missing["normal_recurrence_reason"]


def test_junction_key_uses_grch38_directly_and_lifts_grch37():
    assert junction_key(_row(), None) == (("1", 1000, "2", 2000), None)
    old = _row(genome_build="GRCh37", five_prime_chromosome="1", three_prime_chromosome="2")
    key, note = junction_key(old, lambda chrom, pos: (chrom, pos + 10))
    assert key == ("1", 1010, "2", 2010) and "lifted" in note
    assert junction_key(old, lambda chrom, pos: None)[0] is None
    assert junction_key(old, None)[0] is None


# --- end to end ----------------------------------------------------------------------


def test_annotate_rows_fills_every_column_and_marks_unshared_rows(tmp_path):
    table = _tcga_table(tmp_path / "tcga.tsv", [_tcga_row("GTEX-1-A", "normal")])
    rows = [
        _row(),
        _row(reconstitution_class="no_shared_family", three_prime_gene="ZZZ"),
    ]
    annotated = annotate_rows(
        rows,
        sequences={"GRCh38": FakeSequences(regions=_genome(False))},
        repeat_tracks={"GRCh38": None},
        repeat_unavailable={"GRCh38": "tracks not downloaded"},
        tcga_gtex_table=table,
    )
    assert all(set(FLAG_COLUMNS) <= set(row) for row in annotated)
    assert annotated[0]["normal_recurrence_flag"] == FLAGGED
    assert annotated[0]["paralog_sequence_flag"] == NOT_ASSESSED
    assert annotated[0]["breakpoint_repeat_reason"] == "tracks not downloaded"
    assert annotated[1]["genomic_similarity_flag"] == NOT_APPLICABLE
    assert annotated[1]["normal_recurrence_flag"] == NOT_APPLICABLE


# --- Ensembl client ------------------------------------------------------------------


class FakeResponse:
    def __init__(self, payload, status=200):
        self.payload, self.status_code, self.headers = payload, status, {}

    def json(self):
        return self.payload

    def raise_for_status(self):
        if self.status_code >= 400:
            raise RuntimeError(self.status_code)


class FakeSession:
    def __init__(self, responses):
        self.responses = responses
        self.calls = []

    def get(self, url, params=None, headers=None, timeout=None):
        self.calls.append((url, params))
        for fragment, response in self.responses.items():
            if fragment in url:
                return response
        return FakeResponse({}, 404)


def test_ensembl_source_caches_per_build_and_parses_homology(tmp_path):
    session = FakeSession(
        {
            "/sequence/region/human/7:10..12:1": FakeResponse({"seq": "ACG"}),
            "/homology/symbol/human/ERG": FakeResponse(
                {
                    "data": [
                        {
                            "homologies": [
                                {
                                    "source": {"id": "ENSG_ERG", "seq": "ATGAAA"},
                                    "target": {"id": "ENSG_FLI1", "seq": "ATGCCC"},
                                }
                            ]
                        }
                    ]
                }
            ),
            "/map/human/GRCh37/21:5..5:1/GRCh38": FakeResponse(
                {"mappings": [{"mapped": {"seq_region_name": "21", "start": 9, "end": 9}}]}
            ),
        }
    )
    source = EnsemblSequenceSource("GRCh38", tmp_path, session=session)
    assert source.region("chr7", 10, 12) == "ACG"
    assert session.calls[0][1] == {"coord_system_version": "GRCh38"}
    assert source.paralog_cds("ERG") == {
        "ERG:compara_canonical:ENSG_ERG": "ATGAAA",
        "ERG:paralog:ENSG_FLI1": "ATGCCC",
    }
    assert source.transcript_cds("ENST_missing") is None

    reread = EnsemblSequenceSource("GRCh38", tmp_path, session=FakeSession({}))
    assert reread.region("7", 10, 12) == "ACG"
    assert (tmp_path / "ensembl_sequence" / "GRCh38" / "_provenance.json").exists()

    with pytest.raises(ValueError):
        source.liftover_to_grch38("21", 5)
    old = EnsemblSequenceSource("GRCh37", tmp_path, session=session)
    assert old.liftover_to_grch38("21", 5) == ("21", 9)
    assert old.host == "https://grch37.rest.ensembl.org"
