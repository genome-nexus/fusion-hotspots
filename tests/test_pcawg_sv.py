"""PCAWG data validation and physical breakend-to-transcript strand logic."""

import gzip
import io
import tarfile
from pathlib import Path

import pytest

from cfh.ingestion.pcawg_sv import (
    HEADER,
    Gene,
    GeneIndex,
    pair_genes,
    parse_bedpe,
    read_archive,
    read_gencode,
    sample_index,
)
from cfh.normalization.event_normalizer import normalize

BASE = [
    "1",
    "99",
    "100",
    "2",
    "199",
    "200",
    "SVMERGE1",
    "12",
    "+",
    "-",
    "TRA",
    "SNOWMAN_BRASS_dRANGER_DELLY",
]
SAMPLE = {
    "aliquot_id": "aliquot",
    "donor_unique_id": "BRCA-US::donor",
    "donor_wgs_exclusion_white_gray": "Whitelist",
    "dcc_project_code": "BRCA-US",
    "icgc_specimen_id": "SP1",
    "submitter_specimen_id": "specimen",
    "icgc_donor_id": "DO1",
    "submitter_donor_id": "donor",
    "library_strategy": "WGS",
    "dcc_specimen_type": "Primary tumour - solid tissue",
}


def bedpe(values=None):
    return io.StringIO("\t".join(HEADER) + "\n" + "\t".join(values or BASE) + "\n")


def event(**updates):
    return {**parse_bedpe(bedpe(), "test")[0], **SAMPLE, **updates}


def gene(symbol, chromosome, strand, start=0, end=1000):
    return Gene(f"ENSG_{symbol}", symbol, chromosome, start, end, strand, "protein_coding")


def test_parse_preserves_callers_and_source_support():
    row = parse_bedpe(bedpe(), "sample.bedpe")[0]
    assert row["caller_support"] == 4
    assert row["pe_support"] == 12
    assert row["start1"] == 99 and row["end1"] == 100
    assert row["source_line"] == 2
    assert "Tumor_Split_Read_Count" not in row


@pytest.mark.parametrize(
    ("column", "value"),
    [
        ("start1", "x"),
        ("start2", "-1"),
        ("end1", "99"),
        ("end2", "201"),
        ("chrom1", "chr1"),
        ("chrom2", "unknown"),
        ("strand1", "."),
        ("strand2", "?"),
        ("pe_support", "-1"),
        ("pe_support", "1.5"),
        ("svmethod", "BRASS_UNKNOWN"),
        ("svmethod", "BRASS_BRASS"),
        ("svclass", "INV"),
        ("sv_id", ""),
        ("svclass", "DEL"),  # different chromosomes
    ],
)
def test_malformed_rows_fail_with_source_line(column, value):
    values = list(BASE)
    values[HEADER.index(column)] = value
    with pytest.raises(ValueError, match="bad.bedpe:2:"):
        parse_bedpe(bedpe(values), "bad.bedpe")


@pytest.mark.parametrize(
    "text",
    [
        "wrong\theader\n",
        "\t".join(HEADER) + "\nshort\n",
        "\t".join(HEADER) + "\n" + "\t".join(BASE) + "\textra\n",
    ],
)
def test_invalid_structure(text):
    with pytest.raises(ValueError):
        parse_bedpe(io.StringIO(text), "bad")


def test_duplicate_sv_ids_fail():
    text = bedpe().getvalue() + "\t".join(BASE) + "\n"
    with pytest.raises(ValueError, match="duplicate SV id"):
        parse_bedpe(io.StringIO(text), "duplicate")


@pytest.mark.parametrize("strand1", ["+", "-"])
@pytest.mark.parametrize("strand2", ["+", "-"])
@pytest.mark.parametrize("gene1_strand", ["+", "-"])
@pytest.mark.parametrize("gene2_strand", ["+", "-"])
def test_all_sixteen_strand_combinations(strand1, strand2, gene1_strand, gene2_strand):
    row = event(strand1=strand1, strand2=strand2)
    pairs, _ = pair_genes(
        row, GeneIndex([gene("A", "1", gene1_strand), gene("B", "2", gene2_strand)])
    )
    first_five, second_five = strand1 == gene1_strand, strand2 == gene2_strand
    assert len(pairs) == int(first_five != second_five)
    if pairs:
        pair = pairs[0]
        assert pair["Site1_Hugo_Symbol"] == ("A" if first_five else "B")
        assert pair["Site2_Hugo_Symbol"] == ("B" if first_five else "A")
        assert pair["Site1_Position"] == (100 if first_five else 200)
        assert pair["Site2_Position"] == (200 if first_five else 100)
        assert pair["Tumor_Split_Read_Count"] is None
        assert pair["Tumor_Paired_End_Read_Count"] is None
        assert pair["NCBI_Build"] == "GRCh37"
        import pandas as pd

        normalized = normalize(pd.DataFrame(pairs), None, "pcawg_test")[0]
        assert normalized.Five_prime_gene == pair["Site1_Hugo_Symbol"]
        assert normalized.Three_prime_gene == pair["Site2_Hugo_Symbol"]
        assert normalized.Frame_status == "unknown"


def test_intergenic_and_same_gene_are_audited():
    pairs, audit = pair_genes(event(), GeneIndex([gene("A", "1", "+")]))
    assert pairs == [] and audit["pairing_status"] == "intergenic"
    pairs, audit = pair_genes(event(), GeneIndex([gene("A", "1", "+"), gene("A", "2", "+")]))
    assert pairs == [] and audit["pairing_status"] == "no_sense_pair"


def test_overlapping_gene_pairs_are_all_retained_as_ambiguous():
    pairs, audit = pair_genes(
        event(), GeneIndex([gene("A", "1", "+"), gene("C", "1", "+"), gene("B", "2", "+")])
    )
    assert len(pairs) == 2
    assert {r["Site1_Hugo_Symbol"] for r in pairs} == {"A", "C"}
    assert all(r["gene_pairing_status"] == "ambiguous" for r in pairs)
    assert audit["pairing_status"] == "ambiguous"


def test_half_open_gene_boundary_and_bins():
    index = GeneIndex([gene("A", "1", "+", 99_999, 100_001)])
    assert not index.overlap("1", 99_998, 99_999)
    assert len(index.overlap("1", 99_999, 100_001)) == 1
    assert not index.overlap("1", 100_001, 100_002)


def make_archive(path, *, link=False):
    with tarfile.open(path, "w:gz") as archive:
        member = tarfile.TarInfo("tcga/open/aliquot.pcawg.bedpe.gz")
        if link:
            member.type = tarfile.SYMTYPE
            member.linkname = "/tmp/unsafe"
            archive.addfile(member)
        else:
            payload = gzip.compress(bedpe().getvalue().encode())
            member.size = len(payload)
            archive.addfile(member, io.BytesIO(payload))


def test_archive_mapping_and_counts(tmp_path):
    path = tmp_path / "test.tgz"
    make_archive(path)
    rows, summary = read_archive(path, {"aliquot": SAMPLE})
    assert summary["sv_rows"] == summary["donors"] == summary["bedpe_files"] == 1
    assert rows[0]["icgc_specimen_id"] == "SP1"
    assert rows[0]["donor_wgs_exclusion_white_gray"] == "Whitelist"
    assert sorted(p.name for p in tmp_path.iterdir()) == ["test.tgz"]


@pytest.mark.parametrize("mode", ["link", "unmapped", "oversize", "nested_oversize", "not_tcga"])
def test_archive_fails_closed(tmp_path, monkeypatch, mode):
    import cfh.ingestion.pcawg_sv as parser

    path = tmp_path / "test.tgz"
    make_archive(path, link=mode == "link")
    samples = {"aliquot": dict(SAMPLE)}
    if mode == "unmapped":
        samples = {}
    elif mode == "oversize":
        monkeypatch.setattr(parser, "MAX_UNPACKED", 1)
    elif mode == "nested_oversize":
        # Uncompressed payload is larger than the tar member's gzip payload.
        payload = bedpe().getvalue() + "".join(
            "\t".join([*BASE[:6], f"SVMERGE{i}", *BASE[7:]]) + "\n" for i in range(2, 1000)
        )
        with tarfile.open(path, "w:gz") as archive:
            zipped = gzip.compress(payload.encode())
            member = tarfile.TarInfo("aliquot.pcawg.bedpe.gz")
            member.size = len(zipped)
            archive.addfile(member, io.BytesIO(zipped))
        monkeypatch.setattr(parser, "MAX_UNPACKED", len(zipped) + 10)
    elif mode == "not_tcga":
        samples["aliquot"]["dcc_project_code"] = "BRCA-UK"
    with pytest.raises(ValueError):
        read_archive(path, samples)


def test_sample_sheet_validation(tmp_path):
    path = tmp_path / "samples.tsv"
    keys = list(SAMPLE)
    path.write_text("\t".join(keys) + "\n" + "\t".join(SAMPLE[k] for k in keys) + "\n")
    assert sample_index(path)["aliquot"]["icgc_donor_id"] == "DO1"
    path.write_text(path.read_text().replace("Whitelist", "Maybe"))
    with pytest.raises(ValueError, match="unknown status"):
        sample_index(path)


def test_gencode_transcript_union(tmp_path: Path):
    attrs, comp = tmp_path / "attrs.gz", tmp_path / "comp.gz"
    with gzip.open(attrs, "wt") as handle:
        handle.write("ENSG1\tA\tprotein_coding\t\tENST1\n")
        handle.write("ENSG1\tA\tprotein_coding\t\tENST2\n")
    with gzip.open(comp, "wt") as handle:
        for transcript, start, end in [("ENST1", 10, 30), ("ENST2", 20, 40)]:
            handle.write(
                f"0\t{transcript}\tchr1\t+\t{start}\t{end}\t10\t30\t1\t10,\t30,\t0\tA\tcmpl\tcmpl\t0,\n"
            )
    index = read_gencode(comp, attrs)
    assert [(g.start, g.end) for g in index.overlap("1", 35, 36)] == [(10, 40)]


def test_header_only_bedpe_is_distinct_from_missing_file(tmp_path):
    path = tmp_path / "empty.tgz"
    payload = gzip.compress(("\t".join(HEADER) + "\n").encode())
    with tarfile.open(path, "w:gz") as archive:
        member = tarfile.TarInfo("aliquot.pcawg.bedpe.gz")
        member.size = len(payload)
        archive.addfile(member, io.BytesIO(payload))
    rows, summary = read_archive(path, {"aliquot": SAMPLE})
    assert rows == []
    assert summary["donors"] == 1
    assert summary["empty_bedpe_files"] == 1
    assert summary["file_row_counts"] == {"aliquot": 0}


@pytest.fixture
def screen_driver():
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "pcawg_driver", Path(__file__).parents[1] / "scripts" / "pcawg_tcga_sv_screen.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_successful_batch_omission_is_unassessed_and_cached(tmp_path, monkeypatch, screen_driver):
    from types import SimpleNamespace

    monkeypatch.setattr(screen_driver, "save_summary", lambda **kwargs: None)
    monkeypatch.setattr(
        screen_driver,
        "batch_fetch_canonical_transcripts",
        lambda *args, **kwargs: SimpleNamespace(by_gene_symbol={}),
    )
    source = screen_driver.PcawgTranscripts("GRCh37", cache_dir=tmp_path)
    source.prefetch_canonical(["MISSING"])
    assert source.canonical("MISSING") is None
    assert source._not_found_marker("MISSING").exists()
    assert "GRCh37" in str(source._not_found_marker("MISSING"))


def test_failed_batch_is_unassessed_but_not_persisted_as_notfound(
    tmp_path, monkeypatch, screen_driver
):
    import requests

    summaries = []
    monkeypatch.setattr(screen_driver, "save_summary", lambda **kwargs: summaries.append(kwargs))

    def fail(*args, **kwargs):
        raise requests.ConnectionError("lookup unavailable")

    monkeypatch.setattr(screen_driver, "batch_fetch_canonical_transcripts", fail)
    source = screen_driver.PcawgTranscripts("GRCh37", cache_dir=tmp_path)
    source.prefetch_canonical(["MISSING"])
    assert source.canonical("MISSING") is None
    assert not source._not_found_marker("MISSING").exists()
    assert summaries[-1]["canonical_lookup_errors"] == {"MISSING": "lookup unavailable"}


def test_missing_reference_and_homology_never_pass(tmp_path, monkeypatch, screen_driver):
    import requests

    from cfh.cohort.artefact_flags import EnsemblSequenceSource, genomic_flags, paralog_flag

    def fail(*args, **kwargs):
        raise requests.ConnectionError("lookup unavailable")

    monkeypatch.setattr(EnsemblSequenceSource, "_get", fail)
    source = screen_driver.SafeSequences("GRCh37", tmp_path, max_retries=0)
    row = {
        "junction_context_nt": "ACGT" * 25,
        "five_prime_gene": "A",
        "three_prime_gene": "B",
        "five_prime_chromosome": "1",
        "three_prime_chromosome": "2",
        "five_prime_breakpoint": 1000,
        "three_prime_breakpoint": 2000,
    }
    assert paralog_flag(row, source)["paralog_sequence_flag"] == "not_assessed"
    flags = genomic_flags(row, source, None, 1, "no tracks")
    assert flags["genomic_similarity_flag"] == "not_assessed"
    assert flags["read_evidence_flag"] == "not_assessed"
    assert flags["breakpoint_repeat_flag"] == "not_assessed"


@pytest.mark.parametrize(
    "mode", ["matches", "gene_id_mismatch", "strand_mismatch", "missing_mapping"]
)
def test_canonical_gene_and_strand_must_match_gencode(tmp_path, screen_driver, mode):
    import json
    from types import SimpleNamespace

    directory = tmp_path / "cache"
    directory.mkdir()
    for symbol in ("A", "B"):
        (directory / f"{symbol}.json").write_text(json.dumps({"geneId": f"ENSG_{symbol}"}))
    row = {
        "five_prime_gene": "A",
        "three_prime_gene": "B",
        "five_gene_id": "ENSG_A.1",
        "three_gene_id": "ENSG_B.2",
        "five_gene_strand": "+",
        "three_gene_strand": "-",
    }
    strand = -1 if mode == "strand_mismatch" else 1
    models = {
        "A": SimpleNamespace(exons=[SimpleNamespace(strand=strand)]),
        "B": SimpleNamespace(exons=[SimpleNamespace(strand=-1)]),
    }
    if mode == "gene_id_mismatch":
        row["five_gene_id"] = "ENSG_WRONG"
    if mode == "missing_mapping":
        models["A"] = None
    source = SimpleNamespace(canonical=models.get, _canonical_cache_dir=directory)
    problem = screen_driver.annotation_problem(row, source)
    assert (problem is None) == (mode == "matches")


@pytest.mark.parametrize(
    ("known_count", "unknown", "expected"),
    [
        (1, False, "pass"),
        (1, True, "not_assessed"),
        (2, False, "flagged"),
        (2, True, "flagged"),
    ],
)
def test_donor_burden_counts_unmodelled_shared_calls_once_and_fails_closed(
    screen_driver, known_count, unknown, expected
):
    rows = [
        {
            "donor_unique_id": "D",
            "aliquot_id": "A",
            "sv_id": f"SV{i}",
            "pfam_family": family,
            "reconstitution_class": "no_transcript_model",
        }
        for i in range(known_count)
        for family in ("PF00001", "PF00002")
    ]
    if unknown:
        rows.append(
            {
                "donor_unique_id": "D",
                "aliquot_id": "A",
                "sv_id": "missing",
                "pfam_family": "",
                "reconstitution_class": "no_transcript_model",
            }
        )
    result = screen_driver.donor_burdens(rows)["D"]
    assert result["sample_shared_family_calls"] == known_count
    assert result["sample_burden_flag"] == expected
    assert result["donor_family_lookup_status"] == ("not_assessed" if unknown else "assessed")
