import csv
import json

import pytest
from click.testing import CliRunner

from cfh.cli import main
from cfh.cohort.domain_screen import screen_cohort
from cfh.ingestion.open_fusion_sv_adapter import OPEN_PASSTHROUGH_COLUMNS, open_table_sv_rows
from cfh.ingestion.open_fusion_tables import SCHEMA
from cfh.mapping.junction_domain_classifier import CLEAN_RECONSTITUTION
from test_junction_translation import (
    BUILD,
    OTHER_FAMILY,
    _screen_frame,
    _synthetic_source,
)


def _open_row(**overrides):
    frame = _screen_frame().iloc[0]
    row = {column: "" for column in SCHEMA}
    row.update(
        cohort="TARGET",
        source_table="STAR-Fusion.v1.10.0.TARGETdb.tsv.gz",
        sample_id="TARGET-1-A",
        patient_id="TARGET-1",
        tissue_class="tumor",
        tumor_type="AML",
        gene_5p="AAA",
        gene_3p="BBB",
        chrom_5p="chr1",
        pos_5p=str(frame["Site1_Position"]),
        strand_5p="+",
        chrom_3p="chr2",
        pos_3p=str(frame["Site2_Position"]),
        strand_3p="+",
        genome_build=BUILD,
        junction_reads="12",
        spanning_reads="1e+01",
        ffpm="0.5",
        splice_type="ONLY_REF_SPLICE",
        annotation_flags="{}",
        caller="STAR-Fusion",
    )
    row.update(overrides)
    return row


def _write(path, rows, header=SCHEMA):
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for row in rows:
            writer.writerow([row[column] for column in SCHEMA])
    return path


def test_adapter_maps_open_rows_onto_structural_variant_fields(tmp_path):
    frame = open_table_sv_rows(_write(tmp_path / "t.tsv", [_open_row()]))
    row = frame.iloc[0].to_dict()
    assert row["Site1_Hugo_Symbol"] == "AAA" and row["Site2_Hugo_Symbol"] == "BBB"
    assert row["Event_Info"] == "{AAA:BBB} fusion"
    assert row["Site1_Chromosome"] == "chr1" and row["Site2_Chromosome"] == "chr2"
    assert row["NCBI_Build"] == BUILD
    assert row["Tumor_Split_Read_Count"] == 12
    assert row["Tumor_Paired_End_Read_Count"] == 10
    assert row["Breakpoint_Type"] is None and row["Annotation"] is None
    assert row["Open_tissue_class"] == "tumor"
    assert set(OPEN_PASSTHROUGH_COLUMNS) <= set(frame.columns)


def test_adapter_keeps_rows_without_genes_and_rejects_bad_input(tmp_path):
    frame = open_table_sv_rows(
        _write(tmp_path / "t.tsv", [_open_row(gene_5p="", gene_3p="", pos_5p="", pos_3p="")])
    )
    assert len(frame) == 1
    assert frame.iloc[0]["Event_Info"] == "fusion (genes unavailable)"
    with pytest.raises(ValueError, match="unexpected header"):
        open_table_sv_rows(_write(tmp_path / "h.tsv", [], header=SCHEMA[:-1]))
    with pytest.raises(ValueError, match="nonnegative integer"):
        open_table_sv_rows(_write(tmp_path / "c.tsv", [_open_row(junction_reads="1.5")]))


def test_open_rows_screen_like_the_equivalent_structural_variant_row(tmp_path):
    frame = open_table_sv_rows(_write(tmp_path / "t.tsv", [_open_row()]))
    rows, summary = screen_cohort(
        frame,
        _synthetic_source(),
        cohort="TARGET",
        passthrough_columns=OPEN_PASSTHROUGH_COLUMNS,
    )
    row = next(item for item in rows if item["pfam_family"] == OTHER_FAMILY)
    assert row["reconstitution_class"] == CLEAN_RECONSTITUTION
    assert row["fusion_name"] == "AAA::BBB"
    assert (row["five_prime_chromosome"], row["three_prime_chromosome"]) == ("chr1", "chr2")
    assert row["tumor_split_reads"] == 12 and row["tumor_paired_end_reads"] == 10
    assert row["Open_source_table"] == "STAR-Fusion.v1.10.0.TARGETdb.tsv.gz"
    assert summary["screened_fusion_events"] == 1


def test_open_rows_with_one_gene_or_the_same_gene_are_not_screened(tmp_path):
    frame = open_table_sv_rows(
        _write(tmp_path / "t.tsv", [_open_row(gene_3p=""), _open_row(gene_3p="AAA")])
    )
    rows, summary = screen_cohort(frame, _synthetic_source(), cohort="TARGET")
    assert rows == [] and summary["screened_fusion_events"] == 0


def test_cli_reads_an_open_fusion_table(tmp_path, monkeypatch):
    import cfh.cohort.domain_screen as domain_screen

    table = _write(tmp_path / "t.tsv", [_open_row()])
    seen = {}

    def fake_screen(raw_rows, source, **kwargs):
        seen["rows"] = raw_rows
        seen["kwargs"] = kwargs
        return [], {"class_counts": {}}

    monkeypatch.setattr(domain_screen, "screen_cohort", fake_screen)
    result = CliRunner().invoke(
        main,
        [
            "domain-screen",
            str(table),
            "--cohort",
            "TARGET",
            "--genome-build",
            "GRCh38",
            "--input-format",
            "open-fusion-table",
            "--omit-unshared-rows",
            "--output-dir",
            str(tmp_path / "out"),
        ],
    )
    assert result.exit_code == 0, result.output
    assert seen["rows"].iloc[0]["Site1_Hugo_Symbol"] == "AAA"
    assert seen["kwargs"]["keep_unshared_rows"] is False
    assert seen["kwargs"]["passthrough_columns"] == OPEN_PASSTHROUGH_COLUMNS
    with (tmp_path / "out" / "domain_screen.tsv").open() as handle:
        header = handle.readline().rstrip("\n").split("\t")
    assert header[-len(OPEN_PASSTHROUGH_COLUMNS) :] == list(OPEN_PASSTHROUGH_COLUMNS)
    assert json.loads((tmp_path / "out" / "domain_screen_summary.json").read_text()) == {
        "class_counts": {}
    }
