import csv
import json

import pandas as pd
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


# --- source-format round trips: raw table -> normalize_table -> SV rows -----------

# Minus-strand junction coordinates on chr21 (both '-') for the fixture's ETV5::ETV1 row,
# with each source's own read-count columns.
_MINUS_STRAND = {"LeftBreakpoint": "chr21:41508081:-", "RightBreakpoint": "chr21:38445621:-"}
_SOURCE_COUNTS = {
    # table: (raw count fields, expected junction reads, expected spanning reads)
    "tcga_gtex": ({"Jcount": "17", "Scount": "4"}, 17, 4),
    "target_starfusion": ({"JunctionReadCount": "17", "SpanningFragCount": "4"}, 17, 4),
    # FusionInspector-validated rows report the FI counts, not the STAR-Fusion ones
    "target_consolidated": (
        {
            "JunctionReadCount.FI": "17",
            "SpanningFragCount.FI": "4",
            "JunctionReadCount.starF": "9",
            "SpanningFragCount.starF": "1",
            "FI_validated": "TRUE",
        },
        17,
        4,
    ),
    "depmap_22q2": ({"JunctionReadCount": "17", "SpanningFragCount": "4"}, 17, 4),
}


@pytest.mark.parametrize("table", list(_SOURCE_COUNTS))
def test_minus_strand_source_rows_keep_their_coordinates_through_the_adapter(tmp_path, table):
    from cfh.ingestion.open_fusion_manifest import MANIFEST
    from cfh.ingestion.open_fusion_tables import normalize_table
    from test_open_fusion_tables import fixture_row, metadata, write_fixture

    fields, junction, spanning = _SOURCE_COUNTS[table]
    raw = fixture_row(table)
    raw.update(_MINUS_STRAND)
    raw.update(fields)
    source = tmp_path / MANIFEST[table].filename
    write_fixture(source, table, [raw])
    normalized = tmp_path / "normalized.tsv"
    normalize_table(table, source, normalized, target_sample_info=metadata(tmp_path))

    row = open_table_sv_rows(normalized).iloc[0].to_dict()
    # 5'/3' come from the fusion name, never from coordinate order or strand
    assert (row["Site1_Hugo_Symbol"], row["Site2_Hugo_Symbol"]) == ("ETV5", "ETV1")
    assert row["Event_Info"] == "{ETV5:ETV1} fusion"
    # STAR-Fusion's 1-based junction bases are copied verbatim: no strand flip, no
    # off-by-one, and the 5' position may be larger than the 3' one on the minus strand
    assert (row["Site1_Chromosome"], row["Site1_Position"]) == ("chr21", 41508081)
    assert (row["Site2_Chromosome"], row["Site2_Position"]) == ("chr21", 38445621)
    assert row["Open_strand_5p"] == row["Open_strand_3p"] == "-"
    assert row["Tumor_Split_Read_Count"] == junction
    assert row["Tumor_Paired_End_Read_Count"] == spanning
    assert row["NCBI_Build"] == "GRCh38"


def test_mixed_strand_and_missing_breakpoint_rows_are_kept_as_reported(tmp_path):
    rows = [
        _open_row(strand_5p="+", strand_3p="-", pos_5p="58252506", pos_3p="58275702"),
        _open_row(chrom_5p="", pos_5p="", strand_5p=""),
    ]
    frame = open_table_sv_rows(_write(tmp_path / "t.tsv", rows))
    first, second = (frame.iloc[i].to_dict() for i in range(2))
    assert (first["Site1_Position"], first["Site2_Position"]) == (58252506, 58275702)
    assert (first["Open_strand_5p"], first["Open_strand_3p"]) == ("+", "-")
    # a missing 5' breakpoint stays missing (never 0); the 3' side is untouched
    assert pd.isna(second["Site1_Chromosome"]) and pd.isna(second["Site1_Position"])
    assert second["Site2_Chromosome"] == "chr2"
    assert second["Site2_Position"] == int(_open_row()["pos_3p"])
