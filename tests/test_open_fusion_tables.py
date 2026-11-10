"""Offline fixtures use the actual pinned source headers, independently of the readers."""

import csv
import gzip
import hashlib
import io
import json
from dataclasses import replace
from unittest.mock import patch

import pytest

from cfh.ingestion.open_fusion_manifest import MANIFEST
from cfh.ingestion.open_fusion_tables import SCHEMA, fetch_table, normalize_table, read_table

FIXTURE_HEADERS = {
    "tcga_gtex": "data_class tissue_type tissue_subtype tumor_subtype sample fusion_name Jcount "
    "Scount Tcount SpliceType LeftGene LeftBreakpoint RightGene RightBreakpoint LargeAnchorSupport "
    "FFPM LeftBreakDinuc LeftBreakEntropy RightBreakDinuc RightBreakEntropy Annots "
    "PROT_FUSION_TYPE "
    "tumor_or_normal is_primary cosmic structure_type",
    "target_starfusion": "sample fusion_name JunctionReadCount SpanningFragCount est_J est_S "
    "SpliceType LeftGene LeftBreakpoint RightGene RightBreakpoint LargeAnchorSupport FFPM "
    "LeftBreakDinuc LeftBreakEntropy RightBreakDinuc RightBreakEntropy annots",
    "target_consolidated": "sample fusion_name LeftBreakpoint RightBreakpoint "
    "JunctionReadCount.starF SpanningFragCount.starF FFPM.starF partnerA partnerB "
    "sorted_fusion_name "
    "project participant core_participant_name sample_type tumor_or_normal Cancer.Type "
    "structure_type "
    "annots JunctionReadCount.FI SpanningFragCount.FI FFPM.FI fusion_cluster_att leiden "
    "FI_validated "
    "cosmic contains_ALL_driver_gene is_primary",
    "depmap_22q2": "DepMap_ID FusionName JunctionReadCount SpanningFragCount SpliceType LeftGene "
    "LeftBreakpoint RightGene RightBreakpoint LargeAnchorSupport FFPM LeftBreakDinuc "
    "LeftBreakEntropy "
    "RightBreakDinuc RightBreakEntropy annots CCLE_count",
}
EXPECTED_SCHEMA = (
    "cohort source_table sample_id patient_id tissue_class tumor_type gene_5p gene_3p "
    "transcript_5p "
    "transcript_3p chrom_5p pos_5p strand_5p chrom_3p pos_3p strand_3p genome_build junction_reads "
    "spanning_reads ffpm splice_type annotation_flags caller fi_validated filter_provenance"
).split()


def fixture_row(table):
    row = dict.fromkeys(FIXTURE_HEADERS[table].split(), "NA")
    row.update(
        {
            "LeftBreakpoint": "chr1:123:+",
            "RightBreakpoint": "chr2:456:-",
            "sample": "TARGET-10-ABCDEF-09A",
            "fusion_name": "ETV5::ETV1",
            "LeftGene": "ETV5^ENSG00000000001.1",
            "RightGene": "ETV1^ENSG00000000002.1",
            "JunctionReadCount": "0",
            "SpanningFragCount": "1",
            "FFPM": "0.001",
            "SpliceType": "INCL_NON_REF_SPLICE",
            "annots": '["RED_HERRING","PARALOG"]',
            "LargeAnchorSupport": "NO_LDAS",
        }
    )
    if table == "tcga_gtex":
        row.update(
            data_class="TCGA",
            sample="LUSC-51-4079-TP",
            tissue_type="LUSC",
            tumor_subtype="TP",
            tumor_or_normal="tumor",
            Jcount="0",
            Scount="1",
            Annots='["RED_HERRING","PARALOG"]',
        )
    if table == "target_starfusion":
        row["fusion_name"] = "ETV5--ETV1"
    if table == "target_consolidated":
        row.update(
            {
                "participant": "TARGET-10-ABCDEF",
                "sample_type": "Primary Tumor",
                "tumor_or_normal": "tumor",
                "Cancer.Type": "BALL",
                "partnerA": "ETV1",
                "partnerB": "ETV5",
                "FI_validated": "TRUE",
                "FFPM.FI": "0.002",
                "JunctionReadCount.FI": "2",
                "SpanningFragCount.FI": "3",
                "JunctionReadCount.starF": "0",
                "SpanningFragCount.starF": "1",
                "FFPM.starF": "0.001",
                "fusion_cluster_att": "Artifact",
            }
        )
    if table == "depmap_22q2":
        row.update(
            DepMap_ID="ACH-001113",
            FusionName="ETV5--ETV1",
            LeftGene="ETV5 (ENSG00000000001.1)",
            RightGene="ETV1 (ENSG00000000002.1)",
        )
    return {key: row[key] for key in FIXTURE_HEADERS[table].split()}


def write_fixture(path, table, rows, header=None):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=header or FIXTURE_HEADERS[table].split(),
            delimiter="," if table == "depmap_22q2" else "\t",
        )
        writer.writeheader()
        writer.writerows(rows)


def metadata(tmp_path, sample="TARGET-10-ABCDEF-09A", sample_type="Primary Tumor"):
    path = tmp_path / "metadata.tsv"
    path.write_text(
        "sample\tproject\tparticipant\tsample_type\tcore_participant_name\tCancer.Type\n"
        f"{sample}\tTARGET-ALL-P2\tTARGET-10-ABCDEF\t{sample_type}\tABCDEF\tBALL\n"
    )
    return path


@pytest.mark.parametrize("table", FIXTURE_HEADERS)
def test_exact_schema_and_no_filtering(tmp_path, table):
    source = tmp_path / MANIFEST[table].filename
    write_fixture(source, table, [fixture_row(table)] * 2)
    dest = tmp_path / "normalized.tsv"
    assert normalize_table(table, source, dest, target_sample_info=metadata(tmp_path)) == 2
    with dest.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        assert reader.fieldnames == EXPECTED_SCHEMA == list(SCHEMA)
        rows = list(reader)
    assert len(rows) == 2  # duplicates, low-read, same-family and flagged rows all survive
    row = rows[0]
    assert (row["gene_5p"], row["gene_3p"]) == ("ETV5", "ETV1")
    flags = json.loads(row["annotation_flags"])
    assert "RED_HERRING" in flags.get("Annots", flags.get("annots", ""))
    assert row["transcript_5p"] == row["transcript_3p"] == ""
    assert row["genome_build"] == "GRCh38"
    assert row["source_table"] == MANIFEST[table].filename
    assert "similarity_filter_ran" in row["filter_provenance"]
    assert row["caller"]
    if table == "depmap_22q2":
        assert row["tissue_class"] == "cell_line"
        assert row["patient_id"] == row["tumor_type"] == ""
    if table == "target_consolidated":
        assert row["junction_reads"] == "2"
        assert row["fi_validated"] == "TRUE"
        assert flags["JunctionReadCount.starF"] == "0"
        assert flags["fusion_cluster_att"] == "Artifact"


@pytest.mark.parametrize(
    "cohort,sample,tissue",
    [
        ("GTEx", "GTEX-14BMU-1526-SM-5TDE6", "normal"),
        ("TCGA", "LUSC-51-4079-NT", "normal"),
        ("TCGA", "LUSC-51-4079-TP", "tumor"),
    ],
)
def test_tcga_gtex_tissue_and_patient(tmp_path, cohort, sample, tissue):
    row = fixture_row("tcga_gtex")
    row.update(data_class=cohort, sample=sample, tumor_or_normal=tissue)
    if cohort == "GTEx":
        row.update(tissue_type="Colon", tissue_subtype="Colon - Transverse")
    source = tmp_path / "source.tsv"
    write_fixture(source, "tcga_gtex", [row])
    result = list(read_table("tcga_gtex", source))[0]
    assert result["tissue_class"] == tissue
    assert result["tumor_type"] == ("LUSC" if tissue == "tumor" else "")
    assert result["patient_id"] == ("TCGA-51-4079" if cohort == "TCGA" else "GTEX-14BMU")


@pytest.mark.parametrize("table", ["target_starfusion", "target_consolidated"])
def test_target_normal(tmp_path, table):
    row = fixture_row(table)
    row.update(sample="TARGET-10-ABCDEF-14A")
    if table == "target_consolidated":
        row.update(sample_type="Bone Marrow Normal", tumor_or_normal="normal")
    source = tmp_path / "target.tsv"
    write_fixture(source, table, [row])
    result = list(
        read_table(
            table,
            source,
            target_sample_info=metadata(tmp_path, row["sample"], "Bone Marrow Normal"),
        )
    )[0]
    assert result["tissue_class"] == "normal"
    assert result["patient_id"] == "TARGET-10-ABCDEF"
    assert result["tumor_type"] == ""


@pytest.mark.parametrize("table", FIXTURE_HEADERS)
@pytest.mark.parametrize("change", ["rename", "reorder", "extra", "missing"])
def test_changed_header_leaves_no_output(tmp_path, table, change):
    row = fixture_row(table)
    header = FIXTURE_HEADERS[table].split()
    if change == "rename":
        row["unknown"] = row.pop(header[0])
        header[0] = "unknown"
    elif change == "reorder":
        header[0], header[1] = header[1], header[0]
    elif change == "extra":
        header.append("unknown")
        row["unknown"] = "value"
    else:
        del row[header.pop()]
    source = tmp_path / "source.txt"
    write_fixture(source, table, [row], header)
    dest = tmp_path / "out.tsv"
    with pytest.raises(ValueError, match="header"):
        normalize_table(table, source, dest, target_sample_info=metadata(tmp_path))
    assert not dest.exists()


def test_late_bad_row_preserves_existing_destination(tmp_path):
    source = tmp_path / "source.tsv"
    good = fixture_row("tcga_gtex")
    bad = {**good, "RightBreakpoint": "broken"}
    write_fixture(source, "tcga_gtex", [good, bad])
    dest = tmp_path / "out.tsv"
    dest.write_text("existing output")
    with pytest.raises(ValueError, match="breakpoint"):
        normalize_table("tcga_gtex", source, dest)
    assert dest.read_text() == "existing output"
    assert sorted(p.name for p in tmp_path.iterdir()) == ["out.tsv", "source.tsv"]


def test_bad_field_count(tmp_path):
    source = tmp_path / "source.tsv"
    write_fixture(source, "tcga_gtex", [fixture_row("tcga_gtex")])
    with source.open("a") as f:
        f.write("too\tfew\n")
    with pytest.raises(ValueError, match="field count"):
        normalize_table("tcga_gtex", source, tmp_path / "out.tsv")
    assert not (tmp_path / "out.tsv").exists()


@pytest.mark.parametrize("missing", [True, False])
def test_target_missing_metadata_fails(tmp_path, missing):
    source = tmp_path / "source.tsv"
    write_fixture(source, "target_starfusion", [fixture_row("target_starfusion")])
    sidecar = None if missing else metadata(tmp_path, sample="different-sample")
    with pytest.raises(ValueError, match="TARGET"):
        normalize_table(
            "target_starfusion", source, tmp_path / "out.tsv", target_sample_info=sidecar
        )
    assert not (tmp_path / "out.tsv").exists()


def test_unknown_build_fails(tmp_path):
    with patch.dict(MANIFEST, {"tcga_gtex": replace(MANIFEST["tcga_gtex"], genome_build="")}):
        with pytest.raises(ValueError, match="genome build"):
            list(read_table("tcga_gtex", tmp_path / "unused.tsv"))


def test_checksum_mismatch(tmp_path):
    directory = tmp_path / "download"
    with patch("cfh.ingestion.open_fusion_tables.urlopen", return_value=io.BytesIO(b"corrupt")):
        with pytest.raises(ValueError, match="SHA-256 mismatch"):
            fetch_table("tcga_gtex", directory)
    assert list(directory.iterdir()) == []


def test_fetch_verified_and_no_directory_reuse(tmp_path):
    payload = b"synthetic bytes"
    spec = replace(MANIFEST["tcga_gtex"], sha256=hashlib.sha256(payload).hexdigest())
    directory = tmp_path / "download"
    with patch.dict(MANIFEST, {"tcga_gtex": spec}):
        with patch("cfh.ingestion.open_fusion_tables.urlopen", return_value=io.BytesIO(payload)):
            assert fetch_table("tcga_gtex", directory).read_bytes() == payload
        with pytest.raises(FileExistsError):
            fetch_table("tcga_gtex", directory)


def test_unvalidated_consolidated_uses_starfusion_counts(tmp_path):
    row = fixture_row("target_consolidated")
    row["FI_validated"] = "FALSE"
    source = tmp_path / "source.tsv"
    write_fixture(source, "target_consolidated", [row])
    result = list(read_table("target_consolidated", source))[0]
    assert result["junction_reads"] == "0"
    assert result["ffpm"] == "0.001"
    assert result["splice_type"] == ""


def test_unavailable_breakpoint_is_preserved(tmp_path):
    row = fixture_row("tcga_gtex")
    row["LeftBreakpoint"] = "NA"
    source = tmp_path / "source.tsv"
    write_fixture(source, "tcga_gtex", [row])
    result = list(read_table("tcga_gtex", source))[0]
    assert result["gene_5p"] == "ETV5"
    assert result["chrom_5p"] == result["pos_5p"] == result["strand_5p"] == ""
    assert result["pos_3p"] == "456"


def test_breakpoint_without_strand_keeps_known_coordinates(tmp_path):
    row = fixture_row("target_starfusion")
    row["LeftBreakpoint"] = "chr14:105856169"
    source = tmp_path / "source.tsv"
    write_fixture(source, "target_starfusion", [row])
    result = list(read_table("target_starfusion", source, target_sample_info=metadata(tmp_path)))[0]
    assert result["chrom_5p"] == "chr14"
    assert result["pos_5p"] == "105856169"
    assert result["strand_5p"] == ""


@pytest.mark.parametrize("sample", ["UNKNOWN-DONOR-1", "GTEX-", "GTEX", "K-562-UNKNOWN"])
def test_non_gtex_id_fails(tmp_path, sample):
    row = fixture_row("tcga_gtex")
    row.update(data_class="GTEx", sample=sample, tumor_or_normal="normal")
    source = tmp_path / "source.tsv"
    write_fixture(source, "tcga_gtex", [row])
    with pytest.raises(ValueError, match="GTEx sample identifier"):
        normalize_table("tcga_gtex", source, tmp_path / "out.tsv")
    assert not (tmp_path / "out.tsv").exists()


@pytest.mark.parametrize("sample", ["K-562", "K-562-SM-26GMQ"])
def test_verified_k562_is_cell_line(tmp_path, sample):
    row = fixture_row("tcga_gtex")
    row.update(data_class="GTEx", sample=sample, tumor_or_normal="normal")
    source = tmp_path / "source.tsv"
    write_fixture(source, "tcga_gtex", [row])
    result = list(read_table("tcga_gtex", source))[0]
    assert result["tissue_class"] == "cell_line"
    assert result["patient_id"] == result["tumor_type"] == ""
    assert "tissue_class_override=K-562_ATCC_CCL-243" in result["filter_provenance"]
    assert json.loads(result["annotation_flags"])["tumor_or_normal"] == "normal"


@pytest.mark.parametrize("field", ["Jcount", "Scount", "FFPM"])
@pytest.mark.parametrize("value", ["not-a-number", "NaN", "nan", "inf", "-inf", "1_000"])
def test_malformed_measurement_fails(tmp_path, field, value):
    row = fixture_row("tcga_gtex")
    row[field] = value
    source = tmp_path / "source.tsv"
    write_fixture(source, "tcga_gtex", [row])
    with pytest.raises(ValueError, match="junction_reads|spanning_reads|ffpm"):
        normalize_table("tcga_gtex", source, tmp_path / "out.tsv")
    assert not (tmp_path / "out.tsv").exists()


@pytest.mark.parametrize("field", ["Jcount", "Scount"])
def test_negative_counts_fail(tmp_path, field):
    row = fixture_row("tcga_gtex")
    row[field] = "-1e-5"
    source = tmp_path / "source.tsv"
    write_fixture(source, "tcga_gtex", [row])
    with pytest.raises(ValueError, match="Negative"):
        list(read_table("tcga_gtex", source))


@pytest.mark.parametrize("value", ["1e+05", "0", "1.0", "20e-1", "NA", ".", ""])
def test_numeric_formats_preserve_source(tmp_path, value):
    row = fixture_row("tcga_gtex")
    row.update(Jcount=value, Scount=value, FFPM=value)
    source = tmp_path / "source.tsv"
    write_fixture(source, "tcga_gtex", [row])
    result = list(read_table("tcga_gtex", source))[0]
    expected = "" if value in {"NA", ".", ""} else value
    for key in ["junction_reads", "spanning_reads", "ffpm"]:
        assert result[key] == expected
    assert json.loads(result["annotation_flags"])["Jcount"] == value


def test_ffpm_float_overflow_fails(tmp_path):
    row = fixture_row("tcga_gtex")
    row["FFPM"] = "1e999"
    source = tmp_path / "source.tsv"
    write_fixture(source, "tcga_gtex", [row])
    with pytest.raises(ValueError, match="Non-finite ffpm"):
        list(read_table("tcga_gtex", source))


def test_lf_output(tmp_path):
    source = tmp_path / "source.tsv"
    write_fixture(source, "tcga_gtex", [fixture_row("tcga_gtex")])
    dest = tmp_path / "out.tsv"
    normalize_table("tcga_gtex", source, dest)
    payload = dest.read_bytes()
    assert b"\r" not in payload
    assert payload.count(b"\n") == 2
    assert payload.endswith(b"\n")


def test_unknown_tcga_sample_code_raises(tmp_path):
    row = fixture_row("tcga_gtex")
    row["sample"] = "LUSC-51-4079-UNKNOWN"
    source = tmp_path / "source.tsv"
    write_fixture(source, "tcga_gtex", [row])
    with pytest.raises(ValueError, match="TCGA sample identifier"):
        list(read_table("tcga_gtex", source))


@pytest.mark.parametrize("table", ["target_starfusion", "target_consolidated"])
def test_unknown_target_sample_type_raises(tmp_path, table):
    row = fixture_row(table)
    if table == "target_consolidated":
        row["sample_type"] = "Unknown Cancer Material"
    source = tmp_path / "source.tsv"
    write_fixture(source, table, [row])
    with pytest.raises(ValueError, match="Unknown TARGET sample type"):
        list(
            read_table(
                table,
                source,
                target_sample_info=metadata(tmp_path, sample_type="Unknown Cancer Material"),
            )
        )


@pytest.mark.parametrize("table", ["tcga_gtex", "target_consolidated"])
@pytest.mark.parametrize("value", ["normal", "unknown"])
def test_tumor_or_normal_conflict_raises(tmp_path, table, value):
    row = fixture_row(table)
    row["tumor_or_normal"] = value
    source = tmp_path / "source.tsv"
    write_fixture(source, table, [row])
    with pytest.raises(ValueError, match="tissue metadata|tumor_or_normal"):
        list(read_table(table, source))


def test_tcga_disease_conflict_raises(tmp_path):
    row = fixture_row("tcga_gtex")
    row["tissue_type"] = "BRCA"
    source = tmp_path / "source.tsv"
    write_fixture(source, "tcga_gtex", [row])
    with pytest.raises(ValueError, match="TCGA disease"):
        list(read_table("tcga_gtex", source))


@pytest.mark.parametrize("table", ["tcga_gtex", "target_starfusion", "depmap_22q2"])
def test_oriented_gene_name_conflict_raises(tmp_path, table):
    row = fixture_row(table)
    row["LeftGene"], row["RightGene"] = row["RightGene"], row["LeftGene"]
    source = tmp_path / "source.tsv"
    write_fixture(source, table, [row])
    with pytest.raises(ValueError, match="oriented genes"):
        list(read_table(table, source, target_sample_info=metadata(tmp_path)))


def test_fully_unavailable_fusion_row_survives(tmp_path):
    row = fixture_row("tcga_gtex")
    for key in ["fusion_name", "LeftGene", "RightGene", "LeftBreakpoint", "RightBreakpoint"]:
        row[key] = "NA"
    source = tmp_path / "source.tsv"
    write_fixture(source, "tcga_gtex", [row])
    result = list(read_table("tcga_gtex", source))[0]
    assert result["gene_5p"] == result["gene_3p"] == ""
    assert result["pos_5p"] == result["pos_3p"] == ""


def test_depmap_provenance_is_config_derived():
    assert "provenance_basis=22Q2_saved_config" in MANIFEST["depmap_22q2"].filter_provenance


def test_target_filter_citation_matches_caller_version():
    for table in ["target_starfusion", "target_consolidated"]:
        assert "v1.10.0/util/STAR-Fusion.filter" in MANIFEST[table].filter_evidence


@pytest.mark.parametrize(
    "subtype", ["Cells - EBV-transformed lymphocytes", "Cells - Transformed fibroblasts"]
)
def test_gtex_transformed_cells_are_cell_lines(tmp_path, subtype):
    row = fixture_row("tcga_gtex")
    row.update(
        data_class="GTEx",
        sample="GTEX-14BMU-1526-SM-5TDE6",
        tumor_or_normal="normal",
        tissue_subtype=subtype,
    )
    source = tmp_path / "source.tsv"
    write_fixture(source, "tcga_gtex", [row])
    result = list(read_table("tcga_gtex", source))[0]
    assert result["tissue_class"] == "cell_line"
    assert result["patient_id"] == "GTEX-14BMU"
    assert result["tumor_type"] == ""


@pytest.mark.parametrize("table", FIXTURE_HEADERS)
def test_wrong_fusion_delimiter_raises(tmp_path, table):
    row = fixture_row(table)
    key = "FusionName" if table == "depmap_22q2" else "fusion_name"
    row[key] = row[key].replace("::", "|").replace("--", "|")
    source = tmp_path / "source.tsv"
    write_fixture(source, table, [row])
    with pytest.raises(ValueError, match="fusion name"):
        list(read_table(table, source, target_sample_info=metadata(tmp_path)))


def test_consolidated_partner_conflict_raises(tmp_path):
    row = fixture_row("target_consolidated")
    row["partnerA"] = "UNKNOWN"
    source = tmp_path / "source.tsv"
    write_fixture(source, "target_consolidated", [row])
    with pytest.raises(ValueError, match="partner-identity"):
        list(read_table("target_consolidated", source))


def test_consolidated_distinct_names_at_same_coordinates_survive(tmp_path):
    first = fixture_row("target_consolidated")
    second = {**first, "fusion_name": "ETV5::ETV4", "partnerA": "ETV4", "partnerB": "ETV5"}
    source = tmp_path / "source.tsv"
    write_fixture(source, "target_consolidated", [first, second])
    result = list(read_table("target_consolidated", source))
    assert [row["gene_3p"] for row in result] == ["ETV1", "ETV4"]


@pytest.mark.parametrize("column", ["tissue_type", "tissue_subtype"])
@pytest.mark.parametrize(
    "label",
    [
        "Cells - EBV-transformed lymphocytes",
        "Cells - Transformed fibroblasts",
        "Cells - Cultured fibroblasts",
        "cElLs - cUlTuReD fIbRoBlAsTs",
        "CELLS_-_EBV-transformed_lymphocytes",
        "cells_transformed_fibroblasts",
        "cells_cultured_fibroblasts",
    ],
)
def test_explicit_gtex_cell_label_variants(tmp_path, column, label):
    row = fixture_row("tcga_gtex")
    row.update(
        data_class="GTEx",
        sample="GTEX-14BMU-1526-SM-5TDE6",
        tumor_or_normal="normal",
        tissue_type="Skin",
        tissue_subtype="Skin",
    )
    row[column] = label
    source = tmp_path / "source.tsv"
    write_fixture(source, "tcga_gtex", [row])
    result = list(read_table("tcga_gtex", source))[0]
    assert result["tissue_class"] == "cell_line"
    assert result["tumor_type"] == ""
    assert result["patient_id"] == "GTEX-14BMU"
    assert json.loads(result["annotation_flags"])[column] == label


@pytest.mark.parametrize("column", ["tissue_type", "tissue_subtype"])
@pytest.mark.parametrize("sample", ["GTEX-14BMU-1526-SM-5TDE6", "K-562-SM-26GMQ"])
@pytest.mark.parametrize("label", ["Cells - Unknown lineage", "cELls_unknown", "CellsUnknown"])
def test_unknown_cells_labels_never_default_to_normal(tmp_path, column, sample, label):
    row = fixture_row("tcga_gtex")
    row.update(
        data_class="GTEx",
        sample=sample,
        tumor_or_normal="normal",
        tissue_type="Skin",
        tissue_subtype="Cells - Cultured fibroblasts",
    )
    row[column] = label  # unknown in either column must win over a known cell label
    source = tmp_path / "source.tsv"
    write_fixture(source, "tcga_gtex", [row])
    with pytest.raises(ValueError, match="Unknown GTEx Cells label"):
        normalize_table("tcga_gtex", source, tmp_path / "out.tsv")
    assert not (tmp_path / "out.tsv").exists()


def test_unknown_noncell_gtex_pair_requires_mapping_update(tmp_path):
    row = fixture_row("tcga_gtex")
    row.update(
        data_class="GTEx",
        sample="GTEX-14BMU-1526-SM-5TDE6",
        tumor_or_normal="normal",
        tissue_type="Unknown",
        tissue_subtype="Unknown",
    )
    source = tmp_path / "source.tsv"
    write_fixture(source, "tcga_gtex", [row])
    with pytest.raises(ValueError, match="Unknown GTEx tissue pair"):
        list(read_table("tcga_gtex", source))


@pytest.mark.parametrize(
    "sample_type",
    [
        "Cell Line",
        "Cell Lines",
        "Primary Xenograft Tissue",
        "Cell Line Derived Xenograft Tissue",
        "cell_line",
        "CELL_LINES",
        "Primary Tumor Xenograft",
        "XENOGRAFT_TISSUE",
    ],
)
@pytest.mark.parametrize("table", ["target_starfusion", "target_consolidated"])
def test_target_cell_lines_and_xenografts_never_become_tumor(tmp_path, table, sample_type):
    row = fixture_row(table)
    if table == "target_consolidated":
        row["sample_type"] = sample_type
    source = tmp_path / "source.tsv"
    write_fixture(source, table, [row])
    result = list(
        read_table(table, source, target_sample_info=metadata(tmp_path, sample_type=sample_type))
    )[0]
    assert result["tissue_class"] == "cell_line"
    assert result["tumor_type"] == ""
    assert "tissue_class_override=TARGET_cell_line_or_xenograft" in result["filter_provenance"]


@pytest.mark.parametrize("value", ["-1", "-0.5", "-1e-999"])
def test_negative_ffpm_rejected_before_float_underflow(tmp_path, value):
    row = fixture_row("tcga_gtex")
    row["FFPM"] = value
    source = tmp_path / "source.tsv"
    write_fixture(source, "tcga_gtex", [row])
    with pytest.raises(ValueError, match="Negative ffpm"):
        list(read_table("tcga_gtex", source))


@pytest.mark.parametrize("field", ["Jcount", "Scount"])
@pytest.mark.parametrize("value", ["1.25", ".5", "1e-5"])
def test_fractional_reads_rejected(tmp_path, field, value):
    row = fixture_row("tcga_gtex")
    row[field] = value
    source = tmp_path / "source.tsv"
    write_fixture(source, "tcga_gtex", [row])
    with pytest.raises(ValueError, match="Non-integer"):
        list(read_table("tcga_gtex", source))


@pytest.mark.parametrize("value", ["1.25", ".5", "1e-5"])
def test_fractional_ffpm_is_valid(tmp_path, value):
    row = fixture_row("tcga_gtex")
    row["FFPM"] = value
    source = tmp_path / "source.tsv"
    write_fixture(source, "tcga_gtex", [row])
    assert list(read_table("tcga_gtex", source))[0]["ffpm"] == value


def test_gene_conflict_reports_identity_without_dumping_row(tmp_path):
    row = fixture_row("tcga_gtex")
    row.update(LeftGene="WRONG^ENSG123", Annots="SENTINEL_NOT_FOR_EXCEPTION")
    source = tmp_path / "source.tsv"
    write_fixture(source, "tcga_gtex", [row])
    with pytest.raises(ValueError) as error:
        list(read_table("tcga_gtex", source))
    message = str(error.value)
    assert row["sample"] in message
    assert row["fusion_name"] in message
    assert "SENTINEL_NOT_FOR_EXCEPTION" not in message


@pytest.mark.parametrize("table", ["target_starfusion", "target_consolidated", "depmap_22q2"])
def test_na_name_placeholder_not_accepted_outside_tcga_gtex(tmp_path, table):
    row = fixture_row(table)
    row["FusionName" if table == "depmap_22q2" else "fusion_name"] = "NA"
    source = tmp_path / "source.tsv"
    write_fixture(source, table, [row])
    with pytest.raises(ValueError, match="Unavailable fusion name"):
        list(read_table(table, source, target_sample_info=metadata(tmp_path)))
