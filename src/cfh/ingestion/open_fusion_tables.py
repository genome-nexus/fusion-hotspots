"""Strict, lossless-row readers for pinned open fusion candidate tables."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import os
import re
import tempfile
from collections.abc import Iterator
from decimal import Decimal
from pathlib import Path
from urllib.request import urlopen

from cfh.ingestion.open_fusion_manifest import MANIFEST, TARGET_SAMPLE_INFO, SourceTable

SCHEMA = tuple(
    "cohort source_table sample_id patient_id tissue_class tumor_type gene_5p gene_3p "
    "transcript_5p transcript_3p chrom_5p pos_5p strand_5p chrom_3p pos_3p strand_3p "
    "genome_build junction_reads spanning_reads ffpm splice_type annotation_flags caller "
    "fi_validated filter_provenance".split()
)
HEADERS = {
    "tcga_gtex": tuple(
        "data_class tissue_type tissue_subtype tumor_subtype sample fusion_name Jcount Scount "
        "Tcount SpliceType LeftGene LeftBreakpoint RightGene RightBreakpoint LargeAnchorSupport "
        "FFPM LeftBreakDinuc LeftBreakEntropy RightBreakDinuc RightBreakEntropy Annots "
        "PROT_FUSION_TYPE tumor_or_normal is_primary cosmic structure_type".split()
    ),
    "target_starfusion": tuple(
        "sample fusion_name JunctionReadCount SpanningFragCount est_J est_S SpliceType LeftGene "
        "LeftBreakpoint RightGene RightBreakpoint LargeAnchorSupport FFPM LeftBreakDinuc "
        "LeftBreakEntropy RightBreakDinuc RightBreakEntropy annots".split()
    ),
    "target_consolidated": tuple(
        "sample fusion_name LeftBreakpoint RightBreakpoint JunctionReadCount.starF "
        "SpanningFragCount.starF FFPM.starF partnerA partnerB sorted_fusion_name project "
        "participant core_participant_name sample_type tumor_or_normal Cancer.Type structure_type "
        "annots JunctionReadCount.FI SpanningFragCount.FI FFPM.FI fusion_cluster_att leiden "
        "FI_validated cosmic contains_ALL_driver_gene is_primary".split()
    ),
    "depmap_22q2": tuple(
        "DepMap_ID FusionName JunctionReadCount SpanningFragCount SpliceType LeftGene "
        "LeftBreakpoint RightGene RightBreakpoint LargeAnchorSupport FFPM LeftBreakDinuc "
        "LeftBreakEntropy RightBreakDinuc RightBreakEntropy annots CCLE_count".split()
    ),
}
_METADATA_HEADER = (
    "sample",
    "project",
    "participant",
    "sample_type",
    "core_participant_name",
    "Cancer.Type",
)


# Observed pairs from the pinned TCGA+GTEx table, not an open-ended normal default.
_GTEX_TISSUE_CLASSES = {
    ("Adipose Tissue", "Adipose - Subcutaneous"): "normal",
    ("Adipose Tissue", "Adipose - Visceral (Omentum)"): "normal",
    ("Adrenal Gland", "Adrenal Gland"): "normal",
    ("Bladder", "Bladder"): "normal",
    ("Blood", "Cells - EBV-transformed lymphocytes"): "cell_line",
    ("Blood", "Whole Blood"): "normal",
    ("Blood Vessel", "Artery - Aorta"): "normal",
    ("Blood Vessel", "Artery - Coronary"): "normal",
    ("Blood Vessel", "Artery - Tibial"): "normal",
    ("Brain", "Brain - Amygdala"): "normal",
    ("Brain", "Brain - Anterior cingulate cortex (BA24)"): "normal",
    ("Brain", "Brain - Caudate (basal ganglia)"): "normal",
    ("Brain", "Brain - Cerebellar Hemisphere"): "normal",
    ("Brain", "Brain - Cerebellum"): "normal",
    ("Brain", "Brain - Cortex"): "normal",
    ("Brain", "Brain - Frontal Cortex (BA9)"): "normal",
    ("Brain", "Brain - Hippocampus"): "normal",
    ("Brain", "Brain - Hypothalamus"): "normal",
    ("Brain", "Brain - Nucleus accumbens (basal ganglia)"): "normal",
    ("Brain", "Brain - Putamen (basal ganglia)"): "normal",
    ("Brain", "Brain - Spinal cord (cervical c-1)"): "normal",
    ("Brain", "Brain - Substantia nigra"): "normal",
    ("Breast", "Breast - Mammary Tissue"): "normal",
    ("Cervix Uteri", "Cervix - Ectocervix"): "normal",
    ("Cervix Uteri", "Cervix - Endocervix"): "normal",
    ("Colon", "Colon - Sigmoid"): "normal",
    ("Colon", "Colon - Transverse"): "normal",
    ("Esophagus", "Esophagus - Gastroesophageal Junction"): "normal",
    ("Esophagus", "Esophagus - Mucosa"): "normal",
    ("Esophagus", "Esophagus - Muscularis"): "normal",
    ("Fallopian Tube", "Fallopian Tube"): "normal",
    ("Heart", "Heart - Atrial Appendage"): "normal",
    ("Heart", "Heart - Left Ventricle"): "normal",
    ("Kidney", "Kidney - Cortex"): "normal",
    ("Liver", "Liver"): "normal",
    ("Lung", "Lung"): "normal",
    ("Muscle", "Muscle - Skeletal"): "normal",
    ("Nerve", "Nerve - Tibial"): "normal",
    ("Ovary", "Ovary"): "normal",
    ("Pancreas", "Pancreas"): "normal",
    ("Pituitary", "Pituitary"): "normal",
    ("Prostate", "Prostate"): "normal",
    ("Salivary Gland", "Minor Salivary Gland"): "normal",
    ("Skin", "Cells - Transformed fibroblasts"): "cell_line",
    ("Skin", "Skin - Not Sun Exposed (Suprapubic)"): "normal",
    ("Skin", "Skin - Sun Exposed (Lower leg)"): "normal",
    ("Small Intestine", "Small Intestine - Terminal Ileum"): "normal",
    ("Spleen", "Spleen"): "normal",
    ("Stomach", "Stomach"): "normal",
    ("Testis", "Testis"): "normal",
    ("Thyroid", "Thyroid"): "normal",
    ("Uterus", "Uterus"): "normal",
    ("Vagina", "Vagina"): "normal",
}
# Explicit aliases: case and separators are normalized before exact matching.
_GTEX_CELL_LABELS = {
    "cells ebv transformed lymphocytes",
    "cells transformed fibroblasts",
    "cells cultured fibroblasts",
}


def _tissue_label(value: str) -> str:
    return " ".join(value.replace("_", " ").replace("-", " ").casefold().split())


_GTEX_CANONICAL_CLASSES = {
    (_tissue_label(kind), _tissue_label(subtype)): tissue
    for (kind, subtype), tissue in _GTEX_TISSUE_CLASSES.items()
}


def _gtex_labels(tissue_type: str, tissue_subtype: str) -> tuple[str, str]:
    labels = (_tissue_label(tissue_type), _tissue_label(tissue_subtype))
    cell_labels = [label for label in labels if label.startswith("cells")]
    # Check both columns before accepting either: a known label cannot hide an unknown one.
    for label in cell_labels:
        if label not in _GTEX_CELL_LABELS:
            raise ValueError(f"Unknown GTEx Cells label: {label!r}")
    return labels


def _gtex_tissue(tissue_type: str, tissue_subtype: str) -> str:
    labels = _gtex_labels(tissue_type, tissue_subtype)
    if any(label in _GTEX_CELL_LABELS for label in labels):
        return "cell_line"
    try:
        return _GTEX_CANONICAL_CLASSES[labels]
    except KeyError:
        raise ValueError(f"Unknown GTEx tissue pair: {labels!r}") from None


def fetch_table(table: str, directory: Path) -> Path:
    """Download into a newly created directory; never reuse or overwrite a download.

    ``target_sample_info`` fetches the required raw TARGET metadata sidecar.
    Failed downloads leave no file that could be mistaken for verified data.
    """
    spec = TARGET_SAMPLE_INFO if table == "target_sample_info" else MANIFEST[table]
    directory.mkdir(parents=True, exist_ok=False)
    target = directory / spec.filename
    partial = directory / (spec.filename + ".part")
    try:
        digest = hashlib.sha256()
        with urlopen(spec.url, timeout=120) as response, partial.open("xb") as output:
            while block := response.read(1024 * 1024):
                digest.update(block)
                output.write(block)
        if digest.hexdigest() != spec.sha256:
            raise ValueError(f"SHA-256 mismatch for {spec.filename}: {digest.hexdigest()}")
        partial.replace(target)
    finally:
        partial.unlink(missing_ok=True)
    return target


def _rows(path: Path, header: tuple[str, ...], delimiter: str) -> Iterator[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter=delimiter, strict=True)
        actual = next(reader, None)
        if actual != list(header):
            raise ValueError(f"Unknown or changed header in {path}: {actual!r}")
        for line, values in enumerate(reader, 2):
            if len(values) != len(header):
                raise ValueError(f"Wrong field count in {path}:{line}")
            yield dict(zip(header, values, strict=True))


def _available(value: str) -> str:
    return "" if value in {"NA", "NaN", "."} else value


def _target_metadata(path: Path) -> dict[str, dict[str, str]]:
    result = {}
    for row in _rows(path, _METADATA_HEADER, "\t"):
        if row["sample"] in result:
            raise ValueError(f"Duplicate TARGET metadata sample: {row['sample']}")
        result[row["sample"]] = row
    return result


def _target_tissue(sample_type: str) -> str:
    # Culture and xenograft samples must never enter the tumor/normal tissue panels.
    label = _tissue_label(sample_type)
    if re.search(r"\bcell lines?\b", label) or "xenograft" in label:
        return "cell_line"
    # Exact published sample types: an unfamiliar type must not silently become tumor.
    if sample_type in {"Bone Marrow Normal", "Solid Tissue Normal", "Blood Derived Normal"}:
        return "normal"
    if sample_type in {
        "Primary Tumor",
        "Recurrent Tumor",
        "Metastatic",
        "Primary Blood Derived Cancer - Bone Marrow",
        "Primary Blood Derived Cancer - Peripheral Blood",
        "Recurrent Blood Derived Cancer - Bone Marrow",
        "Recurrent Blood Derived Cancer - Peripheral Blood",
        "Blood Derived Cancer - Bone Marrow, Post-treatment",
        "Blood Derived Cancer - Peripheral Blood, Post-treatment",
    }:
        return "tumor"
    raise ValueError(f"Unknown TARGET sample type: {sample_type!r}")


def _number(value: str, field: str) -> str:
    # Only explicit unavailable tokens are missing; NaN/Inf are invalid measurements.
    if value in {"", "NA", "."}:
        return ""
    if not re.fullmatch(r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?", value):
        raise ValueError(f"Invalid {field}: {value!r}")
    number = Decimal(value)
    if number < 0:
        raise ValueError(f"Negative {field}: {value!r}")
    if field == "ffpm":
        if not math.isfinite(float(value)):
            raise ValueError(f"Non-finite {field}: {value!r}")
    elif number != number.to_integral_value():
        raise ValueError(f"Non-integer {field}: {value!r}")
    return value


def _fusion_genes(row: dict[str, str], table: str) -> tuple[str, str]:
    name = row["FusionName" if table == "depmap_22q2" else "fusion_name"]
    delimiter = "--" if table in {"target_starfusion", "depmap_22q2"} else "::"
    if not _available(name):
        if table == "tcga_gtex":
            return "", ""
        raise ValueError(f"Unavailable fusion name for {table}: {name!r}")
    genes = name.split(delimiter)
    if len(genes) != 2 or not all(genes):
        raise ValueError(f"Invalid oriented fusion name for {table}: {name!r}")
    return genes[0], genes[1]


def _normalize(
    row: dict[str, str], table: str, spec: SourceTable, metadata: dict[str, dict[str, str]]
) -> dict[str, str]:
    out = dict.fromkeys(SCHEMA, "")
    sample = row.get("sample", row.get("DepMap_ID", ""))
    if not sample:
        raise ValueError("Missing sample identifier")
    out.update(
        source_table=spec.filename,
        sample_id=sample,
        genome_build=spec.genome_build,
        filter_provenance=spec.filter_provenance,
    )
    if table == "tcga_gtex":
        cohort = row["data_class"]
        if cohort not in {"TCGA", "GTEx"}:
            raise ValueError(f"Unknown data_class: {cohort}")
        if row["tumor_or_normal"] not in {"normal", "tumor"}:
            raise ValueError(f"Unknown tumor_or_normal: {row['tumor_or_normal']!r}")
        if cohort == "GTEx":
            donor = re.fullmatch(r"GTEX-([^-]+)(-.*)?", sample)
            if donor:
                tissue = _gtex_tissue(row["tissue_type"], row["tissue_subtype"])
                patient = f"GTEX-{donor[1]}"
                if tissue == "cell_line":
                    out["filter_provenance"] += "; tissue_class_override=GTEx_cultured_cells"
            elif re.fullmatch(r"K-562(?:-SM-[A-Za-z0-9]+)?", sample):
                _gtex_labels(row["tissue_type"], row["tissue_subtype"])
                # ATCC CCL-243: verified malignant line, never a normal donor.
                tissue, patient = "cell_line", ""
                out["filter_provenance"] += "; tissue_class_override=K-562_ATCC_CCL-243"
            else:
                raise ValueError(f"Unknown GTEx sample identifier: {sample}")
        else:
            match = re.fullmatch(r"([^-]+)-([^-]+-[^-]+)-(TP|NT|TM|TR|TB|TBM|TAP|TAM)", sample)
            if match is None:
                raise ValueError(f"Unknown Broad TCGA sample identifier: {sample}")
            if match[1] != row["tissue_type"]:
                raise ValueError(f"Conflicting TCGA disease metadata: {sample}")
            tissue = "normal" if match[3] == "NT" else "tumor"
            patient = f"TCGA-{match[2]}"
        if tissue != "cell_line" and row["tumor_or_normal"] != tissue:
            raise ValueError(f"Conflicting tissue metadata: {sample}")
        out.update(
            cohort=cohort,
            tissue_class=tissue,
            patient_id=patient,
            tumor_type=_available(row["tissue_type"]) if tissue == "tumor" else "",
        )
    elif table.startswith("target_"):
        info = metadata[sample] if table == "target_starfusion" else row
        tissue = _target_tissue(info["sample_type"])
        if table == "target_consolidated":
            if row["tumor_or_normal"] not in {"normal", "tumor"} or (
                tissue != "cell_line" and row["tumor_or_normal"] != tissue
            ):
                raise ValueError(f"Conflicting TARGET tissue metadata: {sample}")
        if tissue == "cell_line":
            out["filter_provenance"] += "; tissue_class_override=TARGET_cell_line_or_xenograft"
        out.update(
            cohort="TARGET",
            patient_id=_available(info["participant"]),
            tissue_class=tissue,
            tumor_type=_available(info["Cancer.Type"]) if tissue == "tumor" else "",
        )
    else:
        out.update(cohort="DepMap", tissue_class="cell_line")

    if table == "target_consolidated":
        # partnerA/B are ALPHABETICALLY sorted upstream, not oriented partners.
        out["gene_5p"], out["gene_3p"] = _fusion_genes(row, table)
        if sorted((out["gene_5p"], out["gene_3p"])) != sorted(
            (_available(row["partnerA"]), _available(row["partnerB"]))
        ):
            raise ValueError(
                f"Consolidated TARGET partner-identity mismatch: "
                f"sample={sample!r}, fusion={row['fusion_name']!r}"
            )
        validated = row["FI_validated"]
        if validated not in {"TRUE", "FALSE"}:
            raise ValueError(f"Invalid FI_validated: {validated}")
        suffix = ".FI" if validated == "TRUE" else ".starF"
        out["fi_validated"] = validated
        out["caller"] = "STAR-Fusion v1.10.0 + FusionInspector v2.4.0"
        out["filter_provenance"] += f"; normalized_evidence={suffix[1:]}"
    else:
        for side, key in (("5p", "LeftGene"), ("3p", "RightGene")):
            # ENSG identifiers are gene IDs, never transcript IDs.
            out[f"gene_{side}"] = _available(row[key]).split("^")[0].split(" (")[0]
        genes = _fusion_genes(row, table)
        if genes != (out["gene_5p"], out["gene_3p"]):
            raise ValueError(
                f"Fusion name disagrees with oriented genes: "
                f"sample={sample!r}, fusion={row.get('fusion_name', row.get('FusionName'))!r}"
            )
        suffix = ""
        out["caller"] = {
            "tcga_gtex": "STAR-Fusion v1.7",
            "target_starfusion": "STAR-Fusion v1.10.0",
            "depmap_22q2": "STAR-Fusion",
        }[table]
    for side, key in (("5p", "LeftBreakpoint"), ("3p", "RightBreakpoint")):
        if not _available(row[key]):
            continue
        match = re.fullmatch(r"([^:]+):([1-9][0-9]*)(?::([+-]))?", row[key])
        if not match:
            raise ValueError(f"Invalid breakpoint: {row[key]!r}")
        out[f"chrom_{side}"], out[f"pos_{side}"], out[f"strand_{side}"] = (
            value or "" for value in match.groups()
        )
    for field, key in (
        ("junction_reads", "JunctionReadCount"),
        ("spanning_reads", "SpanningFragCount"),
        ("ffpm", "FFPM"),
    ):
        if table == "tcga_gtex":
            key = {"JunctionReadCount": "Jcount", "SpanningFragCount": "Scount"}.get(key, key)
        out[field] = _number(row[key + suffix], field)
    out["splice_type"] = _available(row.get("SpliceType", ""))
    # Preserve original annotations and all ancillary flags/metrics verbatim, including
    # both callers' counts in the consolidated table. JSON avoids delimiter ambiguity.
    out["annotation_flags"] = json.dumps(row, ensure_ascii=False, sort_keys=True)
    return out


def read_table(
    table: str, path: Path, *, target_sample_info: Path | None = None
) -> Iterator[dict[str, str]]:
    """Yield every source row, without filtering or deduplication.

    Raw TARGET requires the pinned sample-info sidecar. Missing metadata fails.
    Use ``normalize_table`` for atomic file output (no partial TSV on any error).
    """
    spec = MANIFEST[table]
    if not spec.genome_build or not spec.build_evidence:
        raise ValueError(f"Unusable table without documented genome build: {table}")
    metadata = {}
    if table == "target_starfusion":
        if target_sample_info is None:
            raise ValueError("Raw TARGET requires TARGETdb.sample_info.tsv")
        metadata = _target_metadata(target_sample_info)
    for row in _rows(path, HEADERS[table], "," if table == "depmap_22q2" else "\t"):
        if table == "target_starfusion" and row["sample"] not in metadata:
            raise ValueError(f"Missing TARGET metadata: {row['sample']}")
        yield _normalize(row, table, spec, metadata)


def normalize_table(
    table: str, source: Path, destination: Path, *, target_sample_info: Path | None = None
) -> int:
    """Write exactly SCHEMA as TSV atomically; return the number of input rows."""
    count = 0
    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", newline="", dir=destination.parent, delete=False
        ) as handle:
            temporary = Path(handle.name)
            writer = csv.DictWriter(handle, fieldnames=SCHEMA, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row in read_table(table, source, target_sample_info=target_sample_info):
                writer.writerow(row)
                count += 1
        os.replace(temporary, destination)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    return count
