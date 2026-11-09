"""Strict, lossless-row readers for pinned open fusion candidate tables."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import re
import tempfile
from collections.abc import Iterator
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
    if "Normal" in sample_type:
        return "normal"
    if any(word in sample_type for word in ("Tumor", "Cancer", "Metastatic")):
        return "tumor"
    raise ValueError(f"Unknown TARGET sample type: {sample_type!r}")


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
        tissue = "normal" if cohort == "GTEx" or sample.endswith("-NT") else "tumor"
        if row["tumor_or_normal"] != tissue:
            raise ValueError(f"Conflicting tissue metadata: {sample}")
        # Broad uses disease-prefix TCGA IDs, e.g. LUSC-51-4079-NT.
        pattern = r"[^-]+-([^-]+-[^-]+)-(?:TP|NT|TM|TR|TB|TBM|TAP|TAM)$"
        match = re.fullmatch(pattern, sample) if cohort == "TCGA" else None
        if cohort == "TCGA" and match is None:
            raise ValueError(f"Unknown Broad TCGA sample identifier: {sample}")
        out.update(
            cohort=cohort,
            tissue_class=tissue,
            patient_id=f"TCGA-{match[1]}" if match else "-".join(sample.split("-")[:2]),
            tumor_type=_available(row["tissue_type"]) if tissue == "tumor" else "",
        )
    elif table.startswith("target_"):
        info = metadata[sample] if table == "target_starfusion" else row
        tissue = _target_tissue(info["sample_type"])
        if table == "target_consolidated" and row["tumor_or_normal"] != tissue:
            raise ValueError(f"Conflicting TARGET tissue metadata: {sample}")
        out.update(
            cohort="TARGET",
            patient_id=_available(info["participant"]),
            tissue_class=tissue,
            tumor_type=_available(info["Cancer.Type"]),
        )
    else:
        out.update(cohort="DepMap", tissue_class="cell_line")

    if table == "target_consolidated":
        # partnerA/B are ALPHABETICALLY sorted upstream, not oriented partners.
        genes = row["fusion_name"].split("::")
        if len(genes) != 2 or not all(genes):
            raise ValueError(f"Invalid oriented fusion_name: {row['fusion_name']}")
        out["gene_5p"], out["gene_3p"] = genes
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
        out[field] = _available(row[key + suffix])
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
            writer = csv.DictWriter(handle, fieldnames=SCHEMA, delimiter="\t")
            writer.writeheader()
            for row in read_table(table, source, target_sample_info=target_sample_info):
                writer.writerow(row)
                count += 1
        os.replace(temporary, destination)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    return count
