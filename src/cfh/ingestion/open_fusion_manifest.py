"""Pinned open candidate tables; hashes cover downloaded bytes (including gzip).

These are candidate nomination resources, never an unfiltered fusion census.
"""

from dataclasses import dataclass

BROAD_COMMIT = "e8f4613f5caa8c2571c5ddf99483d7b1dde6647b"
BROAD = f"https://raw.githubusercontent.com/broadinstitute/FusionInspectorPaper/{BROAD_COMMIT}"
PAPER = "https://pmc.ncbi.nlm.nih.gov/articles/PMC10261907/"
DEPMAP_CODE = (
    "https://github.com/broadinstitute/depmap_omics/blob/"
    "80f54c2b34a6cc7a88b8898961cba2eee27d9dd8/data/22Q2/RNAconfig/"
)


@dataclass(frozen=True)
class SourceTable:
    filename: str
    url: str
    version: str
    sha256: str
    genome_build: str
    build_evidence: str
    filter_provenance: str
    filter_evidence: str


_BROAD_FILTER = (
    "similarity_filter_ran; annotation_filter_off; min_FFPM=0; caller_evidence_filters_ran"
)
_BROAD_EVIDENCE = (
    f"{PAPER} STAR Methods: STAR-Fusion prediction of fusions in TCGA and GTEx; "
    "and Application of FusionInspector to cancer and normal transcriptomes. "
    "Commands disable annotation filtering, not sequence similarity filtering. "
    "https://github.com/STAR-Fusion/STAR-Fusion/blob/v1.7.0/util/STAR-Fusion.filter "
    "runs FusionFilter/blast_and_promiscuity_filter.pl separately from annotation filtering."
    " TARGET v1.10 implementation: "
    "https://github.com/STAR-Fusion/STAR-Fusion/blob/v1.10.0/util/STAR-Fusion.filter"
)

MANIFEST = {
    "tcga_gtex": SourceTable(
        "TCGA_n_GTEx.STAR-Fusion.v1.7.fusion_preds.tsv.gz",
        f"{BROAD}/data/TCGA_n_GTEx.STAR-Fusion.v1.7.fusion_preds.tsv.gz",
        "STAR-Fusion v1.7; TCGA v11 / GTEx v8",
        "db6e61750b33def456dde0d9430ecd923985bf931295f0039b2b446f59dfbcf4",
        "GRCh38",
        f"{PAPER} STAR Methods, STAR-Fusion prediction of fusions in TCGA and GTEx: "
        "CTAT GRCh38_gencode_v22_CTAT_lib_Sept032019.",
        _BROAD_FILTER,
        _BROAD_EVIDENCE,
    ),
    "target_starfusion": SourceTable(
        "STAR-Fusion.v1.10.0.TARGETdb.tsv.gz",
        f"{BROAD}/data/STAR-Fusion.v1.10.0.TARGETdb.tsv.gz",
        "STAR-Fusion v1.10.0",
        "6814bec6dc4209eec68876b8ac887d20e73d308a191c8c12b9c127838edd1709",
        "GRCh38",
        f"{PAPER} STAR Methods, Application of FusionInspector to cancer and normal "
        "transcriptomes: GRCh38_gencode_v22_CTAT_lib_Apr032020; TARGET processed "
        "as above with STAR-Fusion v1.10.0 then FusionInspector v2.4.0.",
        _BROAD_FILTER,
        _BROAD_EVIDENCE,
    ),
    "target_consolidated": SourceTable(
        "TARGET_STAR-Fusion-v1.10.0_FusionInspector-v2.4.0.fusions.tsv.gz",
        f"{BROAD}/data/TARGET_STAR-Fusion-v1.10.0_FusionInspector-v2.4.0.fusions.tsv.gz",
        "STAR-Fusion v1.10.0 + FusionInspector v2.4.0",
        "9a91c974ca28881ae765294c8b405c39bb78ea2cd87096e6b1d50b7718a532f3",
        "GRCh38",
        f"{PAPER} STAR Methods, Application of FusionInspector to cancer and normal "
        "transcriptomes: TARGET processed with GRCh38_gencode_v22_CTAT_lib_Apr032020. "
        f"{BROAD}/7.TARGET_cancer_db_FI_application/7.5.Consolidate_STARF_and_FusionInspector/"
        "Consolidate_STARF_and_FusionInspector.Rmd joins STAR-Fusion and FI breakpoints.",
        _BROAD_FILTER + "; FI_validation_evidence_filters_ran; STARF_FI_full_join",
        _BROAD_EVIDENCE,
    ),
    "depmap_22q2": SourceTable(
        "CCLE_fusions_unfiltered.csv",
        "https://ndownloader.figshare.com/files/34989934",
        "DepMap 22Q2 Public v2; STAR-Fusion (release README: v1.6.0)",
        "5e69e36e8092097b3ce549d8608249419468164afbbb8823b86e579f152e5d06",
        "GRCh38",
        f"{DEPMAP_CODE}inputs_hg38_STAR_fusion.json specifies "
        "workspace.hg38_star_fusion_ctat_files_v33 and hg38 STAR index. "
        "Release https://api.figshare.com/v2/articles/19700056/versions/2 links "
        "README https://ndownloader.figshare.com/files/34990054 (Fusions: GENCODE v29). "
        "README reports v1.6.0, saved 22Q2 config reports v1.9.1/v33: build agrees, "
        "per-row caller version cannot be established.",
        _BROAD_FILTER + "; depmap_postfilters_off; red_herring_flags_retained"
        "; provenance_basis=22Q2_saved_config",
        "https://ndownloader.figshare.com/files/34990054 Fusions and "
        f"CCLE_fusions_unfiltered.csv sections; {DEPMAP_CODE}hg38_STAR_fusion.wdl "
        "passes --no_annotation_filter --min_FFPM 0 without disabling similarity filtering. "
        "Filtering provenance is config-derived; the conflicting caller versions mean "
        "the exact historical per-row invocation is not established.",
    ),
}

# Required sidecar for raw TARGET: do not guess tissue or cancer type from project codes.
TARGET_SAMPLE_INFO = SourceTable(
    "TARGETdb.sample_info.tsv",
    f"{BROAD}/data/TARGETdb.sample_info.tsv",
    f"FusionInspectorPaper {BROAD_COMMIT}",
    "0fe08262d39a8e69cb6fa79a7dacb5f97a3173aa0a8d8485fb340248ff8462aa",
    "",
    "Sample metadata only; no genomic coordinates.",
    "",
    "",
)
