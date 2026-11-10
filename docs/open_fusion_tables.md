## Open fusion tables: candidate nomination only

These tables can **nominate candidates only; they cannot estimate prevalence**.
STAR-Fusion sequence-similarity filtering already ran, which can remove paralog
pairs. Disabling annotation filtering did not disable similarity or minimum
alignment-evidence filters. DepMap's “unfiltered” label refers to its later
postfilters. `open_fusion_manifest.py` pins source URLs, versions, SHA-256 hashes,
build evidence and upstream filtering for the four tables and the raw TARGET
sample metadata sidecar. All four coordinate builds are GRCh38: the
[FusionInspector paper's STAR Methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC10261907/)
names the Sept032019 CTAT library for TCGA+GTEx and Apr032020 for TARGET; the
[saved DepMap 22Q2 configuration](https://github.com/broadinstitute/depmap_omics/blob/80f54c2b34a6cc7a88b8898961cba2eee27d9dd8/data/22Q2/RNAconfig/inputs_hg38_STAR_fusion.json)
uses hg38 CTAT resources. DepMap's release README reports STAR-Fusion 1.6.0,
whereas its saved configuration reports 1.9.1; we retain `STAR-Fusion` without
asserting a per-row version. DepMap filtering provenance is explicitly config-derived.
TARGET filter evidence cites the v1.10.0 implementation as well as the paper.

```python
from pathlib import Path
from cfh.ingestion.open_fusion_manifest import MANIFEST
from cfh.ingestion.open_fusion_tables import fetch_table, normalize_table

root = Path("open_candidates")
metadata = fetch_table("target_sample_info", root / "target_metadata")
for table in MANIFEST:
    source = fetch_table(table, root / table)  # directory must not exist yet
    normalize_table(table, source, root / f"{table}.tsv", target_sample_info=metadata)
```

Downloads verify the compressed/source bytes before exposing the final filename.
A mismatch raises an error and removes the incomplete file. Readers require exact
headers, reject malformed rows, and write TSVs atomically so errors leave no
partial output. `read_table(table, path, target_sample_info=...)` also provides an
iterator; consumers of that iterator must discard earlier rows if iteration fails.

No rows are filtered or deduplicated, including normals, low counts, same-family
pairs, non-primary isoforms and flagged rows. `annotation_flags` is JSON containing
the original row, preserving annotation strings, ancillary flags and both sets of
TARGET counts verbatim. Missing source values (`NA`, `NaN`, `.`) become empty
normalized fields. For numeric measurements only empty, `NA` and `.` are missing;
NaN/Inf, malformed numbers, negative or fractional read counts, and negative or
nonfinite FFPM fail loudly. The real tables have only integral read counts.
Integral scientific notation such as `1e+05` is accepted for reads without
changing source text; fractional nonnegative FFPM remains valid.
TSVs use LF line endings. No transcript identifiers or junction/read sequences occur in
these four tables; all carry junction and spanning **counts**, not read sequences.
ENSG gene identifiers are not substituted for missing transcript identifiers.
Coordinates retain the source's genomic breakpoint positions and strands. Missing
coordinates or strands remain empty (including strandless TARGET breakpoints).

TCGA+GTEx follows the
[upstream decoration rule](https://github.com/broadinstitute/FusionInspectorPaper/blob/e8f4613f5caa8c2571c5ddf99483d7b1dde6647b/0.Decorate_input_files/Decorate_input_files.Rmd):
GTEx donor IDs must match `GTEX-[^-]+(-.*)?`. The 53 observed GTEx tissue-type/subtype pairs are explicitly pinned in the reader.
New noncell pairs require a mapping update. The allowlisted labels `Cells -
EBV-transformed lymphocytes`, `Cells - Transformed fibroblasts`, and `Cells -
Cultured fibroblasts` in either column are `cell_line`, retaining the real donor ID.
Case, hyphens, underscores and whitespace are normalized before exact lookup;
there is no fuzzy match. Any label beginning with `Cells` in either column that
is outside that allowlist raises, even if the other column has a known cell label.
Cultured fibroblasts are classified as `cell_line` because culture can introduce
acquired rearrangements. This policy keeps cultured material out of the normal
tissue panel. The local pinned copy has only
the EBV-transformed and transformed-fibroblast labels (4,066 and 3,158 rows).
The complete observed pair/count table is recorded in PR #151.
Known noncell GTEx pairs and TCGA IDs ending `-NT` are `normal`; recognized other
TCGA sample codes are `tumor`.
The sole verified non-donor-ID exception is `K-562` or `K-562-SM-<alphanumeric>`:
[ATCC CCL-243](https://www.atcc.org/products/ccl-243) identifies this as a malignant
leukemia line. It is `cell_line`, with no invented patient ID. These explicit
cell-line overrides retain the upstream tissue label in `annotation_flags` and
record the override in `filter_provenance`. Other non-GTEX IDs fail loudly.
The reader checks agreement with `tumor_or_normal`, except for the recorded
cell-line overrides. TCGA `tumor_type` uses `tissue_type` (e.g. LUSC), checked
against the sample disease prefix; `tumor_subtype` contains codes such as TP.
Broad's disease-prefix IDs (e.g. LUSC-51-4079-TP/NT) become patient TCGA-51-4079,
matching standard TCGA participant barcodes. GTEx participant IDs retain the first
two barcode components. TARGET uses the published participant, cancer type and
sample type; `Normal` sample types are normal and explicit tumor/cancer/metastatic
types are tumor using an explicit allowlist. Cell-line and xenograft labels
(case-insensitive, with separator normalization) are always `cell_line`, including
when the source classifies a xenograft as tumor; that override is recorded. None
occur in the real TARGET metadata's 11 sample types, all of which are covered.
Other unknown types fail. Across
all cohorts, normal rows have empty `tumor_type`; their original cancer-type
metadata remains in `annotation_flags`. DepMap
rows are `cell_line`; patient and tumor type are unavailable and remain empty.

The overlapping TARGET tables remain distinct. Consolidated TARGET `fusion_name`
supplies the 5′/3′ gene order; alphabetically sorted `partnerA/B` support only a
partner-identity check. They do not validate orientation or genomic locus alignment.
Validated rows use FI counts/FFPM, others STAR-Fusion counts/FFPM, recorded in
`filter_provenance`; `caller` identifies both methods and `fi_validated` preserves
the original boolean. Both count sets remain in `annotation_flags`. Patient IDs
support later deduplication across tables; this reader performs none. The real
TARGET metadata has no USI under multiple participant IDs, so patient IDs remain
published participant IDs.

Fusion names use `::` for TCGA+GTEx and consolidated TARGET, and `--` for raw
TARGET and DepMap. Readers validate this format and agreement with oriented
LeftGene/RightGene; consolidated partner fields are checked only for partner identity. Four TCGA+GTEx placeholder
rows have unavailable names, genes and breakpoints and remain present. Different
gene annotations at identical sample/breakpoints remain distinct records.
Breakpoint containment in a public GENCODE gene span is not an ingestion rule:
[STAR-Fusion v1.10](https://github.com/STAR-Fusion/STAR-Fusion/blob/v1.10.0/util/STAR-Fusion.map_chimeric_reads_to_genes)
permits breakpoints outside annotated exons;
[CTAT](https://github.com/TrinityCTAT/ctat-genome-lib-builder/blob/master/util/revise_gencode_annotations.pl)
also modifies reference annotations (including IG superloci). Ingestion preserves upstream locus assignments; independent gene-model
and breakpoint validation belongs downstream. It does not certify containment.
