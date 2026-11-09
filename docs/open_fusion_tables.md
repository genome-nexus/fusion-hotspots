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
asserting a per-row version.

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
normalized fields. No transcript identifiers or junction/read sequences occur in
these four tables; all carry junction and spanning **counts**, not read sequences.
ENSG gene identifiers are not substituted for missing transcript identifiers.
Coordinates retain the source's genomic breakpoint positions and strands. Missing
coordinates or strands remain empty (including strandless TARGET breakpoints).

TCGA+GTEx follows the
[upstream decoration rule](https://github.com/broadinstitute/FusionInspectorPaper/blob/e8f4613f5caa8c2571c5ddf99483d7b1dde6647b/0.Decorate_input_files/Decorate_input_files.Rmd):
GTEx and TCGA IDs ending `-NT` are `normal`; other TCGA rows are `tumor`.
The reader checks agreement with `tumor_or_normal`. TCGA `tumor_type` uses
`tissue_type` (e.g. LUSC); `tumor_subtype` contains sample codes such as TP.
Broad's disease-prefix IDs (e.g. LUSC-51-4079-TP/NT) become patient TCGA-51-4079,
matching standard TCGA participant barcodes. GTEx participant IDs retain the first
two barcode components. TARGET uses the published participant, cancer type and
sample type; `Normal` sample types are normal and explicit tumor/cancer/metastatic
types are tumor, matching the upstream annotation. Unknown types fail. DepMap
rows are `cell_line`; patient and tumor type are unavailable and remain empty.

The overlapping TARGET tables remain distinct. Consolidated TARGET `fusion_name`
is oriented 5′ to 3′; alphabetically sorted `partnerA/B` are not used for orientation.
Validated rows use FI counts/FFPM, others STAR-Fusion counts/FFPM, recorded in
`filter_provenance`; `caller` identifies both methods and `fi_validated` preserves
the original boolean. Both count sets remain in `annotation_flags`. Patient IDs
support later deduplication across tables; this reader performs none.
