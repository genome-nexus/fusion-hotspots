# fusion-hotspots (cfh)

[![CI](https://github.com/genome-nexus/fusion-hotspots/actions/workflows/ci.yml/badge.svg)](https://github.com/genome-nexus/fusion-hotspots/actions/workflows/ci.yml)

`cfh` detects gene-fusion breakpoint and domain-retention hotspots from
structural-variant calls in cancer sequencing cohorts. Milestone 1 targets
**BRAF** fusions in the public MSK-IMPACT 50k cohort on cBioPortal, using a
design that generalizes to other genes (already proven on RET plus a TCGA
PanCancer Atlas holdout) without code changes.

## What this is

Gene fusions can retain, disrupt, or lose functional protein domains
depending on exactly where the breakpoint falls. `cfh` builds a reproducible
pipeline that:

1. **Ingests** structural-variant (SV) and clinical data from cBioPortal,
   either from a downloaded study archive or the cBioPortal REST API.
2. **Normalizes** raw SV rows into a typed `FusionEvent` model, classifying
   event type (fusion / deletion / inversion / translocation), reading
   frame, 5'/3' orientation, and tumor-type/OncoTree provenance (joined
   from clinical data) without ever guessing when the source data is
   ambiguous.
3. **Maps** each event's breakpoint onto transcript exons and protein
   domains (via RefSeq annotations, Genome Nexus's canonical-transcript
   and exon-coordinate data as the primary generic fallback, and an
   Ensembl REST / UniProt+InterPro cross-check), producing a
   `FusionFeature` describing which domains are retained, disrupted, or
   lost, plus each domain's retained amino-acid interval, retained fraction,
   and truncation state when its boundaries are known.
4. **Runs pluggable hotspot-detection algorithms** against the normalized
   events and features, each producing a structured `AlgorithmResult`.

Gene-specific biology (canonical transcript, protein accession, domain
boundaries) lives in a per-gene YAML config (see
`src/cfh/genes/configs/`), never hardcoded into the generic
ingestion/normalization/mapping code, so adding a new gene is a config
change, not a code change.

## Install

Requires Python 3.10+.

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

## Quickstart

```python
from cfh.genes.registry import load_gene_config
from cfh.ingestion.archive_reader import load_sv_dataframe
from cfh.ingestion.clinical_parser import parse_clinical_sample
from cfh.normalization.event_normalizer import normalize

gene_config = load_gene_config("braf")

raw_sv = load_sv_dataframe("path/to/msk_impact_50k_study")
clinical = parse_clinical_sample("path/to/msk_impact_50k_study/data_clinical_sample.txt")

events = normalize(raw_sv, clinical, cohort="msk_impact_50k_2026")
print(f"Normalized {len(events)} fusion events for {gene_config.gene_symbol}")
```

Run the test suite (excluding tests that hit real external services):

```bash
pytest -m "not network"
```

Run the live cBioPortal/Genome Nexus benchmark. Each invocation creates a
timestamped directory under `runs/` with a manifest, event-level TSV, structured
JSON, Markdown report, discrepancy table, and SVG visualizations. Before
committing a new run, prune any older run of the same type (see
[CONTRIBUTING.md](CONTRIBUTING.md#committing-run-artifacts)):

```bash
cfh real-benchmark BRAF msk_impact_50k_2026
```

The study ID is an argument to the shared pipeline. For example, the original
MSK-IMPACT publication cohort can be run with:

```bash
cfh real-benchmark BRAF msk_impact_2017
```

To run the complete registered-algorithm orchestrator for a configured gene and
study, use:

```bash
cfh analyze BRAF msk_impact_50k_2026
```

To apply Benjamini-Hochberg FDR correction to the inferential p-values in two
or more existing run artifacts, pass their run directories (or `results.json`
files) to the offline comparison command:

```bash
cfh compare-genes RUN_DIR [RUN_DIR ...] --output adjusted-p-values.tsv
```

The output records the gene, study, algorithm, test, raw p-value, BH-adjusted
q-value, and whether the result is significant at `q < 0.05`. All p-values
collected by one invocation form a single correction family.

RET uses the same command and live ingestion/mapping path:

```bash
cfh analyze RET msk_impact_50k_2026
```

The TCGA PanCancer Atlas studies expose structural variants under the same
`<study_id>_structural_variants` profile convention. The BRAF holdout uses the
thyroid carcinoma study, which contains BRAF fusion calls:

```bash
cfh real-benchmark BRAF thca_tcga_pan_can_atlas_2018
```

All 32 PanCancer Atlas study IDs with `data_sv.txt`/`meta_sv.txt` in the public
cBioPortal Datahub are: `acc`, `blca`, `brca`, `cesc`, `chol`, `coadread`,
`dlbc`, `esca`, `gbm`, `hnsc`, `kich`, `kirc`, `kirp`, `laml`, `lgg`, `lihc`,
`luad`, `lusc`, `meso`, `ov`, `paad`, `pcpg`, `prad`, `sarc`, `skcm`, `stad`,
`tgct`, `thca`, `thym`, `ucec`, `ucs`, and `uvm`, each followed by
`_tcga_pan_can_atlas_2018`.

These commands make unauthenticated requests to both public services.

## Repository layout

```
src/cfh/model/         FusionEvent, FusionFeature, AlgorithmResult schemas
src/cfh/genes/         Per-gene configuration registry (e.g. braf.yaml)
src/cfh/algorithms/    Hotspot-detection algorithm plugin interface + registry
src/cfh/ingestion/     cBioPortal archive/API ingestion
src/cfh/normalization/ Raw SV rows -> FusionEvent
src/cfh/mapping/       Transcript/exon/domain mapping -> FusionFeature
```

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md).

## License

Apache License 2.0 — see [LICENSE](LICENSE).

### Offline cross-cohort CMH comparison

`cfh compare-cohorts RUN_DIR [RUN_DIR ...] --output cmh.json` reads the saved
`summary.frame_domain_contingency_table` for one gene across distinct cohorts.
It performs the uncorrected CMH chi-square test (one degree of freedom) and
reports the Mantel-Haenszel common odds ratio, the original tables, informative
strata, and sample-ID overlap. No ingestion or live algorithm execution occurs.
The optional output file contains the same JSON printed to stdout. Counts must
be nonnegative integers; empty/fixed-margin strata contribute no information.
An entirely uninformative input returns null statistic/p-value; an undefined
odds ratio is null and an infinite odds ratio is the string `"infinity"`.

Reproduce the BRAF comparison using exactly the three committed gene runs:

```bash
cfh compare-cohorts \
  runs/braf_msk-impact-50k-2026_20260909T181926Z \
  runs/braf_msk-impact-2017_20260905T012645Z \
  runs/braf_thca-tcga-pan-can-atlas-2018_20260909T034447Z \
  --output /tmp/braf-cmh.json
```

Rows are domain retained/not retained, and columns are in-frame protein fusion/other:

| Cohort | Saved table | Informative for CMH |
| --- | --- | --- |
| MSK IMPACT 50k 2026 | `[[142, 21], [9, 6]]` | Yes |
| MSK IMPACT 2017 | `[[31, 2], [2, 5]]` | Yes |
| TCGA thyroid 2018 | `[[9, 0], [6, 0]]` | No |

The nominal result is **CMH χ²(1) = 21.36486008, p = 3.796664805 × 10⁻⁶,
common OR = 7.455270793** (alpha = 0.05; no continuity correction). Both
informative cohorts have a positive association. **This does not confirm
concordance across all three cohorts:** thyroid has no comparison-column
observations, and the MSK artifacts share 34 sample IDs, violating independence
between strata. The p-value is therefore a nominal calculation, not valid
evidence of independent replication. The saved tables are intentionally not
deduplicated or recomputed.

CMH tests conditional association; it does not test equality of cohort odds
ratios or demonstrate significance within every cohort. It assumes independent
observations within and between strata, consistent table orientation, and
compatible domain definitions/counting units. The saved kinase boundaries are
458–712 aa in MSK and 457–712 aa in thyroid; metadata is retained in the report.
Sample IDs cannot rule out patient overlap or repeated observations. See the
[StratifiedTable methodological notes](https://www.statsmodels.org/v0.12.2/generated/statsmodels.stats.contingency_tables.StratifiedTable.html)
for independence assumptions. The implementation uses the standard formula and
SciPy's chi-square survival function, with no additional dependency; its
numerical regression uses the published
[Berkeley admissions CMH example](https://www.markirwin.net/stat149/Lecture/Lecture8.pdf).
