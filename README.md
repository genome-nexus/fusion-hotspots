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

The correction family includes retention/disruption Fisher and permutation
tests, corrected cutpoint/window scans, pair enrichment, confidence Welch
tests, both expression comparisons, and each mutation/CNA comparator test.
Composite scores are rankings and are not included as p-values.

Cohort scans fetch expression and comparator inputs when those algorithms
are requested and the gene/study configuration supplies the required metadata.
A failed or malformed comparator fetch skips that target; a successful fetch
with zero calls remains an observed zero. Direct algorithm callers should use
`comparator_alterations=None` for unavailable data and an explicit list
(including `[]`) for available data. For partial availability across targets,
pass `comparator_availability` keyed by `comparator_target_key(target)`.

The live pipeline supplies the same genomic null-model client and requested
permutation budget to retention and disruption, unless per-algorithm parameters
override them. New disruption p-values and expanded-family q-values can therefore
differ from historical saved reports; those reports are not updated automatically.

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

## Published viewer (GitHub Pages)

Every run under `runs/` can be rendered as a backend-less, static-site HTML
viewer (see `src/cfh/reporting/html_viewer.py`; also produced locally by
`cfh analyze`/`cfh cohort-scan --html`). Once GitHub Pages is enabled for
this repository (see below), `.github/workflows/pages-deploy.yml` publishes
a live, shareable copy of that viewer automatically: on every push to `main`
it regenerates the viewer bundle for each committed run directory — offline,
from the already-committed `results.json`/`summary.json`/`*.svg` artifacts,
with no network access — and deploys the result plus a landing page linking
to every run to GitHub Pages. It can also be run manually from the Actions
tab (`workflow_dispatch`).

**One-time manual step required (this PR cannot do this itself):** a repo
admin must enable Pages once, under **Settings > Pages > Source > GitHub
Actions**. Until that setting is flipped, the workflow's deploy step will
fail even though the workflow file itself is already merged.

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
  runs/braf_msk-impact-50k-2026_20260910T145723Z \
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

### Observation units and assay eligibility

Live ingestion preserves patient IDs supplied by the SV and clinical APIs. Run
summaries report event, sample, known-patient, and missing-patient counts. Domain
and spatial tests still use events: these counts do **not** establish patient
independence, and repeated biopsies require a patient-level design before their
p-values can support replication.

`frequency` with `dedup_by_patient=True` counts each patient–partner combination
once, retaining different partners in the same patient. Missing patient identity
falls back to sample, then event identity, with a warning. `count_unit` identifies
the denominator; the legacy `Event_count` column counts those units when deduplication
is enabled, and `Raw_event_count` preserves the input count.

Expression domain comparisons use one measurement per sample. Samples with both
retained and non-retained fusion events are excluded from that split and counted
in its diagnostics. Repeated biopsies remain separate samples and are warned
about when patient identity reveals them.

Live mutation/CNA co-occurrence intersects the study sample list with documented
coverage for the target SV gene and comparator gene/assay. Gene-panel membership
is fetched through the [cBioPortal API](https://www.cbioportal.org/api/v3/api-docs).
Unprofiled samples and unknown coverage are excluded; missing panel metadata is
conservatively treated as unknown, including studies that do not document their
genome-wide coverage. An unavailable eligibility source skips the comparison.
Each tested target records its eligible sample count, a SHA-256 digest of the
sorted eligible sample IDs, and exclusion counts. This follows the distinction between assayed and off-panel genes described in the
[cBioPortal profiling FAQ](https://docs.cbioportal.org/user-guide/faq/).

The live fusion-positive/negative expression comparison likewise requires SV gene
coverage and measured expression. If coverage is unavailable, the domain split
can still run on observed fusion samples, but no fusion-negative group is inferred.
Offline callers supplying only `cohort_sample_ids` are explicitly asserting that
those samples form an eligible universe. Co-occurrence calls outside that universe
are excluded; they never enlarge its denominator. These comparisons remain
unadjusted for tumor type and should be run in prespecified tissue strata when
that confounding matters.

### Interpretation and validation

Mechanistic summaries separate nominal Fisher evidence from the matching
FDR-adjusted result. Standalone runs have unknown cross-gene FDR; cohort scans
attach the corresponding retention/disruption Fisher q-value, rather than a
minimum q-value from another algorithm. Recurrent counter-pattern partners are
a descriptive review flag, not a statistical subcluster test. Curated mechanism
notes provide biological context, not functional validation of individual events.

The composite score is a prioritization heuristic, not a driver probability.
Cutpoint proximity requires a significant scan-corrected cutpoint p-value. Other
components can share gene-level evidence across partners, so independent ranking
validation remains necessary.

Run deterministic synthetic calibration with:

```bash
python -m cfh.stats.calibration --replicates 100 --n-patients 24 \
  --n-permutations 999 --seed 42 --family-size 100 > calibration.json
```

The benchmark reports rejection rates and Wilson intervals under an independent
null, a planted label-associated window, a mapping pile-up null, and a repeated-
patient stress case. It calibrates **label separation**, not excess breakpoint
density or driver classification. The permutation-resolution diagnostic reports
whether the selected budget can resolve small p-values relevant to a correction
family; it does not establish FDR control. Small smoke runs only verify execution.
Actual replication and composite-versus-recurrence performance require independent,
patient-disjoint cohorts with suitable coverage and separately reviewed labels.

See [CALIBRATION.md](CALIBRATION.md) for the held-out ranking input schema,
permutation-budget controls, and strict patient-disjoint cohort comparison.
