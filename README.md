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
collected by one invocation form a single correction family. See
[why Benjamini-Hochberg](#why-benjamini-hochberg) for the choice of correction.

The correction family includes retention/disruption Fisher and permutation
tests, pair enrichment, both expression comparisons, and each mutation/CNA
comparator test. Composite scores are
rankings and are not included as p-values. The `confidence_stats` Welch
read-support test is reported per gene but excluded from the family: it
reflects call quality, not breakpoint selection (#114). The cutpoint and window
scans are excluded too: they permute retained/lost labels that breakpoint
position already determines, so they find a clean split at the domain edge even
in uniformly random breakpoints. Their boundaries remain descriptive output. Cohort
scans correct the junction-concentration test as a family of its own
(`junction_concentration_family.tsv`, `junction_hotspot_candidates` in
`summary.json`); it never makes a gene FDR-significant, because on capture-panel
data its null favors significance (see the panel caveat below).

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

## Key features

Each registered algorithm writes one `AlgorithmResult` (`Summary`, `Tables`,
`Warnings`) per gene. Optional algorithms are no-ops unless the gene config or
study supplies their inputs.

| Algorithm | Question it answers | Output to read | Runs when |
| --- | --- | --- | --- |
| `frequency` | Which partners occur, and how often? | `Partner_gene_counts` (events, samples, patients) | Always |
| `domain_retention` | Are in-frame fusions enriched for retaining the key domain? | Fisher p/OR, permutation p, 2×2 table | Always |
| `domain_disruption` | Are in-frame fusions enriched for *losing* a configured domain? | Same as retention | `disruption_required_domains` set |
| `exon_retention` | Is a target exon retained? | Retention rate per exon | Exon hint configured |
| `confidence_stats` | In-frame rate (Wilson CI) and read support, retained vs not | MLE/CI per group, Welch t-test | Always (gene-agnostic defaults) |
| `junction_concentration` | Does one exon hold more of the gene's fusion junctions than random breakage predicts? | Top junction (observed vs random-breakage share), permutation p; see the panel caveat | Always |
| `cutpoint_detection` | Where is the protein position that best separates retained from lost? | Inferred cutpoint (aa); descriptive, its p-value is not a selection test | Always |
| `window_detection` | Is there an internal protein window that separates them? | Best window; descriptive, like the cutpoint | Always |
| `genomic_position_recurrence` | Do breakpoints pile up at the same DNA coordinate, or only after exon clamping? | Genomic clusters vs protein positions | Always |
| `joint_partner` | Is an ordered pair (e.g. EML4→ALK) enriched over a marginal-independence null? | Pair Fisher test | `gene_pair` configs |
| `mutation_cooccurrence` | Are fusions co-occurring with, or exclusive of, a comparator mutation/CNA? | Fisher p, direction, eligible-sample counts | `mutual_exclusivity_targets` set |
| `expression_association` | Does target-gene mRNA differ by fusion status or domain status? | Welch/Mann-Whitney per comparison | Study has an mRNA profile |
| `mechanistic_interpretation` | Which events run counter to a nominally supported retention/disruption pattern? | Counter-pattern events, recurrent-partner flag, Fisher q in scans | Always |
| `composite_score` | Which partners should be reviewed first? | Ranked partner table (heuristic, not a probability) | Always |

Cross-run tools:

- `cfh cohort-scan`: every recurrently altered gene, with auto-configured genes
  and BH-FDR across all genes' p-values.
- `cfh compare-genes`: BH-FDR across saved runs.
- `cfh compare-cohorts`: CMH across cohorts; `--require-patient-disjoint`
  enforces independence.
- `python -m cfh.stats.calibration`: synthetic null/power stress tests.
- `python -m cfh.stats.independent_validation`: held-out ranking comparison.

## Analyzing a cohort scan

```bash
cfh cohort-scan msk_impact_50k_2026 --output-dir runs --no-pdf   # --quiet hides per-gene progress
```

The run writes `runs/cohort-scan_<study>_<timestamp>/`:

| File | Use it for |
| --- | --- |
| `cohort_scan/summary.tsv` / `summary.json` | One row per gated gene: status, counts, frame and domain-retention percentages, retention/disruption Fisher p and q, minimum q, top composite partner |
| `cohort_scan/fdr_family.tsv` | Every hypothesis in the BH family: raw p, q, permutations used, and `resolution_limited` when a permutation p-value sits at its floor `1/(n+1)` without reaching significance |
| `cohort_scan/summary.md` | FDR-significant genes, honorable mentions (top raw p among non-significant genes), warnings |
| `cohort_scan/manhattan.svg` | Genome-wide view of per-gene significance |
| `cohort_scan/gene_reports/<gene>/` | Full per-gene run for significant, honorable-mention, and curated genes: `results.json`, event-level `results.tsv`, `report.md`, `outliers.tsv`, figures |
| `viewer/index.html` | Static HTML browser for the above (on by default; `--no-html` skips it) |

Permutation tests start at `--n-permutations-small` (100) and move to
`--n-permutations` (1,000) unless clearly non-significant (p > 0.1). A result still
at its floor is rerun with ten times the permutations up to `--n-permutations-max`
(100,000), because a p-value of 1/101 cannot pass BH across a family of thousands
of tests however strong the signal. A `resolution_limited` row in `fdr_family.tsv`
means the budget, not the data, decided that call.

Suggested reading order:

0. **Check the positive controls first.** `cfh eval-controls runs/<scan>` scores
   the scan against textbook fusion biology listed with sources in
   `src/cfh/benchmarks/positive_controls.yaml`: kinase or ETS-domain retention for
   ALK, RET, ROS1, BRAF, FGFR2, NTRK1, NTRK3 and ERG at genome-wide FDR, junction
   concentration for the same genes, their literature junction exons and partners, and four co-occurrence pairs (ERG–SPOP,
   ALK–EGFR, ALK–KRAS, ERG–PTEN). It writes `positive_controls.tsv`. If these fail,
   the scan's silence on other genes says little.
1. **Check the gate and failures.** Compare `genes_after_gating` with `status == "ok"`.
   Genes the Genome Nexus batch call cannot configure are retried with the
   single-gene endpoint, accepted only when its Ensembl gene ID matches HGNC's
   (`transcript_source` records which path was used). Non-coding genes
   (`config_source == "non_coding"`) and still-unresolved genes are listed
   separately in the warnings.
2. **Identify which test made a gene significant.** `fdr_significant` and
   `min_fdr_adjusted_q_value` use the smallest q across *every* test for that gene.
   For a domain claim, read `retention_fisher_q_value` / `disruption_fisher_q_value`
   instead. To list the driving tests for a gene with a saved report:

   ```python
   from pathlib import Path
   from cfh.gene_comparison import collect_p_values

   run = Path("runs/<scan>/cohort_scan/gene_reports/ntrk3")
   for row in sorted(collect_p_values([run]), key=lambda r: r["raw_p"])[:5]:
       print(row["algorithm"], row["test"], row["raw_p"])
   ```

   The `confidence_stats` Welch test (tumor read support, retained vs
   non-retained) is not in this family and never makes a gene significant. It
   reflects call quality, not breakpoint selection, and can rest on very small
   groups (check `n_b` in the gene's results).
3. **Check the counting units.** Each gene report's `summary` lists `total_fusions`,
   `sample_count`, and `known_patient_count`. Before any inferential test, each
   patient counts once per gene (`inference_counting_unit: patient`): repeated
   observations of one fusion (same partner, role, and protein junction) collapse to
   one, and of a patient's distinct fusions of the gene one representative is kept
   (determinate frame first, then earliest sample). `collapsed_repeat_observation_count`
   and `merged_distinct_fusion_count` show how many of each were dropped. Descriptive
   counts, the events table, and `frequency` keep every event. A patient with fusions
   in several genes still contributes to each gene's tests (see
   [CALIBRATION.md](CALIBRATION.md)).
4. **Read the domain evidence.** `domain_retention` gives the frame × domain 2×2
   table, Fisher and permutation p. The permutation null samples breakpoints across
   the transcript and classifies each from the domain's coordinates and the observed
   5'/3' roles (`permutation_null_classifier: domain_coordinates`); features lacking
   coordinates fall back to the nearest observed label. `cutpoint_detection` / `window_detection` locate
   the boundary and compare it with Pfam coordinates. Their labels are a function of
   breakpoint position and role, so a clean split at the domain boundary exists by
   construction; their p-values say there are enough events on both sides, not
   that breakpoints are selected, and they are not in the FDR family. The
   permutation retention test is the selection test. `genomic_position_recurrence`
   separates true DNA hotspots from exon-boundary clamping.
5. **Read the junction evidence, with the panel caveat.** `junction_concentration`
   asks whether one exon holds more of the gene's junctions than breakpoints drawn
   uniformly along the gene would put there, keeping each fusion's 5'/3' role and
   using the same exon mapping as the observed calls. `top_junction` gives the
   observed share and the random-breakage share; breakpoints outside the canonical
   transcript are left out (`n_outside_transcript`). Cohort scans BH-correct it across
   genes as a separate family and list the genes below q < 0.05 as
   `junction_hotspot_candidates`; it never sets `fdr_significant`. See
   [the panel caveat](#panel-capture-and-the-junction-concentration-test) before
   reading it on panel data.
6. **Review exceptions.** In `mechanistic_interpretation`, look at counter-pattern
   events and recurrent-partner flags. The HTML viewer shows curated mechanism notes
   beside them.
7. **Prioritize partners.** Read `composite_score` together with `frequency`.
   Cutpoint proximity contributes only when the corrected cutpoint p < 0.05, and
   gene-level components are shared by every partner of a gene.
8. **Consider confounding.** The Fisher and permutation tests pool tumor types.
   `domain_retention` and `domain_disruption` also report `tumor_type_stratified`: a
   Cochran-Mantel-Haenszel test across OncoTree strata with the Mantel-Haenszel common
   odds ratio, and the per-stratum tables in `frame_domain_contingency_tables_by_tumor_type`.
   The cohort summary carries `retention_cmh_p_value` / `retention_cmh_q_value` (BH across
   genes as a separate family, not used for `fdr_significant`). A pooled association that
   disappears after stratification is likely tumor-type composition. CMH is two-sided,
   while the pooled Fisher test is one-sided.
9. **Check fragility.** The in-frame vs other comparison rests on the table's small
   cells, usually the few patients whose fusion loses the domain, so a small p-value can
   hang on a handful of them. `fisher_fragility` and
   `tumor_type_stratified.fragility` give the fewest patients whose removal from the
   cells supporting the association raises p to 0.05 or above (greedy search, so an
   upper bound on the true minimum), with the strata and cells removed. The cohort
   summary carries `retention_fisher_fragility_index`, `retention_cmh_fragility_index`
   and `disruption_fisher_fragility_index`. An index of 1-3 means the result stands or
   falls with those patients; the breakpoint-redraw permutation test uses every patient
   and is the sturdier evidence.

### Co-occurrence discovery

`cfh cohort-scan <study> --discover-cooccurrence` additionally tests every analyzed
fusion gene against every panel gene × alteration type (mutation, amplification,
deep deletion), with no curated targets needed. It writes
`cohort_scan/cooccurrence_discovery.tsv` / `.json` and a summary section.

- **Assay eligibility:** each pair uses only samples whose SV panel covers the fusion
  gene and whose mutation or CNA panel covers the comparator gene.
- **Tumor type:** the test is Cochran-Mantel-Haenszel stratified by OncoTree code;
  `mh_common_odds_ratio` < 1 means mutual exclusivity. Pooled odds ratios are shown
  only for comparison.
- **Sparse pairs:** pairs failing the Mantel-Fleiss criterion (expected co-occurrence
  too close to its bounds for the chi-square approximation) get no p-value.
- **Multiplicity:** q-values are BH-adjusted across all tested pairs as their own
  family; they never affect `fdr_significant` in the main summary.
- **Hypermutation:** `--hypermutated-min-mutations N` excludes samples with at least
  `N` panel mutations.

The TSV keeps pairs with q < 0.25 plus each fusion gene's five smallest-p pairs.
Copy-number co-occurrence with a gene near the fusion locus (e.g. CDK12 fusions with
ERBB2 amplification on 17q12) usually reflects a rearrangement inside the same
amplicon, not two independent events. Samples are the counting unit.

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

### Panel capture and the junction-concentration test

The junction-concentration null assumes a breakpoint anywhere in the gene would be
detected. MSK-IMPACT and other hybrid-capture panels sequence only selected
introns of each fusion gene, chosen because known fusions break there. A
breakpoint in an intron the panel does not capture is never called, so the
observed junctions are concentrated at the captured introns partly by assay
design. What that does to the results:

- **The p-value overstates the evidence.** On panel data a significant
  `junction_concentration` result means the junctions are more concentrated than
  whole-gene random breakage, which the assay itself guarantees to some degree.
  It does not show, on its own, that the concentration is biological.
- **The design is circular.** Introns were chosen for capture because fusions
  were known to break there, so the panel finds them there again.
- **Genes with uncaptured hotspots are missed.** Fusions at an intron the panel
  does not tile are absent, which can hide a real hotspot or, for genes the panel
  barely tiles, leave too few events to test.
- **The domain-retention test inherits a milder form.** Its null also draws
  breakpoints across the whole gene, so which retention outcomes are detectable
  is shaped by capture too.

In the 2026-10-22 MSK-IMPACT scan, every textbook control gene was concentrated,
but so were tumor suppressors such as APC, ARID1A, KMT2C, NOTCH1, SPEN and KDM6A,
some with under 20 events. APC's top junction is its large last exon (90% observed
against 8% under random breakage). That pattern fits exon-biased detection rather
than fusion biology, and is why the test is its own family.

The repository does not include the panel's bait coordinates. With them, the null
could redraw breakpoints only within captured regions. Until then, read
`junction_concentration` on panel data as "consistent with" a hotspot and confirm
it in a cohort sequenced without intron selection (whole genome or RNA, such as
the TCGA replication runs). `panel_capture_caveat` in each result repeats this.

### Why Benjamini-Hochberg

A cohort scan is a screen: it proposes candidate genes for follow-up, so the error
rate that matters is the share of reported genes that are wrong (FDR), not the
chance of any single false positive (family-wise error). On the 2026-10-26
MSK-IMPACT run (874 tests), the alternatives compare as follows:

| Method | Controls | FDR-significant genes |
|---|---|---|
| Benjamini-Hochberg (used) | FDR, under independence or positive dependence | 17 |
| Storey q-value | FDR, estimating the null share π0 | 17 (π0 estimate is 1.0) |
| Benjamini-Yekutieli | FDR, under any dependence | 15 |
| Holm / Bonferroni | family-wise error | 15 |

- **Storey gains nothing here.** It gains power by estimating that some tests are
  non-null, but 77% of p-values exceed 0.5 and 46% equal 1, so π0 is estimated
  at 1 and the q-values equal BH's. That is typical of a genome-wide screen in
  which most genes are not drivers.
- **BH's dependence assumption is the one to defend.** Tests are not independent:
  Fisher and permutation tests on a gene share its patients, and a fusion such as
  EML4-ALK is counted under both genes. Both kinds of dependence push p-values in
  the same direction, the positive regression dependence under which BH is proven
  to control FDR. Benjamini-Yekutieli drops that assumption at the cost of
  dividing α by about 7.35 (the harmonic sum over 874 tests).
- **The choice only matters at the margin.** Most significant tests sit at the
  permutation floor and survive every method. ALK and NAB2 do not: ALK's
  retention q is 0.017 under BH but it fails under Benjamini-Yekutieli and Holm,
  so ALK's genome-wide significance rests on BH's dependence assumption.
- **BH is conservative for discrete tests.** Fisher tests on tables with few
  patients can reach only a handful of p-values, many of them 1, and BH counts
  them as full tests. Discrete-aware procedures (Tarone's screen, Gilbert 2005,
  Heyse 2011) drop tests that could never reach significance and would give
  borderline genes more room; they are not implemented.

### Observation units and assay eligibility

Live ingestion preserves patient IDs supplied by the SV and clinical APIs. Run
summaries report event, sample, known-patient, and missing-patient counts. Domain,
spatial, and junction tests count each patient once per gene (see the reading
order above); `patient_independence_established` is true when every event had a
patient ID. Patients with fusions in several genes still contribute to each gene's
test, so the cross-gene family is not fully patient-independent.

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
