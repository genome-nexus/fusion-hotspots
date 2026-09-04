# Genome-wide fusion-hotspot analysis of msk_impact_50k_2026

## Abstract

This manuscript synthesizes a genome-wide fusion-hotspot cohort scan of msk_impact_50k_2026: 3919 genes carried at least one structural-variant record in the cohort, of which 544 passed the >= 5-distinct-patient recurrence gate. All 544 gated genes were attempted: 4 using hand-curated gene configs and 520 auto-configured, with 20 gated-in gene(s) left unresolvable; 523 of the 544 attempted genes were successfully analyzed with the full registered algorithm suite. 4 genes reached genome-wide Benjamini-Hochberg FDR significance (q < 0.05): NTRK3 (58 events, 94.8% in-frame, 94.8% domain-retained, Fisher p=0.0695489, q=8.1653e-09); ROS1 (122 events, 59.0% in-frame, 87.7% domain-retained, Fisher p=0.214047, q=0.000221281); ETV6 (90 events, 71.1% in-frame, 75.6% domain-retained, Fisher p=5.12326e-06, q=0.00233279); FLI1 (118 events, 61.9% in-frame, 5.9% domain-retained, Fisher p=0.988837, q=0.00820622). 15 additional genes form a highly ranked non-FDR-significant tier flagged for targeted follow-up (see Honorable mentions, below).

## Methods

Structural-variant records were retrieved from the msk_impact_50k_2026 cBioPortal study and gated to genes with at least 5 distinct patient(s) carrying a structural-variant record (544 of 3919 genes passed the gate). All 544 gated genes were attempted: 4 using hand-curated gene configs and 520 auto-configured from Genome Nexus canonical-transcript/Pfam-domain data, with 20 gated-in gene(s) left unresolvable; 523 of the 544 attempted genes were successfully analyzed. Each successfully analyzed gene was run through the algorithm suite recorded for this scan (composite_score, confidence_stats, cutpoint_detection, domain_disruption, domain_retention, exon_retention, frequency, joint_partner). Domain-retention and domain-disruption significance were assessed per gene with Fisher's exact test and a breakpoint-position permutation test; the resulting p-values across the 359 genes that produced at least one computable p-value were jointly corrected with Benjamini-Hochberg false-discovery-rate correction at q < 0.05.

## Results

### Genome-wide summary

Genome-wide summary of 359 scanned genes with an FDR-adjusted q-value, ranked left-to-right by composite evidence score; the dashed line marks the q=0.05 significance threshold, with 4 genes above it.

![Genome-wide fusion-hotspot summary plot](manhattan.svg)

### FDR-significant and honorable-mention genes

| gene_symbol | tier | n_events_analyzed | in_frame_percent | domain_retention_percent | fisher_p_value | min_fdr_adjusted_q_value |
|---|---|---|---|---|---|---|
| NTRK3 | FDR-significant | 58 | 94.83 | 94.83 | 0.06955 | 8.165e-09 |
| ROS1 | FDR-significant | 122 | 59.02 | 87.7 | 0.214 | 0.0002213 |
| ETV6 | FDR-significant | 90 | 71.11 | 75.56 | 5.123e-06 | 0.002333 |
| FLI1 | FDR-significant | 118 | 61.86 | 5.932 | 0.9888 | 0.008206 |
| RET | Honorable mention | 194 | 75.26 | 92.27 | 0.0004197 | 0.09669 |
| FGFR2 | Honorable mention | 136 | 79.41 | 84.56 | 0.0004247 | 0.09669 |
| ALK | Honorable mention | 272 | 83.09 | 95.96 | 0.001917 | 0.2081 |
| EGFR | Honorable mention | 55 | 50.91 | 74.55 | 0.01145 | 0.2335 |
| BRAF | Honorable mention | 179 | 84.36 | 91.06 | 0.01337 | 0.2646 |
| FGFR3 | Honorable mention | 152 | 48.68 | 93.42 | 0.01948 | 0.2081 |
| CDKN2B | Honorable mention | 12 | 16.67 | 16.67 | 0.02222 | 0.2081 |
| IKBKE | Honorable mention | 12 | 25 | 33.33 | 0.02424 | 0.3211 |
| TMPRSS2 | Honorable mention | 867 | 29.99 | 6.113 | 0.0292 | 0.2081 |
| NTRK1 | Honorable mention | 78 | 41.03 | 82.05 | 0.04826 | 0.2081 |
| INPPL1 | Honorable mention | 14 | 28.57 | 50 | 0.04895 | 0.4687 |
| CDK12 | Honorable mention | 43 | 16.28 | 37.21 | 0.0677 | 0.2081 |
| NAB2 | Honorable mention | 72 | 34.72 | 68.06 | 0.1126 | 0.2081 |
| KDM5A | Honorable mention | 10 | 30 | 50 | 0.119 | 0.898 |
| FH | Honorable mention | 16 | 6.25 | 12.5 | 0.1333 | 0.6863 |

### Gene highlights

#### NTRK3 (FDR-significant)

NTRK3 was analyzed across 58 fusion events, 94.8% in-frame and 94.8% domain-retained. Domain-retention Fisher's exact test on NTRK3's own fusion events alone gives p=0.0695489 (not statistically significant at alpha=0.05). Among 359 genes scanned genome-wide in this cohort run, NTRK3's Benjamini-Hochberg-adjusted q=8.1653e-09 (reaches genome-wide FDR significance).

![NTRK3 key figure](gene_reports/ntrk3/visualizations/fusion_schematic.svg)

Full per-gene detail: [gene_reports/ntrk3/report.md](gene_reports/ntrk3/report.md)

#### ROS1 (FDR-significant)

ROS1 was analyzed across 122 fusion events, 59.0% in-frame and 87.7% domain-retained. Domain-retention Fisher's exact test on ROS1's own fusion events alone gives p=0.214047 (not statistically significant at alpha=0.05). Among 359 genes scanned genome-wide in this cohort run, ROS1's Benjamini-Hochberg-adjusted q=0.000221281 (reaches genome-wide FDR significance).

![ROS1 key figure](gene_reports/ros1/visualizations/fusion_schematic.svg)

Full per-gene detail: [gene_reports/ros1/report.md](gene_reports/ros1/report.md)

#### ETV6 (FDR-significant)

ETV6 was analyzed across 90 fusion events, 71.1% in-frame and 75.6% domain-retained. Domain-retention Fisher's exact test on ETV6's own fusion events alone gives p=5.12326e-06 (statistically significant at alpha=0.05). The Sterile alpha motif (SAM)/Pointed domain appears to be required for retention. Among 359 genes scanned genome-wide in this cohort run, ETV6's Benjamini-Hochberg-adjusted q=0.00233279 (reaches genome-wide FDR significance).

![ETV6 key figure](gene_reports/etv6/visualizations/fusion_schematic.svg)

Full per-gene detail: [gene_reports/etv6/report.md](gene_reports/etv6/report.md)

#### FLI1 (FDR-significant)

FLI1 was analyzed across 118 fusion events, 61.9% in-frame and 5.9% domain-retained. Domain-retention Fisher's exact test on FLI1's own fusion events alone gives p=0.988837 (not statistically significant at alpha=0.05). Among 359 genes scanned genome-wide in this cohort run, FLI1's Benjamini-Hochberg-adjusted q=0.00820622 (reaches genome-wide FDR significance).

![FLI1 key figure](gene_reports/fli1/visualizations/fusion_schematic.svg)

Full per-gene detail: [gene_reports/fli1/report.md](gene_reports/fli1/report.md)

#### RET (Honorable mention, Curated gene config)

RET was analyzed across 194 fusion events, 75.3% in-frame and 92.3% domain-retained. Domain-retention Fisher's exact test on RET's own fusion events alone gives p=0.000419666 (statistically significant at alpha=0.05). The Protein kinase domain appears to be required for retention. Among 359 genes scanned genome-wide in this cohort run, RET's Benjamini-Hochberg-adjusted q=0.096687 (does not reach genome-wide FDR significance). Did not survive genome-wide multiple-testing correction (FDR-adjusted q-value at or above the significance threshold), but ranks highly by raw p-value among the non-FDR-significant genes and may warrant targeted follow-up. This is NOT a claim of statistical significance.

![RET key figure](gene_reports/ret/visualizations/fusion_schematic.svg)

Full per-gene detail: [gene_reports/ret/report.md](gene_reports/ret/report.md)

#### FGFR2 (Honorable mention)

FGFR2 was analyzed across 136 fusion events, 79.4% in-frame and 84.6% domain-retained. Domain-retention Fisher's exact test on FGFR2's own fusion events alone gives p=0.000424687 (statistically significant at alpha=0.05). The Protein tyrosine and serine/threonine kinase appears to be required for retention. Among 359 genes scanned genome-wide in this cohort run, FGFR2's Benjamini-Hochberg-adjusted q=0.096687 (does not reach genome-wide FDR significance). Did not survive genome-wide multiple-testing correction (FDR-adjusted q-value at or above the significance threshold), but ranks highly by raw p-value among the non-FDR-significant genes and may warrant targeted follow-up. This is NOT a claim of statistical significance.

![FGFR2 key figure](gene_reports/fgfr2/visualizations/fusion_schematic.svg)

Full per-gene detail: [gene_reports/fgfr2/report.md](gene_reports/fgfr2/report.md)

#### ALK (Honorable mention, Curated gene config)

ALK was analyzed across 272 fusion events, 83.1% in-frame and 96.0% domain-retained. Domain-retention Fisher's exact test on ALK's own fusion events alone gives p=0.00191661 (statistically significant at alpha=0.05). The Protein kinase domain appears to be required for retention. Among 359 genes scanned genome-wide in this cohort run, ALK's Benjamini-Hochberg-adjusted q=0.208073 (does not reach genome-wide FDR significance). Did not survive genome-wide multiple-testing correction (FDR-adjusted q-value at or above the significance threshold), but ranks highly by raw p-value among the non-FDR-significant genes and may warrant targeted follow-up. This is NOT a claim of statistical significance.

![ALK key figure](gene_reports/alk/visualizations/fusion_schematic.svg)

Full per-gene detail: [gene_reports/alk/report.md](gene_reports/alk/report.md)

#### EGFR (Honorable mention)

EGFR was analyzed across 55 fusion events, 50.9% in-frame and 74.5% domain-retained. Domain-retention Fisher's exact test on EGFR's own fusion events alone gives p=0.0114539 (statistically significant at alpha=0.05). The Protein tyrosine and serine/threonine kinase appears to be required for retention. Among 359 genes scanned genome-wide in this cohort run, EGFR's Benjamini-Hochberg-adjusted q=0.233522 (does not reach genome-wide FDR significance). Did not survive genome-wide multiple-testing correction (FDR-adjusted q-value at or above the significance threshold), but ranks highly by raw p-value among the non-FDR-significant genes and may warrant targeted follow-up. This is NOT a claim of statistical significance.

![EGFR key figure](gene_reports/egfr/visualizations/fusion_schematic.svg)

Full per-gene detail: [gene_reports/egfr/report.md](gene_reports/egfr/report.md)

#### BRAF (Honorable mention, Curated gene config)

BRAF was analyzed across 179 fusion events, 84.4% in-frame and 91.1% domain-retained. Domain-retention Fisher's exact test on BRAF's own fusion events alone gives p=0.0133676 (statistically significant at alpha=0.05). The Protein kinase domain appears to be required for retention. Among 359 genes scanned genome-wide in this cohort run, BRAF's Benjamini-Hochberg-adjusted q=0.264639 (does not reach genome-wide FDR significance). The RAS-binding domain and Cysteine-rich domain appear to require loss or disruption rather than retention. Did not survive genome-wide multiple-testing correction (FDR-adjusted q-value at or above the significance threshold), but ranks highly by raw p-value among the non-FDR-significant genes and may warrant targeted follow-up. This is NOT a claim of statistical significance.

![BRAF key figure](gene_reports/braf/visualizations/fusion_schematic.svg)

Full per-gene detail: [gene_reports/braf/report.md](gene_reports/braf/report.md)

#### FGFR3 (Honorable mention)

FGFR3 was analyzed across 152 fusion events, 48.7% in-frame and 93.4% domain-retained. Domain-retention Fisher's exact test on FGFR3's own fusion events alone gives p=0.0194824 (statistically significant at alpha=0.05). The Protein tyrosine and serine/threonine kinase appears to be required for retention. Among 359 genes scanned genome-wide in this cohort run, FGFR3's Benjamini-Hochberg-adjusted q=0.208073 (does not reach genome-wide FDR significance). Did not survive genome-wide multiple-testing correction (FDR-adjusted q-value at or above the significance threshold), but ranks highly by raw p-value among the non-FDR-significant genes and may warrant targeted follow-up. This is NOT a claim of statistical significance.

![FGFR3 key figure](gene_reports/fgfr3/visualizations/fusion_schematic.svg)

Full per-gene detail: [gene_reports/fgfr3/report.md](gene_reports/fgfr3/report.md)

#### CDKN2B (Honorable mention)

CDKN2B was analyzed across 12 fusion events, 16.7% in-frame and 16.7% domain-retained. Domain-retention Fisher's exact test on CDKN2B's own fusion events alone gives p=0.0222222 (statistically significant at alpha=0.05). Among 359 genes scanned genome-wide in this cohort run, CDKN2B's Benjamini-Hochberg-adjusted q=0.208073 (does not reach genome-wide FDR significance). Did not survive genome-wide multiple-testing correction (FDR-adjusted q-value at or above the significance threshold), but ranks highly by raw p-value among the non-FDR-significant genes and may warrant targeted follow-up. This is NOT a claim of statistical significance.

![CDKN2B key figure](gene_reports/cdkn2b/visualizations/fusion_schematic.svg)

Full per-gene detail: [gene_reports/cdkn2b/report.md](gene_reports/cdkn2b/report.md)

#### IKBKE (Honorable mention)

IKBKE was analyzed across 12 fusion events, 25.0% in-frame and 33.3% domain-retained. Domain-retention Fisher's exact test on IKBKE's own fusion events alone gives p=0.0242424 (statistically significant at alpha=0.05). Among 359 genes scanned genome-wide in this cohort run, IKBKE's Benjamini-Hochberg-adjusted q=0.321091 (does not reach genome-wide FDR significance). Did not survive genome-wide multiple-testing correction (FDR-adjusted q-value at or above the significance threshold), but ranks highly by raw p-value among the non-FDR-significant genes and may warrant targeted follow-up. This is NOT a claim of statistical significance.

![IKBKE key figure](gene_reports/ikbke/visualizations/fusion_schematic.svg)

Full per-gene detail: [gene_reports/ikbke/report.md](gene_reports/ikbke/report.md)

#### TMPRSS2 (Honorable mention)

TMPRSS2 was analyzed across 867 fusion events, 30.0% in-frame and 6.1% domain-retained. Domain-retention Fisher's exact test on TMPRSS2's own fusion events alone gives p=0.0292045 (statistically significant at alpha=0.05). The Trypsin appears to be required for retention. Among 359 genes scanned genome-wide in this cohort run, TMPRSS2's Benjamini-Hochberg-adjusted q=0.208073 (does not reach genome-wide FDR significance). Did not survive genome-wide multiple-testing correction (FDR-adjusted q-value at or above the significance threshold), but ranks highly by raw p-value among the non-FDR-significant genes and may warrant targeted follow-up. This is NOT a claim of statistical significance.

![TMPRSS2 key figure](gene_reports/tmprss2/visualizations/fusion_schematic.svg)

Full per-gene detail: [gene_reports/tmprss2/report.md](gene_reports/tmprss2/report.md)

#### NTRK1 (Honorable mention, Curated gene config)

NTRK1 was analyzed across 78 fusion events, 41.0% in-frame and 82.1% domain-retained. Domain-retention Fisher's exact test on NTRK1's own fusion events alone gives p=0.0482628 (statistically significant at alpha=0.05). The Protein kinase domain appears to be required for retention. Among 359 genes scanned genome-wide in this cohort run, NTRK1's Benjamini-Hochberg-adjusted q=0.208073 (does not reach genome-wide FDR significance). Did not survive genome-wide multiple-testing correction (FDR-adjusted q-value at or above the significance threshold), but ranks highly by raw p-value among the non-FDR-significant genes and may warrant targeted follow-up. This is NOT a claim of statistical significance.

![NTRK1 key figure](gene_reports/ntrk1/visualizations/fusion_schematic.svg)

Full per-gene detail: [gene_reports/ntrk1/report.md](gene_reports/ntrk1/report.md)

#### INPPL1 (Honorable mention)

INPPL1 was analyzed across 14 fusion events, 28.6% in-frame and 50.0% domain-retained. Domain-retention Fisher's exact test on INPPL1's own fusion events alone gives p=0.048951 (statistically significant at alpha=0.05). Among 359 genes scanned genome-wide in this cohort run, INPPL1's Benjamini-Hochberg-adjusted q=0.468663 (does not reach genome-wide FDR significance). Did not survive genome-wide multiple-testing correction (FDR-adjusted q-value at or above the significance threshold), but ranks highly by raw p-value among the non-FDR-significant genes and may warrant targeted follow-up. This is NOT a claim of statistical significance.

![INPPL1 key figure](gene_reports/inppl1/visualizations/fusion_schematic.svg)

Full per-gene detail: [gene_reports/inppl1/report.md](gene_reports/inppl1/report.md)

#### CDK12 (Honorable mention)

CDK12 was analyzed across 43 fusion events, 16.3% in-frame and 37.2% domain-retained. Domain-retention Fisher's exact test on CDK12's own fusion events alone gives p=0.0677006 (not statistically significant at alpha=0.05). Among 359 genes scanned genome-wide in this cohort run, CDK12's Benjamini-Hochberg-adjusted q=0.208073 (does not reach genome-wide FDR significance). Did not survive genome-wide multiple-testing correction (FDR-adjusted q-value at or above the significance threshold), but ranks highly by raw p-value among the non-FDR-significant genes and may warrant targeted follow-up. This is NOT a claim of statistical significance.

![CDK12 key figure](gene_reports/cdk12/visualizations/fusion_schematic.svg)

Full per-gene detail: [gene_reports/cdk12/report.md](gene_reports/cdk12/report.md)

#### NAB2 (Honorable mention)

NAB2 was analyzed across 72 fusion events, 34.7% in-frame and 68.1% domain-retained. Domain-retention Fisher's exact test on NAB2's own fusion events alone gives p=0.11258 (not statistically significant at alpha=0.05). Among 359 genes scanned genome-wide in this cohort run, NAB2's Benjamini-Hochberg-adjusted q=0.208073 (does not reach genome-wide FDR significance). Did not survive genome-wide multiple-testing correction (FDR-adjusted q-value at or above the significance threshold), but ranks highly by raw p-value among the non-FDR-significant genes and may warrant targeted follow-up. This is NOT a claim of statistical significance.

![NAB2 key figure](gene_reports/nab2/visualizations/fusion_schematic.svg)

Full per-gene detail: [gene_reports/nab2/report.md](gene_reports/nab2/report.md)

#### KDM5A (Honorable mention)

KDM5A was analyzed across 10 fusion events, 30.0% in-frame and 50.0% domain-retained. Domain-retention Fisher's exact test on KDM5A's own fusion events alone gives p=0.119048 (not statistically significant at alpha=0.05). Among 359 genes scanned genome-wide in this cohort run, KDM5A's Benjamini-Hochberg-adjusted q=0.898036 (does not reach genome-wide FDR significance). Did not survive genome-wide multiple-testing correction (FDR-adjusted q-value at or above the significance threshold), but ranks highly by raw p-value among the non-FDR-significant genes and may warrant targeted follow-up. This is NOT a claim of statistical significance.

![KDM5A key figure](gene_reports/kdm5a/visualizations/fusion_schematic.svg)

Full per-gene detail: [gene_reports/kdm5a/report.md](gene_reports/kdm5a/report.md)

#### FH (Honorable mention)

FH was analyzed across 16 fusion events, 6.2% in-frame and 12.5% domain-retained. Domain-retention Fisher's exact test on FH's own fusion events alone gives p=0.133333 (not statistically significant at alpha=0.05). Among 359 genes scanned genome-wide in this cohort run, FH's Benjamini-Hochberg-adjusted q=0.686319 (does not reach genome-wide FDR significance). Did not survive genome-wide multiple-testing correction (FDR-adjusted q-value at or above the significance threshold), but ranks highly by raw p-value among the non-FDR-significant genes and may warrant targeted follow-up. This is NOT a claim of statistical significance.

![FH key figure](gene_reports/fh/visualizations/fusion_schematic.svg)

Full per-gene detail: [gene_reports/fh/report.md](gene_reports/fh/report.md)

## Discussion

- Cross-gene Benjamini-Hochberg FDR correction was applied jointly across the 359 genes that produced at least one computable p-value in this scan. This reduces false-positive findings relative to testing each gene in isolation, but is a conservative correction: a real per-gene effect can fail to reach the q < 0.05 threshold once corrected across that full p-value-bearing gene set.
- 4 gene(s) used hand-curated gene configurations, while 520 gene(s) were auto-configured from Genome Nexus canonical-transcript/Pfam-domain data using a kinase/catalytic-keyword heuristic to select the tracked domain; auto-configured domains have not been manually verified the way hand-curated ones have.
- Data were retrieved live from the public msk_impact_50k_2026 cBioPortal study as of 2026-09-04T19:24:30.753442+00:00. Some cBioPortal cohorts (e.g. actively accruing clinical-sequencing panels such as MSK-IMPACT) are updated periodically, so exact counts could shift if this scan is re-run later against such a cohort.
- All findings in this manuscript are computational, hypothesis-generating candidate evidence from bioinformatic analysis of public cohort data, not validated clinical calls; a gene's presence here is not a therapeutic or diagnostic recommendation.

## Appendix: per-gene report index

| gene_symbol | tier | config_source | n_events_analyzed | min_fdr_adjusted_q_value | report |
|---|---|---|---|---|---|
| NTRK3 | FDR-significant | auto | 58 | 8.165e-09 | [gene_reports/ntrk3/report.md](gene_reports/ntrk3/report.md) |
| ROS1 | FDR-significant | auto | 122 | 0.0002213 | [gene_reports/ros1/report.md](gene_reports/ros1/report.md) |
| ETV6 | FDR-significant | auto | 90 | 0.002333 | [gene_reports/etv6/report.md](gene_reports/etv6/report.md) |
| FLI1 | FDR-significant | auto | 118 | 0.008206 | [gene_reports/fli1/report.md](gene_reports/fli1/report.md) |
| RET | Honorable mention; Curated gene config | curated | 194 | 0.09669 | [gene_reports/ret/report.md](gene_reports/ret/report.md) |
| FGFR2 | Honorable mention | auto | 136 | 0.09669 | [gene_reports/fgfr2/report.md](gene_reports/fgfr2/report.md) |
| ALK | Honorable mention; Curated gene config | curated | 272 | 0.2081 | [gene_reports/alk/report.md](gene_reports/alk/report.md) |
| EGFR | Honorable mention | auto | 55 | 0.2335 | [gene_reports/egfr/report.md](gene_reports/egfr/report.md) |
| BRAF | Honorable mention; Curated gene config | curated | 179 | 0.2646 | [gene_reports/braf/report.md](gene_reports/braf/report.md) |
| FGFR3 | Honorable mention | auto | 152 | 0.2081 | [gene_reports/fgfr3/report.md](gene_reports/fgfr3/report.md) |
| CDKN2B | Honorable mention | auto | 12 | 0.2081 | [gene_reports/cdkn2b/report.md](gene_reports/cdkn2b/report.md) |
| IKBKE | Honorable mention | auto | 12 | 0.3211 | [gene_reports/ikbke/report.md](gene_reports/ikbke/report.md) |
| TMPRSS2 | Honorable mention | auto | 867 | 0.2081 | [gene_reports/tmprss2/report.md](gene_reports/tmprss2/report.md) |
| NTRK1 | Honorable mention; Curated gene config | curated | 78 | 0.2081 | [gene_reports/ntrk1/report.md](gene_reports/ntrk1/report.md) |
| INPPL1 | Honorable mention | auto | 14 | 0.4687 | [gene_reports/inppl1/report.md](gene_reports/inppl1/report.md) |
| CDK12 | Honorable mention | auto | 43 | 0.2081 | [gene_reports/cdk12/report.md](gene_reports/cdk12/report.md) |
| NAB2 | Honorable mention | auto | 72 | 0.2081 | [gene_reports/nab2/report.md](gene_reports/nab2/report.md) |
| KDM5A | Honorable mention | auto | 10 | 0.898 | [gene_reports/kdm5a/report.md](gene_reports/kdm5a/report.md) |
| FH | Honorable mention | auto | 16 | 0.6863 | [gene_reports/fh/report.md](gene_reports/fh/report.md) |
