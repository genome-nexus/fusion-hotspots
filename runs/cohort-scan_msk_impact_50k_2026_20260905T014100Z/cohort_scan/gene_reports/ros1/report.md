# ROS1 real-data fusion benchmark: msk_impact_50k_2026

Retrieved from public cBioPortal and Genome Nexus on 2026-09-05.

## Results

- Structural variants returned for ROS1: 207
- Protein-fusion records found: 122
- Protein-fusion records mapped: 118
- Malformed/unmappable fusion records skipped: 4
- In-frame: 72/122 (59.0%)
- PF07714 (1947-2215 aa) retained: 107/122 (87.7%)
- In-frame and PF07714-retained: 67/72
- Fisher exact test (one-sided): odds ratio 2.01, p=0.214047
- Breakpoint-permutation empirical p-value: 0.00990099
- Contingency table `[[retained/in-frame, retained/other], [not-retained/in-frame, not-retained/other]]`: `[[67, 40], [5, 6]]`

### Domain retention and discrepancies

![Domain retention diagram](visualizations/domain_retention_outliers.svg)

*Domain-retention positions for analyzed fusion events; red outlines mark reference discrepancies.*

### Fusion-transcript schematic

![ROS1 fusion-transcript schematic](visualizations/fusion_schematic.svg)

*One row per recurrent partner/breakpoint group, sharing one amino-acid x-axis for ROS1's full protein length; the partner-contributed portion is colored per partner, the retained target-gene portion is colored by domain-retention status, and a red line marks the breakpoint.*

### Intragenic-deletion schematic

![ROS1 intragenic-deletion schematic](visualizations/intragenic_deletion_schematic.svg)

*Same-gene (Site1==Site2==ROS1) intragenic-deletion-style SV records: a retained N-terminal block, a plain connector line for the deleted span, and a resumed C-terminal block.*

## Method

The cBioPortal `msk_impact_50k_2026_structural_variants` structural-variant profile was queried by the configured Entrez gene ID. Fusion-annotated records were adapted to the production SV schema and normalized; when `site2EffectOnFrame=NA`, frame status was resolved from `Event_Info`, not copied into `FusionEvent.Frame_status`.

ROS1 genomic breakpoints were mapped against the Genome Nexus canonical transcript, and retention was classified against its returned PF07714 coordinates. Counts are event-level with no patient deduplication. The Fisher comparison's `other` column combines out-of-frame and unknown-frame events, as pre-specified by the domain-retention algorithm.

For each fusion, breakpoint selection preferred the Genome Nexus canonical transcript's exon-spanned target locus over cBioPortal site labels; malformed rows with no unambiguous target-locus coordinate were skipped and listed in Warnings.

## Full-suite highlights

- Registered algorithms executed: composite_score, confidence_stats, cutpoint_detection, domain_disruption, domain_retention, exon_retention, frequency, joint_partner
- Cutpoint detection: inferred breakpoint 1693 aa; corrected permutation p=0.0959041.
- Top composite score: CD74 (55 events), 0.387383.

## Partners

AFG1L (1), CCDC30 (1), CD74 (55), CEP85L (2), CERT1 (1), COL13A1 (1), EYS (1), EZR (13), FN1 (1), GOLGB1 (1), GOPC (4), GRIK2 (1), IGF1R (1), LARS1 (1), LRIG3 (2), MYH9 (1), NETO1 (1), NKAIN2 (1), NUDCD3 (2), NUS1 (1), PRKN (1), RBPJL (1), SDC4 (17), SLC16A10 (1), SLC34A2 (1), SLC4A4 (2), STX7 (1), TFG (3), TPM3 (2), ZNF157 (1)

## Warnings

- Skipped EVT-P-0016175-T01-IM6-32 (EYS-ROS1): ValueError: could not determine 5'/3' role for ROS1 in EVT-P-0016175-T01-IM6-32; Event_Info='Antisense Fusion'
- Skipped EVT-P-0010101-T01-IM5-48 (NETO1-ROS1): ValueError: could not determine 5'/3' role for ROS1 in EVT-P-0010101-T01-IM5-48; Event_Info='Antisense fusion'
- Skipped EVT-P-0021671-T01-IM6-85 (ROS1-CD74): ValueError: could not determine 5'/3' role for ROS1 in EVT-P-0021671-T01-IM6-85; Event_Info='Antisense Fusion'
- Skipped EVT-P-0006237-T04-IM6-89 (ROS1-CERT1): ValueError: could not determine 5'/3' role for ROS1 in EVT-P-0006237-T04-IM6-89; Event_Info='Antisense Fusion'

## Reference comparison

![Reference comparison](visualizations/reference_comparison.svg)

*Configured reference percentages compared with this run.*

## Interpretation

These values describe the live study named above.
