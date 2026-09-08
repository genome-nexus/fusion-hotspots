# NTRK1 real-data fusion benchmark: msk_impact_50k_2026

Retrieved from public cBioPortal and Genome Nexus on 2026-09-05.

## Results

- Structural variants returned for NTRK1: 118
- Protein-fusion records found: 78
- Protein-fusion records mapped: 76
- Malformed/unmappable fusion records skipped: 2
- In-frame: 32/78 (41.0%)
- PF07714 (512-781 aa) retained: 64/78 (82.1%)
- In-frame and PF07714-retained: 30/32
- Fisher exact test (one-sided): odds ratio 4.41176, p=0.0482628
- Breakpoint-permutation empirical p-value: 0.00990099
- Contingency table `[[retained/in-frame, retained/other], [not-retained/in-frame, not-retained/other]]`: `[[30, 34], [2, 10]]`

### Domain retention and discrepancies

![Domain retention diagram](visualizations/domain_retention_outliers.svg)

*Domain-retention positions for analyzed fusion events; red outlines mark reference discrepancies.*

### Fusion-transcript schematic

![NTRK1 fusion-transcript schematic](visualizations/fusion_schematic.svg)

*One row per recurrent partner/breakpoint group, sharing one amino-acid x-axis for NTRK1's full protein length; the partner-contributed portion is colored per partner, the retained target-gene portion is colored by domain-retention status, and a red line marks the breakpoint.*

### Intragenic-deletion schematic

![NTRK1 intragenic-deletion schematic](visualizations/intragenic_deletion_schematic.svg)

*Same-gene (Site1==Site2==NTRK1) intragenic-deletion-style SV records: a retained N-terminal block, a plain connector line for the deleted span, and a resumed C-terminal block.*

## Method

The cBioPortal `msk_impact_50k_2026_structural_variants` structural-variant profile was queried by the configured Entrez gene ID. Fusion-annotated records were adapted to the production SV schema and normalized; when `site2EffectOnFrame=NA`, frame status was resolved from `Event_Info`, not copied into `FusionEvent.Frame_status`.

NTRK1 genomic breakpoints were mapped against the Genome Nexus canonical transcript, and retention was classified against its returned PF07714 coordinates. Counts are event-level with no patient deduplication. The Fisher comparison's `other` column combines out-of-frame and unknown-frame events, as pre-specified by the domain-retention algorithm.

For each fusion, breakpoint selection preferred the Genome Nexus canonical transcript's exon-spanned target locus over cBioPortal site labels; malformed rows with no unambiguous target-locus coordinate were skipped and listed in Warnings.

## Full-suite highlights

- Registered algorithms executed: composite_score, confidence_stats, cutpoint_detection, domain_disruption, domain_retention, exon_retention, frequency, joint_partner
- Cutpoint detection: inferred breakpoint 501 aa; corrected permutation p=0.019802.
- Top composite score: LMNA (12 events), 0.292333.

## Partners

AFAP1 (1), ANKRD36 (1), ARHGEF11 (2), ATP1A2 (1), BCAN (2), CADM3 (1), COP1 (1), CTRC (2), DDR2 (1), DIAPH1 (1), EML4 (1), EPS15 (1), F11R (1), GON4L (1), IQGAP3 (1), IRF2BP2 (3), KIF21B (1), LMNA (12), LTAP1 (1), MEF2D (1), METTL25B (1), NELFCD (1), NOS1AP (1), P2RY8 (1), PEAR1 (1), PLEKHA6 (3), PRCC (1), RAB25 (1), RABGAP1L (1), SCP2 (1), SEMA4A (1), SHPRH (1), SLAMF6 (3), SMYD2 (1), STK11 (1), TAFA2 (1), TARS2 (1), TPM3 (12), TPR (5), TRIM63 (1), TRPM8 (1), VANGL2 (1), ZBTB7B (1)

## Warnings

- Skipped EVT-P-0022046-T01-IM6-34 (METTL25B-NTRK1): ValueError: could not determine 5'/3' role for NTRK1 in EVT-P-0022046-T01-IM6-34; Event_Info='Antisense Fusion'
- Skipped EVT-P-0008246-T02-IM5-40 (NTRK1-DDR2): ValueError: could not determine 5'/3' role for NTRK1 in EVT-P-0008246-T02-IM5-40; Event_Info='Antisense fusion'

## Reference comparison

![Reference comparison](visualizations/reference_comparison.svg)

*Configured reference percentages compared with this run.*

## Interpretation

LMNA and TPM3 each occur in 12/78 and 12/78 events, respectively; 30/32 in-frame events retain PF07714. This is directionally consistent with LMNA-NTRK1/TPM3-NTRK1 biology, while the all-event in-frame percentage is reported as observed rather than treated as a literature replication.
