# ETV6 real-data fusion benchmark: msk_impact_50k_2026

Retrieved from public cBioPortal and Genome Nexus on 2026-09-05.

## Results

- Structural variants returned for ETV6: 214
- Protein-fusion records found: 90
- Protein-fusion records mapped: 88
- Malformed/unmappable fusion records skipped: 2
- In-frame: 64/90 (71.1%)
- PF02198 (40-123 aa) retained: 68/90 (75.6%)
- In-frame and PF02198-retained: 58/64
- Fisher exact test (one-sided): odds ratio 13.5333, p=5.12326e-06
- Breakpoint-permutation empirical p-value: 0.00990099
- Contingency table `[[retained/in-frame, retained/other], [not-retained/in-frame, not-retained/other]]`: `[[58, 10], [6, 14]]`

### Domain retention and discrepancies

![Domain retention diagram](visualizations/domain_retention_outliers.svg)

*Domain-retention positions for analyzed fusion events; red outlines mark reference discrepancies.*

### Fusion-transcript schematic

![ETV6 fusion-transcript schematic](visualizations/fusion_schematic.svg)

*One row per recurrent partner/breakpoint group, sharing one amino-acid x-axis for ETV6's full protein length; the partner-contributed portion is colored per partner, the retained target-gene portion is colored by domain-retention status, and a red line marks the breakpoint.*

### Intragenic-deletion schematic

![ETV6 intragenic-deletion schematic](visualizations/intragenic_deletion_schematic.svg)

*Same-gene (Site1==Site2==ETV6) intragenic-deletion-style SV records: a retained N-terminal block, a plain connector line for the deleted span, and a resumed C-terminal block.*

## Method

The cBioPortal `msk_impact_50k_2026_structural_variants` structural-variant profile was queried by the configured Entrez gene ID. Fusion-annotated records were adapted to the production SV schema and normalized; when `site2EffectOnFrame=NA`, frame status was resolved from `Event_Info`, not copied into `FusionEvent.Frame_status`.

ETV6 genomic breakpoints were mapped against the Genome Nexus canonical transcript, and retention was classified against its returned PF02198 coordinates. Counts are event-level with no patient deduplication. The Fisher comparison's `other` column combines out-of-frame and unknown-frame events, as pre-specified by the domain-retention algorithm.

For each fusion, breakpoint selection preferred the Genome Nexus canonical transcript's exon-spanned target locus over cBioPortal site labels; malformed rows with no unambiguous target-locus coordinate were skipped and listed in Warnings.

## Full-suite highlights

- Registered algorithms executed: composite_score, confidence_stats, cutpoint_detection, domain_disruption, domain_retention, exon_retention, frequency, joint_partner
- Cutpoint detection: inferred breakpoint 155 aa; corrected permutation p=0.168317.
- Top composite score: NTRK3 (55 events), 0.448858.

## Partners

ABCC9 (1), AEBP2 (1), ANKS1B (1), APOLD1 (1), ARHGAP26 (1), BCL2L14 (3), BORCS5 (3), CCND2 (2), CDKN1B (1), CUX2 (1), FAM234B (1), FGFR1OP2 (1), GALNT18 (1), GNB1 (1), GRIN2B (1), HCFC1 (2), IKZF3 (1), LRP6 (1), NOL4 (1), NTRK3 (55), PDE3A (1), PLEKHG7 (1), PPIL2 (2), PTPRN2 (1), SGO2 (2), SLIT2 (2), SUFU (1)

## Warnings

- Skipped EVT-P-0055462-T01-IM6-7 (CCND2-ETV6): ValueError: could not determine 5'/3' role for ETV6 in EVT-P-0055462-T01-IM6-7; Event_Info='Antisense Fusion'
- Skipped EVT-P-0019604-T01-IM6-8 (CCND2-ETV6): ValueError: could not determine 5'/3' role for ETV6 in EVT-P-0019604-T01-IM6-8; Event_Info='Antisense Fusion'

## Reference comparison

![Reference comparison](visualizations/reference_comparison.svg)

*Configured reference percentages compared with this run.*

## Interpretation

These values describe the live study named above.
