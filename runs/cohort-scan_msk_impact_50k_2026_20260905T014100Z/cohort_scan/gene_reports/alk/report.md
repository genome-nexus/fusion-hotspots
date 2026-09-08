# ALK real-data fusion benchmark: msk_impact_50k_2026

Retrieved from public cBioPortal and Genome Nexus on 2026-09-05.

## Results

- Structural variants returned for ALK: 322
- Protein-fusion records found: 272
- Protein-fusion records mapped: 271
- Malformed/unmappable fusion records skipped: 1
- In-frame: 226/272 (83.1%)
- PF07714 (1117-1382 aa) retained: 261/272 (96.0%)
- In-frame and PF07714-retained: 222/226
- Fisher exact test (one-sided): odds ratio 8.53846, p=0.00191661
- Breakpoint-permutation empirical p-value: 0.00990099
- Contingency table `[[retained/in-frame, retained/other], [not-retained/in-frame, not-retained/other]]`: `[[222, 39], [4, 6]]`

### Domain retention and discrepancies

![Domain retention diagram](visualizations/domain_retention_outliers.svg)

*Domain-retention positions for analyzed fusion events; red outlines mark reference discrepancies.*

### Fusion-transcript schematic

![ALK fusion-transcript schematic](visualizations/fusion_schematic.svg)

*One row per recurrent partner/breakpoint group, sharing one amino-acid x-axis for ALK's full protein length; the partner-contributed portion is colored per partner, the retained target-gene portion is colored by domain-retention status, and a red line marks the breakpoint.*

### Intragenic-deletion schematic

![ALK intragenic-deletion schematic](visualizations/intragenic_deletion_schematic.svg)

*Same-gene (Site1==Site2==ALK) intragenic-deletion-style SV records: a retained N-terminal block, a plain connector line for the deleted span, and a resumed C-terminal block.*

## Method

The cBioPortal `msk_impact_50k_2026_structural_variants` structural-variant profile was queried by the configured Entrez gene ID. Fusion-annotated records were adapted to the production SV schema and normalized; when `site2EffectOnFrame=NA`, frame status was resolved from `Event_Info`, not copied into `FusionEvent.Frame_status`.

ALK genomic breakpoints were mapped against the Genome Nexus canonical transcript, and retention was classified against its returned PF07714 coordinates. Counts are event-level with no patient deduplication. The Fisher comparison's `other` column combines out-of-frame and unknown-frame events, as pre-specified by the domain-retention algorithm.

For each fusion, breakpoint selection preferred the Genome Nexus canonical transcript's exon-spanned target locus over cBioPortal site labels; malformed rows with no unambiguous target-locus coordinate were skipped and listed in Warnings.

## Full-suite highlights

- Registered algorithms executed: composite_score, confidence_stats, cutpoint_detection, domain_disruption, domain_retention, exon_retention, frequency, joint_partner
- Cutpoint detection: inferred breakpoint 972 aa; corrected permutation p=0.00990099.
- Top composite score: EML4 (221 events), 0.496286.

## Partners

ASAP2 (1), CLTC (2), CTNNA1 (1), CTSE (1), DCTN1 (1), DIAPH2 (1), DYSF (1), EMILIN1 (1), EML4 (221), FN1 (1), GCC2 (1), H6PD (1), HMBOX1 (1), KIF5B (5), LRP1B (1), NOTCH1 (1), NRP2 (1), PICALM (3), PLEKHA7 (1), PLEKHH2 (1), PPP1CB (1), RANBP2 (3), SOS1 (1), SPTBN1 (4), SQSTM1 (2), STRN (7), TFG (1), TPM1 (2), TRAF3 (1), TTC27 (1), VCL (1), ZFPM2 (1)

## Warnings

- Skipped EVT-P-0055607-T02-IM6-180 (ALK-SOS1): ValueError: could not determine 5'/3' role for ALK in EVT-P-0055607-T02-IM6-180; Event_Info='Antisense Fusion'

## Reference comparison

![Reference comparison](visualizations/reference_comparison.svg)

*Configured reference percentages compared with this run.*

## Interpretation

EML4 is the recurrent partner (221/272 events), and PF07714 retention is 96.0%. These directions are consistent with the well-known EML4-ALK fusion pattern; this is a cohort-specific live measurement, not a comparison forced to a literature value.
