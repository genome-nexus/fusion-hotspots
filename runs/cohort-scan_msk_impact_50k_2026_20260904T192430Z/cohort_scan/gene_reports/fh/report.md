# FH real-data fusion benchmark: msk_impact_50k_2026

Retrieved from public cBioPortal and Genome Nexus on 2026-09-04.

## Results

- Structural variants returned for FH: 30
- Protein-fusion records found: 16
- Protein-fusion records mapped: 15
- Malformed/unmappable fusion records skipped: 1
- In-frame: 1/16 (6.2%)
- PF00206 (58-389 aa) retained: 2/16 (12.5%)
- In-frame and PF00206-retained: 1/1
- Fisher exact test (one-sided): odds ratio unavailable, p=0.133333
- Breakpoint-permutation empirical p-value: 0.287129
- Contingency table `[[retained/in-frame, retained/other], [not-retained/in-frame, not-retained/other]]`: `[[1, 1], [0, 13]]`

### Domain retention and discrepancies

![Domain retention diagram](visualizations/domain_retention_outliers.svg)

*Domain-retention positions for analyzed fusion events; red outlines mark reference discrepancies.*

### Fusion-transcript schematic

![FH fusion-transcript schematic](visualizations/fusion_schematic.svg)

*One row per recurrent partner/breakpoint group, sharing one amino-acid x-axis for FH's full protein length; the partner-contributed portion is colored per partner, the retained target-gene portion is colored by domain-retention status, and a red line marks the breakpoint.*

## Method

The cBioPortal `msk_impact_50k_2026_structural_variants` structural-variant profile was queried by the configured Entrez gene ID. Fusion-annotated records were adapted to the production SV schema and normalized; when `site2EffectOnFrame=NA`, frame status was resolved from `Event_Info`, not copied into `FusionEvent.Frame_status`.

FH genomic breakpoints were mapped against the Genome Nexus canonical transcript, and retention was classified against its returned PF00206 coordinates. Counts are event-level with no patient deduplication. The Fisher comparison's `other` column combines out-of-frame and unknown-frame events, as pre-specified by the domain-retention algorithm.

For each fusion, breakpoint selection preferred the Genome Nexus canonical transcript's exon-spanned target locus over cBioPortal site labels; malformed rows with no unambiguous target-locus coordinate were skipped and listed in Warnings.

## Full-suite highlights

- Registered algorithms executed: composite_score, confidence_stats, cutpoint_detection, domain_disruption, domain_retention, exon_retention, frequency, joint_partner
- Cutpoint detection: inferred breakpoint 45 aa; corrected permutation p=0.0689311.
- Top composite score: RYR2 (1 events), 0.276127.

## Partners

ARID4B (1), DGKH (5), FMN2 (1), MAFTRR (1), NOTCH1 (1), PDE1C (1), RGS7 (2), RYR2 (1), SLC35F3 (1), TRMT1L (1), unknown (1)

## Warnings

- Skipped EVT-P-0048251-T01-IM6-30 (TRMT1L-FH): ValueError: could not determine 5'/3' role for FH in EVT-P-0048251-T01-IM6-30; Event_Info='Antisense Fusion'

## Reference comparison

![Reference comparison](visualizations/reference_comparison.svg)

*Configured reference percentages compared with this run.*

## Interpretation

These values describe the live study named above.
