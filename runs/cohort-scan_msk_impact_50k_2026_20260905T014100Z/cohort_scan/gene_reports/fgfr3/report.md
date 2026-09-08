# FGFR3 real-data fusion benchmark: msk_impact_50k_2026

Retrieved from public cBioPortal and Genome Nexus on 2026-09-05.

## Results

- Structural variants returned for FGFR3: 195
- Protein-fusion records found: 152
- Protein-fusion records mapped: 151
- Malformed/unmappable fusion records skipped: 1
- In-frame: 74/152 (48.7%)
- PF07714 (472-748 aa) retained: 142/152 (93.4%)
- In-frame and PF07714-retained: 73/74
- Fisher exact test (one-sided): odds ratio 8.46377, p=0.0194824
- Breakpoint-permutation empirical p-value: 0.00990099
- Contingency table `[[retained/in-frame, retained/other], [not-retained/in-frame, not-retained/other]]`: `[[73, 69], [1, 8]]`

### Domain retention and discrepancies

![Domain retention diagram](visualizations/domain_retention_outliers.svg)

*Domain-retention positions for analyzed fusion events; red outlines mark reference discrepancies.*

### Fusion-transcript schematic

![FGFR3 fusion-transcript schematic](visualizations/fusion_schematic.svg)

*One row per recurrent partner/breakpoint group, sharing one amino-acid x-axis for FGFR3's full protein length; the partner-contributed portion is colored per partner, the retained target-gene portion is colored by domain-retention status, and a red line marks the breakpoint.*

### Intragenic-deletion schematic

![FGFR3 intragenic-deletion schematic](visualizations/intragenic_deletion_schematic.svg)

*Same-gene (Site1==Site2==FGFR3) intragenic-deletion-style SV records: a retained N-terminal block, a plain connector line for the deleted span, and a resumed C-terminal block.*

## Method

The cBioPortal `msk_impact_50k_2026_structural_variants` structural-variant profile was queried by the configured Entrez gene ID. Fusion-annotated records were adapted to the production SV schema and normalized; when `site2EffectOnFrame=NA`, frame status was resolved from `Event_Info`, not copied into `FusionEvent.Frame_status`.

FGFR3 genomic breakpoints were mapped against the Genome Nexus canonical transcript, and retention was classified against its returned PF07714 coordinates. Counts are event-level with no patient deduplication. The Fisher comparison's `other` column combines out-of-frame and unknown-frame events, as pre-specified by the domain-retention algorithm.

For each fusion, breakpoint selection preferred the Genome Nexus canonical transcript's exon-spanned target locus over cBioPortal site labels; malformed rows with no unambiguous target-locus coordinate were skipped and listed in Warnings.

## Full-suite highlights

- Registered algorithms executed: composite_score, confidence_stats, cutpoint_detection, domain_disruption, domain_retention, exon_retention, frequency, joint_partner
- Cutpoint detection: inferred breakpoint 677 aa; corrected permutation p=0.00990099.
- Top composite score: TACC3 (131 events), 0.492354.

## Partners

BABAM2 (1), CCDC149 (2), CTBP1 (1), DGKQ (1), GRK4 (1), JAKMIP1 (2), LPCAT2 (2), MAEA (1), NSD2 (3), PTBP1 (1), RBM20 (1), ROBO1 (1), SLIT2 (1), TACC3 (131), TNIP2 (3)

## Warnings

- Skipped EVT-P-0020683-T01-IM6-45 (FGFR3-SLIT2): ValueError: could not determine 5'/3' role for FGFR3 in EVT-P-0020683-T01-IM6-45; Event_Info='Antisense Fusion'

## Reference comparison

![Reference comparison](visualizations/reference_comparison.svg)

*Configured reference percentages compared with this run.*

## Interpretation

These values describe the live study named above.
