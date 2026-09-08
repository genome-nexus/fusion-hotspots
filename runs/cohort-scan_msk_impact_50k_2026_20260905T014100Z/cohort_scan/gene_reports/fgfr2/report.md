# FGFR2 real-data fusion benchmark: msk_impact_50k_2026

Retrieved from public cBioPortal and Genome Nexus on 2026-09-05.

## Results

- Structural variants returned for FGFR2: 195
- Protein-fusion records found: 136
- Protein-fusion records mapped: 134
- Malformed/unmappable fusion records skipped: 2
- In-frame: 108/136 (79.4%)
- PF07714 (481-757 aa) retained: 115/136 (84.6%)
- In-frame and PF07714-retained: 99/108
- Fisher exact test (one-sided): odds ratio 6.875, p=0.000424687
- Breakpoint-permutation empirical p-value: 0.00990099
- Contingency table `[[retained/in-frame, retained/other], [not-retained/in-frame, not-retained/other]]`: `[[99, 16], [9, 10]]`

### Domain retention and discrepancies

![Domain retention diagram](visualizations/domain_retention_outliers.svg)

*Domain-retention positions for analyzed fusion events; red outlines mark reference discrepancies.*

### Fusion-transcript schematic

![FGFR2 fusion-transcript schematic](visualizations/fusion_schematic.svg)

*One row per recurrent partner/breakpoint group, sharing one amino-acid x-axis for FGFR2's full protein length; the partner-contributed portion is colored per partner, the retained target-gene portion is colored by domain-retention status, and a red line marks the breakpoint.*

## Method

The cBioPortal `msk_impact_50k_2026_structural_variants` structural-variant profile was queried by the configured Entrez gene ID. Fusion-annotated records were adapted to the production SV schema and normalized; when `site2EffectOnFrame=NA`, frame status was resolved from `Event_Info`, not copied into `FusionEvent.Frame_status`.

FGFR2 genomic breakpoints were mapped against the Genome Nexus canonical transcript, and retention was classified against its returned PF07714 coordinates. Counts are event-level with no patient deduplication. The Fisher comparison's `other` column combines out-of-frame and unknown-frame events, as pre-specified by the domain-retention algorithm.

For each fusion, breakpoint selection preferred the Genome Nexus canonical transcript's exon-spanned target locus over cBioPortal site labels; malformed rows with no unambiguous target-locus coordinate were skipped and listed in Warnings.

## Full-suite highlights

- Registered algorithms executed: composite_score, confidence_stats, cutpoint_detection, domain_disruption, domain_retention, exon_retention, frequency, joint_partner
- Cutpoint detection: inferred breakpoint 767 aa; corrected permutation p=0.00990099.
- Top composite score: BICC1 (28 events), 0.337329.

## Partners

ABLIM1 (1), ADAM12 (1), ADAM9 (1), AHCYL1 (2), AMOT (1), ANKRD11 (1), ATAD1 (1), ATE1 (3), ATP6V1D (1), BAIAP2L1 (2), BICC1 (28), BICD1 (1), BNIP2 (1), BTBD16 (1), CAT (1), CCDC186 (1), CCDC6 (1), CCDC73 (1), CEACAM7 (1), CEP112 (1), CTNNA3 (3), DDX21 (1), DIABLO (1), FOXP1 (2), GAB2 (1), GRID1 (1), HSPA12A (1), INA (2), KCNH7 (1), KIAA1217 (4), KIF14 (1), KLHL29 (1), L3MBTL3 (1), LAMC1 (1), MARVELD3 (1), MBIP (1), MYH15 (1), MYH9 (1), MYLK (1), MYO18B (1), NID2 (1), NOL4 (4), NRAP (1), NRBF2 (2), NRG3 (2), OPTN (2), PAH (1), PALLD (1), PDZRN4 (1), PKP4 (1), POC1B (1), PPP1R21 (1), PRDX3 (1), PRKAG2 (1), RABGAP1L (1), RNF145 (1), ROCK1 (1), RPAP3 (1), SEPTIN10 (1), SH3GLB1 (1), SHC2 (1), SHTN1 (3), SKI (1), SORBS1 (2), SRRT (1), ST8SIA6 (1), STAT4 (1), TACC2 (8), TBC1D4 (3), TFEC (1), TTC28 (2), UBP1 (1), VCL (1), VPS26A (1), WAC (2), YTHDF3 (1), ZFAT (1), ZMYM4 (1)

## Warnings

- Skipped EVT-P-0013148-T01-IM5-128 (FGFR2-PRKAG2): ValueError: could not determine 5'/3' role for FGFR2 in EVT-P-0013148-T01-IM5-128; Event_Info='Antisense fusion'
- Skipped EVT-P-0034590-T01-IM6-138 (FGFR2-ST8SIA6): ValueError: could not determine 5'/3' role for FGFR2 in EVT-P-0034590-T01-IM6-138; Event_Info='Antisense Fusion'

## Reference comparison

![Reference comparison](visualizations/reference_comparison.svg)

*Configured reference percentages compared with this run.*

## Interpretation

These values describe the live study named above.
