# EGFR real-data fusion benchmark: msk_impact_50k_2026

Retrieved from public cBioPortal and Genome Nexus on 2026-09-05.

## Results

- Structural variants returned for EGFR: 559
- Protein-fusion records found: 55
- Protein-fusion records mapped: 55
- Malformed/unmappable fusion records skipped: 0
- In-frame: 28/55 (50.9%)
- PF07714 (713-965 aa) retained: 41/55 (74.5%)
- In-frame and PF07714-retained: 25/28
- Fisher exact test (one-sided): odds ratio 5.72917, p=0.0114539
- Breakpoint-permutation empirical p-value: 0.227723
- Contingency table `[[retained/in-frame, retained/other], [not-retained/in-frame, not-retained/other]]`: `[[25, 16], [3, 11]]`

### Domain retention and discrepancies

![Domain retention diagram](visualizations/domain_retention_outliers.svg)

*Domain-retention positions for analyzed fusion events; red outlines mark reference discrepancies.*

### Fusion-transcript schematic

![EGFR fusion-transcript schematic](visualizations/fusion_schematic.svg)

*One row per recurrent partner/breakpoint group, sharing one amino-acid x-axis for EGFR's full protein length; the partner-contributed portion is colored per partner, the retained target-gene portion is colored by domain-retention status, and a red line marks the breakpoint.*

### Intragenic-deletion schematic

![EGFR intragenic-deletion schematic](visualizations/intragenic_deletion_schematic.svg)

*Same-gene (Site1==Site2==EGFR) intragenic-deletion-style SV records: a retained N-terminal block, a plain connector line for the deleted span, and a resumed C-terminal block.*

## Method

The cBioPortal `msk_impact_50k_2026_structural_variants` structural-variant profile was queried by the configured Entrez gene ID. Fusion-annotated records were adapted to the production SV schema and normalized; when `site2EffectOnFrame=NA`, frame status was resolved from `Event_Info`, not copied into `FusionEvent.Frame_status`.

EGFR genomic breakpoints were mapped against the Genome Nexus canonical transcript, and retention was classified against its returned PF07714 coordinates. Counts are event-level with no patient deduplication. The Fisher comparison's `other` column combines out-of-frame and unknown-frame events, as pre-specified by the domain-retention algorithm.

For each fusion, breakpoint selection preferred the Genome Nexus canonical transcript's exon-spanned target locus over cBioPortal site labels; malformed rows with no unambiguous target-locus coordinate were skipped and listed in Warnings.

## Full-suite highlights

- Registered algorithms executed: composite_score, confidence_stats, cutpoint_detection, domain_disruption, domain_retention, exon_retention, frequency, joint_partner
- Cutpoint detection: inferred breakpoint 981 aa; corrected permutation p=0.128713.
- Top composite score: SEPTIN14 (18 events), 0.321946.

## Partners

BLTP3B (1), CADPS (1), CDH7 (1), CSF2RA (1), DDC (1), EEA1 (1), ELAPOR2 (1), FADD (2), GARS1 (1), KIF5B (1), LANCL2 (3), NIPSNAP2 (3), NUMA1 (1), PCDH15 (1), PDE1C (1), PKD1L1 (1), PLOD3 (1), RAD51 (2), SCAF4 (1), SEL1L (1), SEPTIN14 (18), TNRC18 (1), TNS3 (1), TUT7 (1), VOPP1 (1), VPS41 (1), VSTM2A (1), VWC2 (2), YIF1B (1), ZNF713 (1), ZPBP (1)

## Reference comparison

![Reference comparison](visualizations/reference_comparison.svg)

*Configured reference percentages compared with this run.*

## Interpretation

These values describe the live study named above.
