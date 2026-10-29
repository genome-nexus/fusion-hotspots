# FH real-data fusion benchmark: msk_impact_50k_2026

Retrieved from public cBioPortal and Genome Nexus on 2026-10-29.

## Results

- Structural variants returned for FH: 30
- Protein-fusion records found: 16
- Protein-fusion records mapped: 15
- Malformed/unmappable fusion records skipped: 1
- In-frame among known-frame events: 1/6 (16.7%)
- Unknown frame status: 9/15
- Lyase (58-389 aa) retained: 2/15 (13.3%)
- In-frame and Lyase-retained: 1/1
- Fisher exact test (one-sided): odds ratio unavailable, p=0.142857
- Breakpoint-permutation empirical p-value: 0.108911
- Contingency table `[[retained/in-frame, retained/other], [not-retained/in-frame, not-retained/other]]`: `[[1, 1], [0, 12]]`

- Mutation/CNA co-occurrence: not computed (FH has no mutual_exclusivity_targets configured; co-occurrence/mutual-exclusivity analysis was skipped.)

### Mechanistic interpretation

**Lyase retention:** nominal support not established (nominal p=0.142857, odds ratio=unavailable; q=1; does not pass FDR). No counter-pattern classification was made.

### Domain retention and discrepancies

![Domain retention diagram](visualizations/domain_retention_outliers.svg)

*Domain-retention positions for analyzed fusion events; red outlines mark reference discrepancies.*

### Fusion-transcript schematic

![FH fusion-transcript schematic](visualizations/fusion_schematic.svg)

*One row per recurrent partner/breakpoint group, sharing one amino-acid x-axis for FH's full protein length; the partner-contributed portion is colored per partner, the retained target-gene portion is colored by domain-retention status, and a red line marks the breakpoint.*

## Method

The cBioPortal `msk_impact_50k_2026_structural_variants` structural-variant profile was queried by the configured Entrez gene ID. Fusion-annotated records were adapted to the production SV schema and normalized; when `site2EffectOnFrame=NA`, frame status was resolved from `Event_Info`, not copied into `FusionEvent.Frame_status`.

FH genomic breakpoints were mapped against the Genome Nexus canonical transcript, and retention was classified against its returned Lyase coordinates. Counts are event-level with no patient deduplication. In-frame percentage uses only events explicitly called in-frame or out-of-frame; unknown-frame events are reported separately and excluded from that denominator. The Fisher comparison's `other` column combines out-of-frame and unknown-frame events, as pre-specified by the domain-retention algorithm.

For each fusion, breakpoint selection preferred the Genome Nexus canonical transcript's exon-spanned target locus over cBioPortal site labels; malformed rows with no unambiguous target-locus coordinate were skipped and listed in Warnings.

## Full-suite highlights

- Registered algorithms executed: composite_score, confidence_stats, cutpoint_detection, domain_disruption, domain_retention, exon_retention, expression_association, frequency, genomic_position_recurrence, joint_partner, junction_concentration, mechanistic_interpretation, mutation_cooccurrence, window_detection
- Cutpoint detection (descriptive): inferred breakpoint 45 aa (exon 2); label-permutation p=0.0709291, not a selection test.
- Genomic-position recurrence: The 2 events sharing protein position 45 aa use 2 distinct genomic positions spanning 101 bp -- the protein-position recurrence is broader than any single genomic hotspot, consistent with independent intronic breakpoints mapped (or clamped) onto the same protein/exon boundary rather than one shared DNA lesion.
- Top composite score: DGKH (5 events), 0.240656.
- Expression association: No mRNA expression data was available for FH in this cohort; expression-association analysis was skipped.

## Partners

ARID4B (1), DGKH (5), FMN2 (1), MAFTRR (1), NOTCH1 (1), PDE1C (1), RGS7 (2), RYR2 (1), SLC35F3 (1), TRMT1L (1), unknown (1)

## Warnings

- Skipped EVT-P-0048251-T01-IM6-30 (TRMT1L-FH): ValueError: could not determine 5'/3' role for FH in EVT-P-0048251-T01-IM6-30; Event_Info='Antisense Fusion'
- Collapsed 1 repeated observation(s) of the same fusion in 1 patient(s) before inferential tests (observation_unit='patient_gene').

## Interpretation

These values describe the live study named above.
