# IKBKE real-data fusion benchmark: msk_impact_50k_2026

Retrieved from public cBioPortal and Genome Nexus on 2026-09-15.

## Results

- Structural variants returned for IKBKE: 28
- Protein-fusion records found: 12
- Protein-fusion records mapped: 11
- Malformed/unmappable fusion records skipped: 1
- In-frame among known-frame events: 3/6 (50.0%)
- Unknown frame status: 5/11
- Protein kinase domain (10-244 aa) retained: 4/11 (36.4%)
- In-frame and Protein kinase domain-retained: 3/3
- Fisher exact test (one-sided): odds ratio unavailable, p=0.0242424
- Breakpoint-permutation empirical p-value: 0.019802
- Contingency table `[[retained/in-frame, retained/other], [not-retained/in-frame, not-retained/other]]`: `[[3, 1], [0, 7]]`

- Mutation/CNA co-occurrence: not computed (IKBKE has no mutual_exclusivity_targets configured; co-occurrence/mutual-exclusivity analysis was skipped.)

### Mechanistic interpretation

**Protein kinase domain retention:** statistically supported (p=0.0242424, odds ratio=unavailable). All 3 in-frame events with a determinate status match it -- no counter-intuitive events observed.

### Domain retention and discrepancies

![Domain retention diagram](visualizations/domain_retention_outliers.svg)

*Domain-retention positions for analyzed fusion events; red outlines mark reference discrepancies.*

### Fusion-transcript schematic

![IKBKE fusion-transcript schematic](visualizations/fusion_schematic.svg)

*One row per recurrent partner/breakpoint group, sharing one amino-acid x-axis for IKBKE's full protein length; the partner-contributed portion is colored per partner, the retained target-gene portion is colored by domain-retention status, and a red line marks the breakpoint.*

### Intragenic-deletion schematic

![IKBKE intragenic-deletion schematic](visualizations/intragenic_deletion_schematic.svg)

*Same-gene (Site1==Site2==IKBKE) intragenic-deletion-style SV records: a retained N-terminal block, a plain connector line for the deleted span, and a resumed C-terminal block.*

## Method

The cBioPortal `msk_impact_50k_2026_structural_variants` structural-variant profile was queried by the configured Entrez gene ID. Fusion-annotated records were adapted to the production SV schema and normalized; when `site2EffectOnFrame=NA`, frame status was resolved from `Event_Info`, not copied into `FusionEvent.Frame_status`.

IKBKE genomic breakpoints were mapped against the Genome Nexus canonical transcript, and retention was classified against its returned Protein kinase domain coordinates. Counts are event-level with no patient deduplication. In-frame percentage uses only events explicitly called in-frame or out-of-frame; unknown-frame events are reported separately and excluded from that denominator. The Fisher comparison's `other` column combines out-of-frame and unknown-frame events, as pre-specified by the domain-retention algorithm.

For each fusion, breakpoint selection preferred the Genome Nexus canonical transcript's exon-spanned target locus over cBioPortal site labels; malformed rows with no unambiguous target-locus coordinate were skipped and listed in Warnings.

## Full-suite highlights

- Registered algorithms executed: composite_score, confidence_stats, cutpoint_detection, domain_disruption, domain_retention, exon_retention, expression_association, frequency, genomic_position_recurrence, joint_partner, mechanistic_interpretation, mutation_cooccurrence, window_detection
- Cutpoint detection: inferred breakpoint 359 aa (exon 10); corrected permutation p=0.39604.
- Genomic-position recurrence: The 2 events sharing protein position 447 aa use 2 distinct genomic positions spanning 2 bp -- the protein-position recurrence is broader than any single genomic hotspot, consistent with independent intronic breakpoints mapped (or clamped) onto the same protein/exon boundary rather than one shared DNA lesion.
- Top composite score: RASAL2 (1 events), 0.315016.
- Expression association: No mRNA expression data was available for IKBKE in this cohort; expression-association analysis was skipped.

## Partners

IL20 (1), NAV1 (2), NFASC (1), RABGAP1L (1), RASAL2 (1), RASSF5 (1), RBBP5 (1), SMYD2 (1), SRGAP2 (1), TCEA3 (1), THBS3 (1)

## Warnings

- Skipped EVT-P-0032949-T01-IM6-20 (IKBKE-SMYD2): ValueError: could not determine 5'/3' role for IKBKE in EVT-P-0032949-T01-IM6-20; Event_Info='Antisense Fusion'

## Interpretation

These values describe the live study named above.
