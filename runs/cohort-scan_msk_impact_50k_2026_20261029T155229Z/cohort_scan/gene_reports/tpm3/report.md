# TPM3 real-data fusion benchmark: msk_impact_50k_2026

Retrieved from public cBioPortal and Genome Nexus on 2026-10-29.

## Results

- Structural variants returned for TPM3: 14
- Protein-fusion records found: 14
- Protein-fusion records mapped: 14
- Malformed/unmappable fusion records skipped: 0
- In-frame among known-frame events: 11/11 (100.0%)
- Unknown frame status: 3/14
- Tropomyosin (49-284 aa) retained: 9/14 (64.3%)
- In-frame and Tropomyosin-retained: 6/11
- Fisher exact test (one-sided): odds ratio 0, p=1
- Breakpoint-permutation empirical p-value: 1.99998e-05
- Contingency table `[[retained/in-frame, retained/other], [not-retained/in-frame, not-retained/other]]`: `[[6, 2], [3, 0]]`

- Mutation/CNA co-occurrence: not computed (TPM3 has no mutual_exclusivity_targets configured; co-occurrence/mutual-exclusivity analysis was skipped.)

### Mechanistic interpretation

**Tropomyosin retention:** nominal support not established (nominal p=1, odds ratio=0; q=1; does not pass FDR). No counter-pattern classification was made.

### Domain retention and discrepancies

![Domain retention diagram](visualizations/domain_retention_outliers.svg)

*Domain-retention positions for analyzed fusion events; red outlines mark reference discrepancies.*

### Fusion-transcript schematic

![TPM3 fusion-transcript schematic](visualizations/fusion_schematic.svg)

*One row per recurrent partner/breakpoint group, sharing one amino-acid x-axis for TPM3's full protein length; the partner-contributed portion is colored per partner, the retained target-gene portion is colored by domain-retention status, and a red line marks the breakpoint.*

## Method

The cBioPortal `msk_impact_50k_2026_structural_variants` structural-variant profile was queried by the configured Entrez gene ID. Fusion-annotated records were adapted to the production SV schema and normalized; when `site2EffectOnFrame=NA`, frame status was resolved from `Event_Info`, not copied into `FusionEvent.Frame_status`.

TPM3 genomic breakpoints were mapped against the Genome Nexus canonical transcript, and retention was classified against its returned Tropomyosin coordinates. Counts are event-level with no patient deduplication. In-frame percentage uses only events explicitly called in-frame or out-of-frame; unknown-frame events are reported separately and excluded from that denominator. The Fisher comparison's `other` column combines out-of-frame and unknown-frame events, as pre-specified by the domain-retention algorithm.

For each fusion, breakpoint selection preferred the Genome Nexus canonical transcript's exon-spanned target locus over cBioPortal site labels; malformed rows with no unambiguous target-locus coordinate were skipped and listed in Warnings.

## Full-suite highlights

- Registered algorithms executed: composite_score, confidence_stats, cutpoint_detection, domain_disruption, domain_retention, exon_retention, expression_association, frequency, genomic_position_recurrence, joint_partner, junction_concentration, mechanistic_interpretation, mutation_cooccurrence, window_detection
- Cutpoint detection (descriptive): inferred breakpoint 259 aa (exon 8/9 boundary); label-permutation p=0.003996, not a selection test.
- Genomic-position recurrence: The 3 events sharing protein position 259 aa use 3 distinct genomic positions spanning 425 bp -- the protein-position recurrence is broader than any single genomic hotspot, consistent with independent intronic breakpoints mapped (or clamped) onto the same protein/exon boundary rather than one shared DNA lesion.
- Top composite score: NTRK1 (12 events), 0.597054.
- Expression association: No mRNA expression data was available for TPM3 in this cohort; expression-association analysis was skipped.

## Partners

NTRK1 (12), ROS1 (2)

## Warnings

- Collapsed 3 repeated observation(s) of the same fusion in 3 patient(s) before inferential tests (observation_unit='patient_gene').

## Interpretation

These values describe the live study named above.
