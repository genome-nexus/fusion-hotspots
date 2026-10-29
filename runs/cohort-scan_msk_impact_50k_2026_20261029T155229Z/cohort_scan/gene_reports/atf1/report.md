# ATF1 real-data fusion benchmark: msk_impact_50k_2026

Retrieved from public cBioPortal and Genome Nexus on 2026-10-29.

## Results

- Structural variants returned for ATF1: 30
- Protein-fusion records found: 30
- Protein-fusion records mapped: 30
- Malformed/unmappable fusion records skipped: 0
- In-frame among known-frame events: 21/25 (84.0%)
- Unknown frame status: 5/30
- bZIP transcription factor (212-268 aa) retained: 27/30 (90.0%)
- In-frame and bZIP transcription factor-retained: 19/21
- Fisher exact test (one-sided): odds ratio 1.0625, p=0.704212
- Breakpoint-permutation empirical p-value: 0.158416
- Contingency table `[[retained/in-frame, retained/other], [not-retained/in-frame, not-retained/other]]`: `[[17, 8], [2, 1]]`

- Mutation/CNA co-occurrence: not computed (ATF1 has no mutual_exclusivity_targets configured; co-occurrence/mutual-exclusivity analysis was skipped.)

### Mechanistic interpretation

**bZIP transcription factor retention:** nominal support not established (nominal p=0.704212, odds ratio=1.0625; q=1; does not pass FDR). No counter-pattern classification was made.

**pKID domain disruption:** nominal support not established (nominal p=0.547619, odds ratio=2.25; q=1; does not pass FDR). No counter-pattern classification was made.

### Domain retention and discrepancies

![Domain retention diagram](visualizations/domain_retention_outliers.svg)

*Domain-retention positions for analyzed fusion events; red outlines mark reference discrepancies.*

### Fusion-transcript schematic

![ATF1 fusion-transcript schematic](visualizations/fusion_schematic.svg)

*One row per recurrent partner/breakpoint group, sharing one amino-acid x-axis for ATF1's full protein length; the partner-contributed portion is colored per partner, the retained target-gene portion is colored by domain-retention status, and a red line marks the breakpoint.*

## Method

The cBioPortal `msk_impact_50k_2026_structural_variants` structural-variant profile was queried by the configured Entrez gene ID. Fusion-annotated records were adapted to the production SV schema and normalized; when `site2EffectOnFrame=NA`, frame status was resolved from `Event_Info`, not copied into `FusionEvent.Frame_status`.

ATF1 genomic breakpoints were mapped against the Genome Nexus canonical transcript, and retention was classified against its returned bZIP transcription factor coordinates. Counts are event-level with no patient deduplication. In-frame percentage uses only events explicitly called in-frame or out-of-frame; unknown-frame events are reported separately and excluded from that denominator. The Fisher comparison's `other` column combines out-of-frame and unknown-frame events, as pre-specified by the domain-retention algorithm.

For each fusion, breakpoint selection preferred the Genome Nexus canonical transcript's exon-spanned target locus over cBioPortal site labels; malformed rows with no unambiguous target-locus coordinate were skipped and listed in Warnings.

## Full-suite highlights

- Registered algorithms executed: composite_score, confidence_stats, cutpoint_detection, domain_disruption, domain_retention, exon_retention, expression_association, frequency, genomic_position_recurrence, joint_partner, junction_concentration, mechanistic_interpretation, mutation_cooccurrence, window_detection
- Cutpoint detection (descriptive): inferred breakpoint 31 aa (exon 2); label-permutation p=0.036963, not a selection test.
- Genomic-position recurrence: The 17 events sharing protein position 65 aa use 17 distinct genomic positions spanning 12207 bp -- the protein-position recurrence is broader than any single genomic hotspot, consistent with independent intronic breakpoints mapped (or clamped) onto the same protein/exon boundary rather than one shared DNA lesion.
- Top composite score: EWSR1 (29 events), 0.454937.
- Expression association: No mRNA expression data was available for ATF1 in this cohort; expression-association analysis was skipped.

## Partners

ERBB3 (1), EWSR1 (29)

## Warnings

- Collapsed 2 repeated observation(s) of the same fusion in 2 patient(s) before inferential tests (observation_unit='patient_gene').

## Interpretation

These values describe the live study named above.
