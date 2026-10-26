# EML4 real-data fusion benchmark: msk_impact_50k_2026

Retrieved from public cBioPortal and Genome Nexus on 2026-10-26.

## Results

- Structural variants returned for EML4: 225
- Protein-fusion records found: 223
- Protein-fusion records mapped: 223
- Malformed/unmappable fusion records skipped: 0
- In-frame among known-frame events: 196/201 (97.5%)
- Unknown frame status: 22/223
- HELP motif (223-298 aa) retained: 119/223 (53.4%)
- In-frame and HELP motif-retained: 95/196
- Fisher exact test (one-sided): odds ratio 0.146617, p=0.999931
- Breakpoint-permutation empirical p-value: 9.9999e-06
- Contingency table `[[retained/in-frame, retained/other], [not-retained/in-frame, not-retained/other]]`: `[[78, 19], [84, 3]]`

- Mutation/CNA co-occurrence: not computed (EML4 has no mutual_exclusivity_targets configured; co-occurrence/mutual-exclusivity analysis was skipped.)

### Mechanistic interpretation

**HELP motif retention:** nominal support not established (nominal p=0.999931, odds ratio=0.146617; q=1; does not pass FDR). No counter-pattern classification was made.

### Domain retention and discrepancies

![Domain retention diagram](visualizations/domain_retention_outliers.svg)

*Domain-retention positions for analyzed fusion events; red outlines mark reference discrepancies.*

### Fusion-transcript schematic

![EML4 fusion-transcript schematic](visualizations/fusion_schematic.svg)

*One row per recurrent partner/breakpoint group, sharing one amino-acid x-axis for EML4's full protein length; the partner-contributed portion is colored per partner, the retained target-gene portion is colored by domain-retention status, and a red line marks the breakpoint.*

## Method

The cBioPortal `msk_impact_50k_2026_structural_variants` structural-variant profile was queried by the configured Entrez gene ID. Fusion-annotated records were adapted to the production SV schema and normalized; when `site2EffectOnFrame=NA`, frame status was resolved from `Event_Info`, not copied into `FusionEvent.Frame_status`.

EML4 genomic breakpoints were mapped against the Genome Nexus canonical transcript, and retention was classified against its returned HELP motif coordinates. Counts are event-level with no patient deduplication. In-frame percentage uses only events explicitly called in-frame or out-of-frame; unknown-frame events are reported separately and excluded from that denominator. The Fisher comparison's `other` column combines out-of-frame and unknown-frame events, as pre-specified by the domain-retention algorithm.

For each fusion, breakpoint selection preferred the Genome Nexus canonical transcript's exon-spanned target locus over cBioPortal site labels; malformed rows with no unambiguous target-locus coordinate were skipped and listed in Warnings.

## Full-suite highlights

- Registered algorithms executed: composite_score, confidence_stats, cutpoint_detection, domain_disruption, domain_retention, exon_retention, expression_association, frequency, genomic_position_recurrence, joint_partner, junction_concentration, mechanistic_interpretation, mutation_cooccurrence, window_detection
- Cutpoint detection (descriptive): inferred breakpoint 223 aa (exon 6/7 boundary); label-permutation p=0.000999001, not a selection test.
- Genomic-position recurrence: The 71 events sharing protein position 223 aa use 71 distinct genomic positions spanning 15835 bp -- the protein-position recurrence is broader than any single genomic hotspot, consistent with independent intronic breakpoints mapped (or clamped) onto the same protein/exon boundary rather than one shared DNA lesion.
- Top composite score: ALK (221 events), 0.757816.
- Expression association: No mRNA expression data was available for EML4 in this cohort; expression-association analysis was skipped.

## Partners

ALK (221), NTRK1 (1), RET (1)

## Warnings

- Collapsed 37 repeated observation(s) of the same fusion in 27 patient(s) before inferential tests (observation_unit='patient_gene').
- Kept one of several distinct EML4 fusions in 2 patient(s), dropping 2 so each patient counts once in inferential tests.

## Interpretation

These values describe the live study named above.
