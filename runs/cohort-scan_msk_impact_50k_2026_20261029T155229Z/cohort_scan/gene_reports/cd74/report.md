# CD74 real-data fusion benchmark: msk_impact_50k_2026

Retrieved from public cBioPortal and Genome Nexus on 2026-10-29.

## Results

- Structural variants returned for CD74: 74
- Protein-fusion records found: 67
- Protein-fusion records mapped: 66
- Malformed/unmappable fusion records skipped: 1
- In-frame among known-frame events: 39/57 (68.4%)
- Unknown frame status: 9/66
- CLIP, MHC2 interacting (17-129 aa) retained: 59/66 (89.4%)
- In-frame and CLIP, MHC2 interacting-retained: 36/39
- Fisher exact test (one-sided): odds ratio 2.27451, p=0.270083
- Breakpoint-permutation empirical p-value: 9.9999e-06
- Contingency table `[[retained/in-frame, retained/other], [not-retained/in-frame, not-retained/other]]`: `[[29, 17], [3, 4]]`

- Mutation/CNA co-occurrence: not computed (CD74 has no mutual_exclusivity_targets configured; co-occurrence/mutual-exclusivity analysis was skipped.)

### Mechanistic interpretation

**CLIP, MHC2 interacting retention:** nominal support not established (nominal p=0.270083, odds ratio=2.27451; q=1; does not pass FDR). No counter-pattern classification was made.

### Domain retention and discrepancies

![Domain retention diagram](visualizations/domain_retention_outliers.svg)

*Domain-retention positions for analyzed fusion events; red outlines mark reference discrepancies.*

### Fusion-transcript schematic

![CD74 fusion-transcript schematic](visualizations/fusion_schematic.svg)

*One row per recurrent partner/breakpoint group, sharing one amino-acid x-axis for CD74's full protein length; the partner-contributed portion is colored per partner, the retained target-gene portion is colored by domain-retention status, and a red line marks the breakpoint.*

## Method

The cBioPortal `msk_impact_50k_2026_structural_variants` structural-variant profile was queried by the configured Entrez gene ID. Fusion-annotated records were adapted to the production SV schema and normalized; when `site2EffectOnFrame=NA`, frame status was resolved from `Event_Info`, not copied into `FusionEvent.Frame_status`.

CD74 genomic breakpoints were mapped against the Genome Nexus canonical transcript, and retention was classified against its returned CLIP, MHC2 interacting coordinates. Counts are event-level with no patient deduplication. In-frame percentage uses only events explicitly called in-frame or out-of-frame; unknown-frame events are reported separately and excluded from that denominator. The Fisher comparison's `other` column combines out-of-frame and unknown-frame events, as pre-specified by the domain-retention algorithm.

For each fusion, breakpoint selection preferred the Genome Nexus canonical transcript's exon-spanned target locus over cBioPortal site labels; malformed rows with no unambiguous target-locus coordinate were skipped and listed in Warnings.

## Full-suite highlights

- Registered algorithms executed: composite_score, confidence_stats, cutpoint_detection, domain_disruption, domain_retention, exon_retention, expression_association, frequency, genomic_position_recurrence, joint_partner, junction_concentration, mechanistic_interpretation, mutation_cooccurrence, window_detection
- Cutpoint detection (descriptive): inferred breakpoint 148 aa (exon 5); label-permutation p=0.43956, not a selection test.
- Genomic-position recurrence: The 39 events sharing protein position 209 aa use 39 distinct genomic positions spanning 1332 bp -- the protein-position recurrence is broader than any single genomic hotspot, consistent with independent intronic breakpoints mapped (or clamped) onto the same protein/exon boundary rather than one shared DNA lesion.
- Top composite score: ROS1 (55 events), 0.658873.
- Expression association: No mRNA expression data was available for CD74 in this cohort; expression-association analysis was skipped.

## Partners

CAMK2A (2), NRG1 (9), PDE6A (1), ROS1 (55)

## Warnings

- Skipped EVT-P-0021671-T01-IM6-64 (ROS1-CD74): ValueError: could not determine 5'/3' role for CD74 in EVT-P-0021671-T01-IM6-64; Event_Info='Antisense Fusion'
- Collapsed 11 repeated observation(s) of the same fusion in 9 patient(s) before inferential tests (observation_unit='patient_gene').
- Kept one of several distinct CD74 fusions in 2 patient(s), dropping 2 so each patient counts once in inferential tests.

## Interpretation

These values describe the live study named above.
