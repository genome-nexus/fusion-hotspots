# PLK2 real-data fusion benchmark: msk_impact_50k_2026

Retrieved from public cBioPortal and Genome Nexus on 2026-10-29.

## Results

- Structural variants returned for PLK2: 22
- Protein-fusion records found: 7
- Protein-fusion records mapped: 7
- Malformed/unmappable fusion records skipped: 0
- In-frame among known-frame events: 2/4 (50.0%)
- Unknown frame status: 3/7
- Protein kinase domain (83-334 aa) retained: 4/7 (57.1%)
- In-frame and Protein kinase domain-retained: 2/2
- Fisher exact test (one-sided): odds ratio unavailable, p=0.285714
- Breakpoint-permutation empirical p-value: 0.237624
- Contingency table `[[retained/in-frame, retained/other], [not-retained/in-frame, not-retained/other]]`: `[[2, 2], [0, 3]]`

- Mutation/CNA co-occurrence: not computed (PLK2 has no mutual_exclusivity_targets configured; co-occurrence/mutual-exclusivity analysis was skipped.)

### Mechanistic interpretation

**Protein kinase domain retention:** nominal support not established (nominal p=0.285714, odds ratio=unavailable; q=1; does not pass FDR). No counter-pattern classification was made.

### Domain retention and discrepancies

![Domain retention diagram](visualizations/domain_retention_outliers.svg)

*Domain-retention positions for analyzed fusion events; red outlines mark reference discrepancies.*

### Fusion-transcript schematic

![PLK2 fusion-transcript schematic](visualizations/fusion_schematic.svg)

*One row per recurrent partner/breakpoint group, sharing one amino-acid x-axis for PLK2's full protein length; the partner-contributed portion is colored per partner, the retained target-gene portion is colored by domain-retention status, and a red line marks the breakpoint.*

## Method

The cBioPortal `msk_impact_50k_2026_structural_variants` structural-variant profile was queried by the configured Entrez gene ID. Fusion-annotated records were adapted to the production SV schema and normalized; when `site2EffectOnFrame=NA`, frame status was resolved from `Event_Info`, not copied into `FusionEvent.Frame_status`.

PLK2 genomic breakpoints were mapped against the Genome Nexus canonical transcript, and retention was classified against its returned Protein kinase domain coordinates. Counts are event-level with no patient deduplication. In-frame percentage uses only events explicitly called in-frame or out-of-frame; unknown-frame events are reported separately and excluded from that denominator. The Fisher comparison's `other` column combines out-of-frame and unknown-frame events, as pre-specified by the domain-retention algorithm.

For each fusion, breakpoint selection preferred the Genome Nexus canonical transcript's exon-spanned target locus over cBioPortal site labels; malformed rows with no unambiguous target-locus coordinate were skipped and listed in Warnings.

## Full-suite highlights

- Registered algorithms executed: composite_score, confidence_stats, cutpoint_detection, domain_disruption, domain_retention, exon_retention, expression_association, frequency, genomic_position_recurrence, joint_partner, junction_concentration, mechanistic_interpretation, mutation_cooccurrence, window_detection
- Cutpoint detection (descriptive): inferred breakpoint 238 aa (exon 5/6 boundary); label-permutation p=0.414585, not a selection test.
- Genomic-position recurrence: No protein-junction position was shared by two or more events with a known genomic breakpoint, so genomic-vs-protein clustering could not be compared.
- Top composite score: PDE4D (4 events), 0.324333.
- Expression association: No mRNA expression data was available for PLK2 in this cohort; expression-association analysis was skipped.

## Partners

MYO15A (1), PDE4D (4), PPP2R3B (1), WWOX (1)

## Interpretation

These values describe the live study named above.
