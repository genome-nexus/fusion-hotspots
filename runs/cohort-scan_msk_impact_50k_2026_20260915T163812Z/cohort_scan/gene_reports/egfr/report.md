# EGFR real-data fusion benchmark: msk_impact_50k_2026

Retrieved from public cBioPortal and Genome Nexus on 2026-09-15.

## Results

- Structural variants returned for EGFR: 559
- Protein-fusion records found: 55
- Protein-fusion records mapped: 55
- Malformed/unmappable fusion records skipped: 0
- In-frame among known-frame events: 28/40 (70.0%)
- Unknown frame status: 15/55
- Protein tyrosine and serine/threonine kinase (713-965 aa) retained: 41/55 (74.5%)
- In-frame and Protein tyrosine and serine/threonine kinase-retained: 25/28
- Fisher exact test (one-sided): odds ratio 5.72917, p=0.0114539
- Breakpoint-permutation empirical p-value: 0.227723
- Contingency table `[[retained/in-frame, retained/other], [not-retained/in-frame, not-retained/other]]`: `[[25, 16], [3, 11]]`

The Protein tyrosine and serine/threonine kinase appears to be required for retention.

- Mutation/CNA co-occurrence: not computed (EGFR has no mutual_exclusivity_targets configured; co-occurrence/mutual-exclusivity analysis was skipped.)

### Mechanistic interpretation

**Protein tyrosine and serine/threonine kinase retention:** statistically supported (p=0.0114539, odds ratio=5.72917), but 3/28 in-frame events (10.7%) show the opposite status.

Spread across distinct partner genes with no partner recurring 2+ times among them -- consistent with background noise or individual passenger events rather than a distinct recurrent subgroup; too weak to infer an alternate mechanism from this cohort alone.

| Event | Sample | Partner | Breakpoint (aa) | Status |
|---|---|---|---:|---|
| EVT-P-0059222-T03-IM7-500 | P-0059222-T03-IM7 | LANCL2 | 875 | disrupted |
| EVT-P-0029159-T01-IM6-516 | P-0029159-T01-IM6 | SEPTIN14 | 1055 | lost |
| EVT-P-0016009-T01-IM6-527 | P-0016009-T01-IM6 | VOPP1 | 249 | lost |

**Receptor L domain, Furin-like cysteine rich region, and Growth factor receptor domain IV disruption:** not statistically significant (p=0.991541, odds ratio=0.293706) -- too weak to say what this gene's fusions require, let alone characterize which events run counter to it.

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

EGFR genomic breakpoints were mapped against the Genome Nexus canonical transcript, and retention was classified against its returned Protein tyrosine and serine/threonine kinase coordinates. Counts are event-level with no patient deduplication. In-frame percentage uses only events explicitly called in-frame or out-of-frame; unknown-frame events are reported separately and excluded from that denominator. The Fisher comparison's `other` column combines out-of-frame and unknown-frame events, as pre-specified by the domain-retention algorithm.

For each fusion, breakpoint selection preferred the Genome Nexus canonical transcript's exon-spanned target locus over cBioPortal site labels; malformed rows with no unambiguous target-locus coordinate were skipped and listed in Warnings.

## Full-suite highlights

- Registered algorithms executed: composite_score, confidence_stats, cutpoint_detection, domain_disruption, domain_retention, exon_retention, expression_association, frequency, genomic_position_recurrence, joint_partner, mechanistic_interpretation, mutation_cooccurrence, window_detection
- Cutpoint detection: inferred breakpoint 981 aa (exon 24); corrected permutation p=0.128713.
- Genomic-position recurrence: The 16 events sharing protein position 982 aa use 15 distinct genomic positions spanning 700 bp -- the protein-position recurrence is broader than any single genomic hotspot, consistent with independent intronic breakpoints mapped (or clamped) onto the same protein/exon boundary rather than one shared DNA lesion.
- Top composite score: SEPTIN14 (18 events), 0.330593.
- Expression association: No mRNA expression data was available for EGFR in this cohort; expression-association analysis was skipped.

## Partners

BLTP3B (1), CADPS (1), CDH7 (1), CSF2RA (1), DDC (1), EEA1 (1), ELAPOR2 (1), FADD (2), GARS1 (1), KIF5B (1), LANCL2 (3), NIPSNAP2 (3), NUMA1 (1), PCDH15 (1), PDE1C (1), PKD1L1 (1), PLOD3 (1), RAD51 (2), SCAF4 (1), SEL1L (1), SEPTIN14 (18), TNRC18 (1), TNS3 (1), TUT7 (1), VOPP1 (1), VPS41 (1), VSTM2A (1), VWC2 (2), YIF1B (1), ZNF713 (1), ZPBP (1)

## Interpretation

These values describe the live study named above.
