# NTRK1 real-data fusion benchmark: thca_tcga_pan_can_atlas_2018

Retrieved from public cBioPortal and Genome Nexus on 2026-10-29.

## Results

- Structural variants returned for NTRK1: 6
- Protein-fusion records found: 6
- Protein-fusion records mapped: 6
- Malformed/unmappable fusion records skipped: 0
- In-frame among known-frame events: 6/6 (100.0%)
- Unknown frame status: 0/6
- Protein kinase domain (512-781 aa) retained: 6/6 (100.0%)
- In-frame and Protein kinase domain-retained: 6/6
- Fisher exact test (one-sided): odds ratio unavailable, p=1
- Breakpoint-permutation empirical p-value: 0.156843
- Contingency table `[[retained/in-frame, retained/other], [not-retained/in-frame, not-retained/other]]`: `[[6, 0], [0, 0]]`

- Mutation/CNA co-occurrence: not computed (NTRK1 has no mutual_exclusivity_targets configured; co-occurrence/mutual-exclusivity analysis was skipped.)

### Mechanistic interpretation

**Curated mechanism:** Combined loss of ligand-binding autoinhibition plus partner-driven dimerization. Wild-type NTRK1/TRKA requires NGF binding to its extracellular Ig-like domains to dimerize and autophosphorylate; engineered TRKA variants lacking that region (the TRKA-III splice variant, and a single P203A substitution within Ig1) are spontaneously dimerizing and transforming even without a fusion partner, indicating the deleted ectodomain is restraining, not merely permissive. No disruption_required_domains target is configured for this gene because no supported Pfam accession marks that specific restraining region -- a data-availability gap, not a mechanism gap. The fusion partner's own coiled-coil (e.g. LMNA's alpha-helical rod, TPM3's) then supplies the second, partner-driven half of the mechanism: forcing ligand-independent dimerization and trans-autophosphorylation of the retained kinase domain.

**Protein kinase domain retention:** nominal support not established (nominal p=1, odds ratio=unavailable; FDR not assessed). No counter-pattern classification was made.

### Domain retention and discrepancies

![Domain retention diagram](visualizations/domain_retention_outliers.svg)

*Domain-retention positions for analyzed fusion events; red outlines mark reference discrepancies.*

### Fusion-transcript schematic

![NTRK1 fusion-transcript schematic](visualizations/fusion_schematic.svg)

*One row per recurrent partner/breakpoint group, sharing one amino-acid x-axis for NTRK1's full protein length; the partner-contributed portion is colored per partner, the retained target-gene portion is colored by domain-retention status, and a red line marks the breakpoint.*

## Method

The cBioPortal `thca_tcga_pan_can_atlas_2018_structural_variants` structural-variant profile was queried by the configured Entrez gene ID. Fusion-annotated records were adapted to the production SV schema and normalized; when `site2EffectOnFrame=NA`, frame status was resolved from `Event_Info`, not copied into `FusionEvent.Frame_status`.

NTRK1 genomic breakpoints were mapped against the Genome Nexus canonical transcript, and retention was classified against its returned Protein kinase domain coordinates. Counts are event-level with no patient deduplication. In-frame percentage uses only events explicitly called in-frame or out-of-frame; unknown-frame events are reported separately and excluded from that denominator. The Fisher comparison's `other` column combines out-of-frame and unknown-frame events, as pre-specified by the domain-retention algorithm.

For each fusion, breakpoint selection preferred the Genome Nexus canonical transcript's exon-spanned target locus over cBioPortal site labels; malformed rows with no unambiguous target-locus coordinate were skipped and listed in Warnings.

## Full-suite highlights

- Registered algorithms executed: composite_score, confidence_stats, cutpoint_detection, domain_disruption, domain_retention, exon_retention, expression_association, frequency, genomic_position_recurrence, joint_partner, junction_concentration, mechanistic_interpretation, mutation_cooccurrence, window_detection
- Cutpoint detection: not determinable (all events share a single outcome class; no separation is possible).
- Genomic-position recurrence: All 5 events sharing protein position 399 aa share the exact same genomic breakpoint -- consistent with a real DNA-level recurrent breakpoint, not just a protein-coordinate coincidence.
- Top composite score: IRF2BP2 (2 events), 0.278585.
- Expression association: NTRK1: fusion-positive-vs-negative expression comparison skipped; no sample had both established SV assay coverage and measured expression, so no fusion-negative group could be defined.

## Partners

IRF2BP2 (2), SQSTM1 (1), SSBP2 (1), TFG (1), TPM3 (1)

## Warnings

- Expression comparison: 498 measured samples excluded from the fusion-positive/negative denominator because SV gene coverage was not established.

## Interpretation

These values describe the live study named above.
