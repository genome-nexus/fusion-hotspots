# ROS1 real-data fusion benchmark: luad_tcga_pan_can_atlas_2018

Retrieved from public cBioPortal and Genome Nexus on 2026-10-29.

## Results

- Structural variants returned for ROS1: 10
- Protein-fusion records found: 10
- Protein-fusion records mapped: 10
- Malformed/unmappable fusion records skipped: 0
- In-frame among known-frame events: 8/8 (100.0%)
- Unknown frame status: 2/10
- Protein kinase domain (1948-2215 aa) retained: 8/10 (80.0%)
- In-frame and Protein kinase domain-retained: 8/8
- Fisher exact test (one-sided): odds ratio unavailable, p=1
- Breakpoint-permutation empirical p-value: 0.152847
- Contingency table `[[retained/in-frame, retained/other], [not-retained/in-frame, not-retained/other]]`: `[[7, 0], [0, 0]]`

- Mutation/CNA co-occurrence: not computed (ROS1 has no mutual_exclusivity_targets configured; co-occurrence/mutual-exclusivity analysis was skipped.)

### Mechanistic interpretation

**Curated mechanism:** Loss of an autoinhibitory extracellular/ligand-binding domain, with retention of an intact kinase domain -- not primarily partner-driven dimerization the way ALK/RET fusions work. The fusion deletes essentially the entire ROS1 ectodomain (the fibronectin type III repeats, YWTD beta-propellers, and N-terminal CATCH domain) while preserving the juxtamembrane segment, kinase domain, and C-terminal tail, yielding ligand-independent constitutive autophosphorylation. The partner (e.g. CD74, a type II integral membrane protein) contributes a topological role rather than a dimerization motif: because CD74's N-terminus is cytoplasmic, the ROS1-encoded transmembrane segment retained at the breakpoint is what returns the chain to the cytosol so the kinase domain can access its substrates -- positioning, not oligomerization, is what the partner supplies here.

**Protein kinase domain retention:** nominal support not established (nominal p=1, odds ratio=unavailable; FDR not assessed). No counter-pattern classification was made.

**Fibronectin type III domain disruption:** nominal support not established (nominal p=1, odds ratio=unavailable; FDR not assessed). No counter-pattern classification was made.

### Domain retention and discrepancies

![Domain retention diagram](visualizations/domain_retention_outliers.svg)

*Domain-retention positions for analyzed fusion events; red outlines mark reference discrepancies.*

### Fusion-transcript schematic

![ROS1 fusion-transcript schematic](visualizations/fusion_schematic.svg)

*One row per recurrent partner/breakpoint group, sharing one amino-acid x-axis for ROS1's full protein length; the partner-contributed portion is colored per partner, the retained target-gene portion is colored by domain-retention status, and a red line marks the breakpoint.*

## Method

The cBioPortal `luad_tcga_pan_can_atlas_2018_structural_variants` structural-variant profile was queried by the configured Entrez gene ID. Fusion-annotated records were adapted to the production SV schema and normalized; when `site2EffectOnFrame=NA`, frame status was resolved from `Event_Info`, not copied into `FusionEvent.Frame_status`.

ROS1 genomic breakpoints were mapped against the Genome Nexus canonical transcript, and retention was classified against its returned Protein kinase domain coordinates. Counts are event-level with no patient deduplication. In-frame percentage uses only events explicitly called in-frame or out-of-frame; unknown-frame events are reported separately and excluded from that denominator. The Fisher comparison's `other` column combines out-of-frame and unknown-frame events, as pre-specified by the domain-retention algorithm.

For each fusion, breakpoint selection preferred the Genome Nexus canonical transcript's exon-spanned target locus over cBioPortal site labels; malformed rows with no unambiguous target-locus coordinate were skipped and listed in Warnings.

## Full-suite highlights

- Registered algorithms executed: composite_score, confidence_stats, cutpoint_detection, domain_disruption, domain_retention, exon_retention, expression_association, frequency, genomic_position_recurrence, joint_partner, junction_concentration, mechanistic_interpretation, mutation_cooccurrence, window_detection
- Cutpoint detection: not determinable (all events share a single outcome class; no separation is possible).
- Genomic-position recurrence: All 5 events sharing protein position 1853 aa share the exact same genomic breakpoint -- consistent with a real DNA-level recurrent breakpoint, not just a protein-coordinate coincidence.
- Top composite score: CD74 (3 events), 0.269351.
- Expression association: ROS1: fusion-positive-vs-negative expression comparison skipped; no sample had both established SV assay coverage and measured expression, so no fusion-negative group could be defined.

## Partners

CD74 (3), CLTC (2), EZR (2), FBXO9 (1), SLC34A2 (2)

## Warnings

- Kept one of several distinct ROS1 fusions in 3 patient(s), dropping 3 so each patient counts once in inferential tests.
- Expression comparison: 510 measured samples excluded from the fusion-positive/negative denominator because SV gene coverage was not established.

## Interpretation

These values describe the live study named above.
