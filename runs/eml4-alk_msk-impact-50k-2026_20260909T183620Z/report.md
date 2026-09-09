# EML4-ALK joint-partner fusion benchmark: msk_impact_50k_2026

Retrieved from public cBioPortal on 2026-09-09.

## Method

Structural-variant records were live-fetched for ALK (the pair member(s) with a curated single-gene config) through the same cBioPortal/Genome Nexus ingestion, normalization, and breakpoint-mapping pipeline used for single-gene analysis (see `run_real_benchmark`), then pooled (deduplicated by event id) and tested via `JointPartnerMode` for whether the configured ordered pair EML4->ALK is enriched relative to a marginal-independence null.

## Results

- Structural variants returned: 322
- Eligible fusion events (determinable 5'/3' orientation): 271
- Observed EML4->ALK count: 218
- Expected under independence: 209.96
- Fisher exact test (one-sided, greater): odds ratio unavailable, p=3.91794e-08
- Enriched (p < 0.05): True

## Warnings

- [ALK] Skipped EVT-P-0055607-T02-IM6-180 (ALK-SOS1): ValueError: could not determine 5'/3' role for ALK in EVT-P-0055607-T02-IM6-180; Event_Info='Antisense Fusion'
