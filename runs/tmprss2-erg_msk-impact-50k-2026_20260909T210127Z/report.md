# TMPRSS2-ERG joint-partner fusion benchmark: msk_impact_50k_2026

Retrieved from public cBioPortal on 2026-09-09.

## Mechanism

Promoter-swap/expression-driven fusion, not domain-retention: neither TMPRSS2 nor ERG is oncogenic alone. TMPRSS2's androgen-responsive promoter drives ERG overexpression; there is no catalytic/kinase domain whose retention determines oncogenicity here, unlike the domain-retention mechanism BRAF/RET/ALK/NTRK1 (and EML4-ALK) fusions rely on.

## Method

Structural-variant records were live-fetched for ERG (the pair member(s) with a curated single-gene config) through the same cBioPortal/Genome Nexus ingestion, normalization, and breakpoint-mapping pipeline used for single-gene analysis (see `run_real_benchmark`), then pooled (deduplicated by event id) and tested via `JointPartnerMode` for whether the configured ordered pair TMPRSS2->ERG is enriched relative to a marginal-independence null.

## Results

- Structural variants returned: 863
- Eligible fusion events (determinable 5'/3' orientation): 775
- Observed TMPRSS2->ERG count: 736
- Expected under independence: 720.81
- Fisher exact test (one-sided, greater): odds ratio unavailable, p=5.44455e-23
- Enriched (p < 0.05): True

## Warnings

- [ERG] Skipped EVT-P-0003294-T02-IM5-14 (ERG-PDE1C): ValueError: could not determine 5'/3' role for ERG in EVT-P-0003294-T02-IM5-14; Event_Info='Antisense fusion'
- [ERG] Skipped EVT-P-0064423-T01-IM7-23 (ERG-TMPRSS2): ValueError: could not determine 5'/3' role for ERG in EVT-P-0064423-T01-IM7-23; Event_Info='Antisense Fusion'
- [ERG] Skipped EVT-P-0005081-T01-IM5-81 (ERG-TMPRSS2): ValueError: could not determine 5'/3' role for ERG in EVT-P-0005081-T01-IM5-81; Event_Info='Antisense fusion'
- [ERG] Skipped EVT-P-0006057-T01-IM5-102 (ERG-TMPRSS2): ValueError: could not determine 5'/3' role for ERG in EVT-P-0006057-T01-IM5-102; Event_Info='Antisense fusion'
- [ERG] Skipped EVT-P-0052874-T01-IM6-105 (ERG-TMPRSS2): ValueError: could not determine 5'/3' role for ERG in EVT-P-0052874-T01-IM6-105; Event_Info='Antisense Fusion'
- [ERG] Skipped EVT-P-0050541-T01-IM6-121 (ERG-TMPRSS2): ValueError: could not determine 5'/3' role for ERG in EVT-P-0050541-T01-IM6-121; Event_Info='Antisense Fusion'
- [ERG] Skipped EVT-P-0050022-T02-IM6-143 (ERG-TMPRSS2): ValueError: could not determine 5'/3' role for ERG in EVT-P-0050022-T02-IM6-143; Event_Info='Antisense Fusion'
- [ERG] Skipped EVT-P-0012322-T01-IM5-215 (ERG-TMPRSS2): ValueError: could not determine 5'/3' role for ERG in EVT-P-0012322-T01-IM5-215; Event_Info='Antisense fusion'
- [ERG] Skipped EVT-P-0013858-T02-IM5-245 (ERG-TMPRSS2): ValueError: could not determine 5'/3' role for ERG in EVT-P-0013858-T02-IM5-245; Event_Info='Antisense fusion'
- [ERG] Skipped EVT-P-0018791-T01-IM6-324 (ERG-TMPRSS2): ValueError: could not determine 5'/3' role for ERG in EVT-P-0018791-T01-IM6-324; Event_Info='Antisense Fusion'
- [ERG] Skipped EVT-P-0018812-T01-IM6-326 (ERG-TMPRSS2): ValueError: could not determine 5'/3' role for ERG in EVT-P-0018812-T01-IM6-326; Event_Info='Antisense Fusion'
- [ERG] Skipped EVT-P-0025793-T01-IM6-424 (ERG-TMPRSS2): ValueError: could not determine 5'/3' role for ERG in EVT-P-0025793-T01-IM6-424; Event_Info='Antisense Fusion'
- [ERG] Skipped EVT-P-0028090-T01-IM6-452 (ERG-TMPRSS2): ValueError: could not determine 5'/3' role for ERG in EVT-P-0028090-T01-IM6-452; Event_Info='Antisense Fusion'
