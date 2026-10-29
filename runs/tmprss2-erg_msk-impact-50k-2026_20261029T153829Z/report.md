# TMPRSS2-ERG joint-partner fusion benchmark: msk_impact_50k_2026

Retrieved from public cBioPortal on 2026-10-29.

## Mechanism

Promoter-swap/expression-driven fusion, not domain-retention: neither TMPRSS2 nor ERG is oncogenic alone. TMPRSS2's androgen-responsive promoter drives ERG overexpression; there is no catalytic/kinase domain whose retention determines oncogenicity here, unlike the domain-retention mechanism BRAF/RET/ALK/NTRK1 (and EML4-ALK) fusions rely on.

## Method

Structural-variant genes were discovered across the study, and their records were fetched in batches from its SV molecular profile. Repeated API records were deduplicated by structural-variant ID (or full payload when absent), then normalized once. The background includes all complete, oriented protein-fusion events, without requiring transcript mapping or curated partner configs. `JointPartnerMode` tests whether TMPRSS2->ERG is enriched relative to an event-level marginal-independence null. Assay and tumor-type heterogeneity and multiple events per patient can violate this null's assumptions.

## Results

- Structural variants returned: 15770
- Eligible fusion events (determinable 5'/3' orientation): 5440
- Observed TMPRSS2->ERG count: 736
- Expected under independence: 110.50
- Fisher exact test (one-sided, greater): odds ratio 2642.86, p=0
- Enriched (p < 0.05): True

## Warnings

- Excluded 10330 of 15770 unique structural-variant records without a complete, oriented protein-fusion annotation.
- The Fisher test counts fusion events; multiple events from one sample or patient may not be independent observations.
