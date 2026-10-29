# EML4-ALK joint-partner fusion benchmark: msk_impact_50k_2026

Retrieved from public cBioPortal on 2026-10-29.

## Mechanism

Domain-retention/partner-driven-oligomerization, the best-characterized instance of the mechanism this benchmark's domain-retention framework is built around: EML4 contributes its N-terminal trimerization domain (TD, also called the HELP domain, historically described as a coiled-coil), which forces constitutive, ligand-independent self-association of the fused ALK intracellular kinase domain and drives trans-autophosphorylation. This substitutes for ALK's native ligand-induced extracellular-domain dimerization -- loss of that extracellular region is permissive (removes ALKAL1/ALKAL2 dependence) but the partner-driven oligomerization, not the loss itself, is the necessary activating event.

## Method

Structural-variant genes were discovered across the study, and their records were fetched in batches from its SV molecular profile. Repeated API records were deduplicated by structural-variant ID (or full payload when absent), then normalized once. The background includes all complete, oriented protein-fusion events, without requiring transcript mapping or curated partner configs. `JointPartnerMode` tests whether EML4->ALK is enriched relative to an event-level marginal-independence null. Assay and tumor-type heterogeneity and multiple events per patient can violate this null's assumptions.

## Results

- Structural variants returned: 15770
- Eligible fusion events (determinable 5'/3' orientation): 5440
- Observed EML4->ALK count: 218
- Expected under independence: 10.60
- Fisher exact test (one-sided, greater): odds ratio 12822.4, p=0
- Enriched (p < 0.05): True

## Warnings

- Excluded 10330 of 15770 unique structural-variant records without a complete, oriented protein-fusion annotation.
- The Fisher test counts fusion events; multiple events from one sample or patient may not be independent observations.
