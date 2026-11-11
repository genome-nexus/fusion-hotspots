## Tier counts (HEURISTIC)

Distinct chimeric-domain records, all cohorts:

| class | likely_intact | uncertain | not_assessed | likely_broken | total |
|---|---|---|---|---|---|
| clean_reconstitution | 1449 | 63 | 81 | 10 | 1603 |
| tandem_duplication_or_insertion | 44 | 270 | 2 | 16 | 332 |
| gap_or_deletion | 61 | 168 | 2 | 70 | 301 |
| ambiguous_junction_register | 608 | 76 | 18 | 10 | 712 |
| **all** | 2162 | 577 | 103 | 106 | 2948 |

Records per cohort and class (a record seen in several cohorts counts in each):

| cohort | class | likely_intact | uncertain | not_assessed | likely_broken | total |
|---|---|---|---|---|---|---|
| depmap_22q2 | ambiguous_junction_register | 207 | 36 | 6 | 4 | 253 |
| depmap_22q2 | clean_reconstitution | 694 | 36 | 38 | 3 | 771 |
| depmap_22q2 | gap_or_deletion | 25 | 61 | 0 | 21 | 107 |
| depmap_22q2 | tandem_duplication_or_insertion | 15 | 85 | 1 | 7 | 108 |
| msk_impact_50k_2026 | tandem_duplication_or_insertion | 0 | 2 | 0 | 0 | 2 |
| target_consolidated | ambiguous_junction_register | 19 | 15 | 3 | 1 | 38 |
| target_consolidated | clean_reconstitution | 118 | 10 | 16 | 5 | 149 |
| target_consolidated | gap_or_deletion | 13 | 35 | 2 | 16 | 66 |
| target_consolidated | tandem_duplication_or_insertion | 16 | 46 | 0 | 4 | 66 |
| target_starfusion | ambiguous_junction_register | 14 | 9 | 2 | 1 | 26 |
| target_starfusion | clean_reconstitution | 70 | 7 | 8 | 5 | 90 |
| target_starfusion | gap_or_deletion | 3 | 5 | 0 | 9 | 17 |
| target_starfusion | tandem_duplication_or_insertion | 1 | 13 | 0 | 1 | 15 |
| tcga_gtex | ambiguous_junction_register | 483 | 43 | 10 | 6 | 542 |
| tcga_gtex | clean_reconstitution | 1086 | 48 | 57 | 7 | 1198 |
| tcga_gtex | gap_or_deletion | 35 | 92 | 0 | 43 | 170 |
| tcga_gtex | tandem_duplication_or_insertion | 22 | 163 | 1 | 8 | 194 |

Screen rows (samples) per cohort and class, by their record's tier:

| cohort | class | likely_intact | uncertain | not_assessed | likely_broken | total |
|---|---|---|---|---|---|---|
| depmap_22q2 | ambiguous_junction_register | 1115 | 679 | 6 | 6 | 1806 |
| depmap_22q2 | clean_reconstitution | 8578 | 746 | 230 | 3 | 9557 |
| depmap_22q2 | gap_or_deletion | 79 | 125 | 0 | 39 | 243 |
| depmap_22q2 | tandem_duplication_or_insertion | 30 | 250 | 24 | 13 | 317 |
| msk_impact_50k_2026 | tandem_duplication_or_insertion | 0 | 2 | 0 | 0 | 2 |
| target_consolidated | ambiguous_junction_register | 73 | 27 | 5 | 1 | 106 |
| target_consolidated | clean_reconstitution | 732 | 30 | 37 | 8 | 807 |
| target_consolidated | gap_or_deletion | 19 | 68 | 4 | 25 | 116 |
| target_consolidated | tandem_duplication_or_insertion | 29 | 177 | 0 | 5 | 211 |
| target_starfusion | ambiguous_junction_register | 46 | 14 | 4 | 1 | 65 |
| target_starfusion | clean_reconstitution | 559 | 16 | 18 | 8 | 601 |
| target_starfusion | gap_or_deletion | 3 | 6 | 0 | 12 | 21 |
| target_starfusion | tandem_duplication_or_insertion | 1 | 92 | 0 | 1 | 94 |
| tcga_gtex | ambiguous_junction_register | 8496 | 349 | 591 | 11 | 9447 |
| tcga_gtex | clean_reconstitution | 22788 | 816 | 998 | 23 | 24625 |
| tcga_gtex | gap_or_deletion | 348 | 200 | 0 | 261 | 809 |
| tcga_gtex | tandem_duplication_or_insertion | 28 | 1247 | 1 | 10 | 1286 |

## Rules fired, missing inputs, metric disagreements (records)

| rule | records |
|---|---|
| B1_below_GA | 57 |
| B2_low_coverage | 52 |
| B3_point_feature_lost | 32 |
| I1_all_pass | 2162 |
| N1_missing | 103 |
| U1_parent_below_GA | 224 |
| U2_bit_drop | 153 |
| U3_partial_coverage | 171 |
| U4_point_feature_altered | 3 |
| U5_indel_in_ss | 200 |
| U6_indel_near_key | 70 |
| U7_junction_outside_hit | 143 |
| U8_domain_split | 251 |

| missing input | records |
|---|---|
| hmm | 1 |
| indel_secondary_structure | 26 |
| key_residues | 120 |

| disagreement | records |
|---|---|
| below_GA_but_key_points_retained | 97 |
| below_GA_but_no_broken_contacts | 82 |
| likely_intact_but_broken_core_contacts | 893 |
| passes_GA_but_indel_in_helix_or_strand | 228 |
| passes_GA_but_key_point_lost | 24 |
| passes_GA_with_junction_residues_in_insert_states | 220 |

## Shortlist for PR 2 structure prediction

One record per fusion+family, ordered by tier, then lowest SCHEMA E (max of the two parent structures), then GA margin, then read support. `ranked`: a row passes all six artefact flags; `fill`: a row passes all but the near-universal genomic_similarity and sample_burden flags (fired flags named); `named_case`: appended named fusions.

| # | fusion | family | class | tier | bits vs GA (parents) | cov | E 5'/3' | indel SS | key residues | max split reads / samples | cohorts | artefact flags (rows flagged) | shortlist |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | KRT82::KRT7 | PF00038 Filament | clean reconstitution | likely_intact | 388.31 / GA 40.0 (5' 343.76, 3' 403.96) | 1.0 | 0/0 | no_indel | swapped 9; retained 7 | 1 / 2 | tcga_gtex | genomic_similarity=1, sample_burden=1 | ranked;flagged:genomic_similarity,sample_burden |
| 2 | YWHAB::YWHAZ | PF00244 14-3-3 | clean reconstitution | likely_intact | 365.32 / GA 33.20 (5' 366.91, 3' 365.88) | 1.0 | 0/0 | no_indel | retained 4 | 1 / 6 | depmap_22q2, tcga_gtex | genomic_similarity=5, sample_burden=5 | ranked;flagged:genomic_similarity,sample_burden |
| 3 | PTN::MDK | PF01091 PTN_MK_C | clean reconstitution | likely_intact | 96.97 / GA 25.29 (5' 98.99, 3' 96.97) | 1.0 | 0/0 | no_indel | swapped 3; retained 1 | 3 / 12 | depmap_22q2, tcga_gtex | genomic_similarity=11, sample_burden=11 | ranked;flagged:genomic_similarity,sample_burden |
| 4 | CRIP1::CRIP2 | PF00412 LIM | clean reconstitution | likely_intact | 54.54 / GA 22.89 (5' 57.47, 3' 50.12) | 0.983 | 2/3 | no_indel | none in window | 5 / 5 | depmap_22q2, tcga_gtex | genomic_similarity=5, sample_burden=5 | ranked;flagged:genomic_similarity,sample_burden |
| 5 | AFM::ALB | PF00273 Serum_albumin | clean reconstitution | likely_intact | 181.12 / GA 34.40 (5' 177.4, 3' 195.4) | 1.0 | 25/26 | no_indel | swapped 11; retained 12 | 2 / 2 | tcga_gtex | genomic_similarity=1, sample_burden=1 | ranked;flagged:genomic_similarity,sample_burden |
| 6 | KLK3::PRSS8 | PF00089 Trypsin | ambiguous junction register | likely_intact | 221.94 / GA 20.60 (5' 227.5, 3' 247.06) | 1.0 | 34/42 | no_indel | retained 6 | 1 / 10 | tcga_gtex | genomic_similarity=9, sample_burden=9 | ranked;flagged:genomic_similarity,sample_burden |
| 7 | LPO::MPO | PF03098 An_peroxidase | clean reconstitution | likely_intact | 609.61 / GA 25.0 (5' 602.76, 3' 627.48) | 0.998 | 40/43 | no_indel | retained 22 | 3 / 8 | tcga_gtex | genomic_similarity=5, sample_burden=5 | ranked;flagged:genomic_similarity,sample_burden |
| 8 | GPC3::GPC4 | PF01153 Glypican | tandem duplication or insertion | uncertain | 585.44 / GA 26.39 (5' 610.0, 3' 788.15) | 0.895 | 0/0 | helix|loop|strand | retained 1 | 5 / 2 | depmap_22q2, tcga_gtex | genomic_similarity=1, sample_burden=1 | ranked;flagged:genomic_similarity,sample_burden |
| 9 | VIM::KRT18 | PF00038 Filament | clean reconstitution | uncertain | 354.58 / GA 40.0 (5' 382.38, 3' 370.36) | 0.994 | 0/0 | no_indel | swapped 12; retained 10 | 1 / 1 | tcga_gtex | none flagged | ranked |
| 10 | COL3A1::COL1A1 | PF01391 Collagen | gap or deletion | uncertain | 47.53 / GA 27.20 (5' 35.91, 3' 45.04) | 0.983 | 0/0 | loop | lost 2; swapped 1; retained 1 | 1 / 1 | tcga_gtex | none flagged | ranked |
| 11 | SPNS2::SPNS3 | PF07690 MFS_1 | tandem duplication or insertion | uncertain | 128.47 / GA 33.40 (5' 134.31, 3' 137.17) | 0.914 | 10/7 | helix|loop | none in window | 1 / 1 | tcga_gtex | none flagged | ranked |
| 12 | PIM3::ULK1 | PF00069 Pkinase | tandem duplication or insertion | uncertain | 147.05 / GA 31.70 (5' 209.86, 3' 240.06) | 0.653 | 45/34 | helix|loop|strand | swapped 1; retained 6 | 2 / 2 | target_consolidated, target_starfusion | none flagged | ranked |
| 13 | GNAL::GNAS | PF00503 G-alpha | clean reconstitution | not_assessed | 398.71 / GA 48.79 (5' 398.62, 3' 393.33) | 0.997 | 0/not_assessed | no_indel | lost not_assessed; altered not_assessed; swapped not_assessed; retained not_assessed | 1 / 1 | tcga_gtex | none flagged | ranked |
| 14 | COL1A1::COL18A1 | PF01391 Collagen | tandem duplication or insertion | likely_broken | 26.18 / GA 27.20 (5' 36.38, 3' 38.22) | 0.867 | 0/0 | loop | swapped 1; retained 3 | 2 / 4 | target_consolidated, target_starfusion, tcga_gtex | genomic_similarity=2, sample_burden=2 | ranked;flagged:genomic_similarity,sample_burden |
| 15 | COL1A1::COL12A1 | PF01391 Collagen | gap or deletion | likely_broken | 25.61 / GA 27.20 (5' 39.54, 3' 35.57) | 0.967 | 0/not_assessed | not_assessed | lost 2; retained 4 | 3 / 59 | target_consolidated, target_starfusion, tcga_gtex | genomic_similarity=56, sample_burden=56 | ranked;flagged:genomic_similarity,sample_burden |
| 16 | GPC6::GPC5 | PF01153 Glypican | gap or deletion | likely_broken | 10.16 / GA 26.39 (5' 840.44, 3' 618.29) | 0.083 | 0/0 | no_indel | none in window | 48 / 19 | depmap_22q2, target_consolidated, target_starfusion, tcga_gtex | genomic_similarity=18, sample_burden=18 | ranked;flagged:genomic_similarity,sample_burden |
| 17 | GPC4::GPC3 | PF01153 Glypican | gap or deletion | likely_broken | 8.07 / GA 26.39 (5' 788.15, 3' 610.0) | 0.086 | 0/0 | no_indel | retained 1 | 25 / 3 | depmap_22q2, tcga_gtex | genomic_similarity=2, sample_burden=2 | ranked;flagged:genomic_similarity,sample_burden |
| 18 | CMTM7::CMTM8 | PF01284 MARVEL | gap or deletion | likely_broken | 29.89 / GA 32.79 (5' 75.67, 3' 81.35) | 0.444 | 8/7 | helix|loop | none in window | 4 / 2 | tcga_gtex | genomic_similarity=1, sample_burden=1 | ranked;flagged:genomic_similarity,sample_burden |
| 19 | MYH10::MYH14 | PF00063 Myosin_head | clean reconstitution | likely_intact | 1013.82 / GA 33.29 (5' 1013.61, 3' 945.32) | 1.0 | 0/0 | no_indel | swapped 1; retained 3 | 1 / 1 | tcga_gtex | genomic_similarity=1, sample_burden=1 | fill;flagged:genomic_similarity,sample_burden |
| 20 | KRT13::KRT5 | PF00038 Filament | ambiguous junction register | likely_intact | 415.35 / GA 40.0 (5' 409.22, 3' 411.45) | 1.0 | 0/0 | no_indel | swapped 9; retained 6 | 7 / 1 | depmap_22q2 | genomic_similarity=1, sample_burden=1 | fill;flagged:genomic_similarity,sample_burden |
| 21 | KRT13::KRT7 | PF00038 Filament | clean reconstitution | likely_intact | 413.58 / GA 40.0 (5' 409.22, 3' 403.96) | 1.0 | 0/0 | no_indel | swapped 10; retained 6 | 4 / 1 | depmap_22q2 | genomic_similarity=1, sample_burden=1 | fill;flagged:genomic_similarity,sample_burden |
| 22 | KRT8::KRT5 | PF00038 Filament | clean reconstitution | likely_intact | 412.91 / GA 40.0 (5' 408.94, 3' 411.45) | 1.0 | 0/0 | no_indel | swapped 9; retained 8 | 8 / 46 | depmap_22q2, tcga_gtex | genomic_similarity=46, sample_burden=46 | fill;flagged:genomic_similarity,sample_burden |
| 23 | KRT19::KRT7 | PF00038 Filament | clean reconstitution | likely_intact | 408.09 / GA 40.0 (5' 407.15, 3' 403.96) | 1.0 | 0/0 | no_indel | swapped 11; retained 8 | 6 / 3 | depmap_22q2 | genomic_similarity=3, sample_burden=3 | fill;flagged:genomic_similarity,sample_burden |
| 24 | KRT5::KRT8 | PF00038 Filament | clean reconstitution | likely_intact | 407.07 / GA 40.0 (5' 411.45, 3' 408.94) | 1.0 | 0/0 | no_indel | swapped 8; retained 9 | 5 / 2 | depmap_22q2 | genomic_similarity=2, sample_burden=2 | fill;flagged:genomic_similarity,sample_burden |
| 25 | KRT5::KRT13 | PF00038 Filament | ambiguous junction register | likely_intact | 405.32 / GA 40.0 (5' 411.45, 3' 409.22) | 0.997 | 0/0 | no_indel | swapped 8; retained 7 | 3 / 1 | depmap_22q2 | genomic_similarity=1, sample_burden=1 | fill;flagged:genomic_similarity,sample_burden |
| 26 | FLI1::ERG | PF00178 Ets | tandem duplication or insertion | uncertain | 106.01 / GA 26.89 (5' 120.58, 3' 122.37) | 0.988 | 0/0 | helix|loop | swapped 1; retained 1 | 0 / 1 | msk_impact_50k_2026 | genomic_similarity=1, sample_burden=1, read_evidence=1 | named_case;flagged:genomic_similarity,sample_burden,read_evidence |
| 27 | STAT5B::STAT5A | PF02864 STAT_bind | tandem duplication or insertion | uncertain | 174.64 / GA 27.0 (5' 175.88, 3' 174.81) | 1.0 | 0/0 | loop|strand | retained 1 | 0 / 1 | msk_impact_50k_2026 | genomic_similarity=1, breakpoint_repeat=1, read_evidence=1 | named_case;flagged:genomic_similarity,breakpoint_repeat,read_evidence |

## Named cases (every record)


**FLI1::ERG** (1 record(s))

| # | fusion | family | class | tier | bits vs GA (parents) | cov | E 5'/3' | indel SS | key residues | max split reads / samples | cohorts | artefact flags (rows flagged) | shortlist |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|  | FLI1::ERG | PF00178 Ets | tandem duplication or insertion | uncertain | 106.01 / GA 26.89 (5' 120.58, 3' 122.37) | 0.988 | 0/0 | helix|loop | swapped 1; retained 1 | 0 / 1 | msk_impact_50k_2026 | genomic_similarity=1, sample_burden=1, read_evidence=1 |  |

- `chim_364c2913cb5e53bc` tandem_duplication_or_insertion: rules `U2_bit_drop;U5_indel_in_ss`; missing `-`; disagreements `passes_GA_but_indel_in_helix_or_strand;passes_GA_with_junction_residues_in_insert_states`; aligned indel +7 (IYYYDKN)/-0 cols; junction residues in HMM inserts 7/9; delta bits 5'/3' -14.57/-16.36; HMMvar 5'/3' -7.57/-9.22; junction DSSP 5' `TTTTSEE` 3' `HHHTTTT`; E_core 5'/3' 0/0

**STAT5B::STAT5A** (1 record(s))

| # | fusion | family | class | tier | bits vs GA (parents) | cov | E 5'/3' | indel SS | key residues | max split reads / samples | cohorts | artefact flags (rows flagged) | shortlist |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|  | STAT5B::STAT5A | PF02864 STAT_bind | tandem duplication or insertion | uncertain | 174.64 / GA 27.0 (5' 175.88, 3' 174.81) | 1.0 | 0/0 | loop|strand | retained 1 | 0 / 1 | msk_impact_50k_2026 | genomic_similarity=1, breakpoint_repeat=1, read_evidence=1 |  |

- `chim_b5debc0f4543c09e` tandem_duplication_or_insertion: rules `U7_junction_outside_hit;U8_domain_split;U5_indel_in_ss`; missing `-`; disagreements `passes_GA_but_indel_in_helix_or_strand`; aligned indel +16 (TKQPPQVLKTQTKFAA)/-0 cols; junction residues in HMM inserts 0/18; delta bits 5'/3' -1.23/-0.16; HMMvar 5'/3' 9.17/10.24; junction DSSP 5' `EEEEEEE` 3' `EEEEPPP`; E_core 5'/3' 0/0

**GPC4::GPC3** (5 record(s))

| # | fusion | family | class | tier | bits vs GA (parents) | cov | E 5'/3' | indel SS | key residues | max split reads / samples | cohorts | artefact flags (rows flagged) | shortlist |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|  | GPC4::GPC3 | PF01153 Glypican | gap or deletion | likely_broken | 345.41 / GA 26.39 (5' 788.15, 3' 610.0) | 0.407 | 8/0 | helix|loop|strand | lost 1; POINT LOST GPC3:Site:358-359(Cleavage) | 1 / 1 | depmap_22q2 | genomic_similarity=1, sample_burden=1 |  |
|  | GPC4::GPC3 | PF01153 Glypican | ambiguous junction register | likely_intact | 618.41 / GA 26.39 (5' 788.15, 3' 610.0) | 0.993 | 48/54 | no_indel | retained 1 | 1 / 1 | tcga_gtex | genomic_similarity=1, sample_burden=1 |  |
|  | GPC4::GPC3 | PF01153 Glypican | gap or deletion | likely_broken | 8.07 / GA 26.39 (5' 788.15, 3' 610.0) | 0.086 | 0/0 | no_indel | retained 1 | 25 / 3 | depmap_22q2, tcga_gtex | genomic_similarity=2, sample_burden=2 |  |
|  | GPC4::GPC3 | PF01153 Glypican | ambiguous junction register | likely_intact | 613.73 / GA 26.39 (5' 788.15, 3' 610.0) | 0.993 | 12/9 | no_indel | retained 1 | 2 / 1 | tcga_gtex | genomic_similarity=1, sample_burden=1, normal_recurrence=1 |  |
|  | GPC4::GPC3 | PF01153 Glypican | gap or deletion | uncertain | 490.71 / GA 26.39 (5' 788.15, 3' 610.0) | 0.831 | 3/4 | no_indel | retained 1 | 4 / 8 | depmap_22q2, target_consolidated, target_starfusion, tcga_gtex | genomic_similarity=8, sample_burden=8 |  |

- `chim_238f5f15a8848add` gap_or_deletion: rules `B3_point_feature_lost`; missing `-`; disagreements `passes_GA_but_key_point_lost;passes_GA_but_indel_in_helix_or_strand`; aligned indel +0 (-)/-242 cols; junction residues in HMM inserts 0/6; delta bits 5'/3' -442.74/-264.59; HMMvar 5'/3' -416.08/-237.97; junction DSSP 5' `HHTTPPP` 3' `HHHHHHH`; E_core 5'/3' 5/0

- `chim_6b54922b94f5d770` ambiguous_junction_register: rules `I1_all_pass`; missing `-`; disagreements `likely_intact_but_broken_core_contacts`; aligned indel +0 (-)/-0 cols; junction residues in HMM inserts 0/6; delta bits 5'/3' -169.74/8.41; HMMvar 5'/3' -169.73/8.38; junction DSSP 5' `HHHHHHH` 3' `HHHHHHH`; E_core 5'/3' 20/28

- `chim_796649504e0c72d7` gap_or_deletion: rules `B1_below_GA;B2_low_coverage`; missing `-`; disagreements `below_GA_but_key_points_retained;below_GA_but_no_broken_contacts`; aligned indel +0 (-)/-0 cols; junction residues in HMM inserts 0/6; delta bits 5'/3' -780.08/-601.93; HMMvar 5'/3' -756.01/-577.9; junction DSSP 5' `E-EE-TT` 3' `TS-SSS-`; E_core 5'/3' 0/0

- `chim_8d26111b3595034a` ambiguous_junction_register: rules `I1_all_pass`; missing `-`; disagreements `likely_intact_but_broken_core_contacts`; aligned indel +0 (-)/-0 cols; junction residues in HMM inserts 0/6; delta bits 5'/3' -174.43/3.73; HMMvar 5'/3' -174.42/3.69; junction DSSP 5' `E-EE-TT` 3' `E--TT-S`; E_core 5'/3' 0/3

- `chim_eec45d4e1f6a4fa1` gap_or_deletion: rules `U2_bit_drop;U8_domain_split`; missing `-`; disagreements `-`; aligned indel +0 (-)/-0 cols; junction residues in HMM inserts 0/6; delta bits 5'/3' -297.45/-119.29; HMMvar 5'/3' -278.95/-100.84; junction DSSP 5' `E-EE-TT` 3' `HHHHHHH`; E_core 5'/3' 0/2

**GPC6::GPC5** (4 record(s))

| # | fusion | family | class | tier | bits vs GA (parents) | cov | E 5'/3' | indel SS | key residues | max split reads / samples | cohorts | artefact flags (rows flagged) | shortlist |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|  | GPC6::GPC5 | PF01153 Glypican | gap or deletion | likely_broken | 10.16 / GA 26.39 (5' 840.44, 3' 618.29) | 0.083 | 0/0 | no_indel | none in window | 48 / 19 | depmap_22q2, target_consolidated, target_starfusion, tcga_gtex | genomic_similarity=18, sample_burden=18 |  |
|  | GPC6::GPC5 | PF01153 Glypican | gap or deletion | likely_broken | 122.1 / GA 26.39 (5' 840.44, 3' 618.29) | 0.157 | 0/0 | no_indel | none in window | 2 / 3 | depmap_22q2, target_consolidated | genomic_similarity=3, sample_burden=3 |  |
|  | GPC6::GPC5 | PF01153 Glypican | gap or deletion | uncertain | 490.85 / GA 26.39 (5' 840.44, 3' 618.29) | 0.488 | 0/0 | no_indel | none in window | 1 / 1 | tcga_gtex | genomic_similarity=1, sample_burden=1 |  |
|  | GPC6::GPC5 | PF01153 Glypican | ambiguous junction register | likely_intact | 630.94 / GA 26.39 (5' 840.44, 3' 618.29) | 0.975 | 17/10 | no_indel | none in window | 1 / 1 | depmap_22q2 | genomic_similarity=1, sample_burden=1 |  |

- `chim_1c4c8a349a54bfef` gap_or_deletion: rules `B1_below_GA;B2_low_coverage`; missing `-`; disagreements `below_GA_but_key_points_retained;below_GA_but_no_broken_contacts`; aligned indel +0 (-)/-0 cols; junction residues in HMM inserts 0/6; delta bits 5'/3' -830.29/-608.14; HMMvar 5'/3' -801.6/-579.4; junction DSSP 5' `E-EE-TT` 3' `-S---S-`; E_core 5'/3' 0/0

- `chim_34676df6784d7737` gap_or_deletion: rules `B2_low_coverage`; missing `-`; disagreements `-`; aligned indel +0 (-)/-0 cols; junction residues in HMM inserts 0/6; delta bits 5'/3' -718.34/-496.2; HMMvar 5'/3' -709.58/-487.38; junction DSSP 5' `HHHHHHH` 3' `-S---S-`; E_core 5'/3' 0/0

- `chim_9acf3992fdb462ea` gap_or_deletion: rules `U2_bit_drop;U3_partial_coverage;U7_junction_outside_hit;U8_domain_split`; missing `-`; disagreements `-`; aligned indel +0 (-)/-0 cols; junction residues in HMM inserts 0/6; delta bits 5'/3' -349.59/-127.44; HMMvar 5'/3' -343.08/-120.88; junction DSSP 5' `HHHHHHH` 3' `-S---S-`; E_core 5'/3' 0/0

- `chim_9af3c328edb433f7` ambiguous_junction_register: rules `I1_all_pass`; missing `-`; disagreements `likely_intact_but_broken_core_contacts`; aligned indel +0 (-)/-0 cols; junction residues in HMM inserts 0/6; delta bits 5'/3' -209.5/12.65; HMMvar 5'/3' -209.52/12.68; junction DSSP 5' `E-EE-TT` 3' `E--TT-S`; E_core 5'/3' 2/3

**GPC3::GPC4** (3 record(s))

| # | fusion | family | class | tier | bits vs GA (parents) | cov | E 5'/3' | indel SS | key residues | max split reads / samples | cohorts | artefact flags (rows flagged) | shortlist |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|  | GPC3::GPC4 | PF01153 Glypican | gap or deletion | uncertain | 302.39 / GA 26.39 (5' 610.0, 3' 788.15) | 0.44 | 15/14 | no_indel | retained 1 | 1 / 1 | depmap_22q2 | genomic_similarity=1, sample_burden=1 |  |
|  | GPC3::GPC4 | PF01153 Glypican | ambiguous junction register | likely_intact | 781.9 / GA 26.39 (5' 610.0, 3' 788.15) | 0.95 | 56/48 | no_indel | retained 1 | 1 / 2 | target_consolidated, target_starfusion | genomic_similarity=2, sample_burden=2 |  |
|  | GPC3::GPC4 | PF01153 Glypican | tandem duplication or insertion | uncertain | 585.44 / GA 26.39 (5' 610.0, 3' 788.15) | 0.895 | 0/0 | helix|loop|strand | retained 1 | 5 / 2 | depmap_22q2, tcga_gtex | genomic_similarity=1, sample_burden=1 |  |

- `chim_63f99b284ba1acb6` gap_or_deletion: rules `U2_bit_drop;U3_partial_coverage;U8_domain_split`; missing `-`; disagreements `-`; aligned indel +0 (-)/-0 cols; junction residues in HMM inserts 0/6; delta bits 5'/3' -307.61/-485.76; HMMvar 5'/3' -207.84/-385.95; junction DSSP 5' `HHHHHHH` 3' `HHHHHHH`; E_core 5'/3' 7/7

- `chim_bb16148d3ebf4abd` ambiguous_junction_register: rules `I1_all_pass`; missing `-`; disagreements `likely_intact_but_broken_core_contacts`; aligned indel +0 (-)/-0 cols; junction residues in HMM inserts 0/6; delta bits 5'/3' 171.9/-6.26; HMMvar 5'/3' 171.89/-6.22; junction DSSP 5' `HHHHHHH` 3' `HHHHHHH`; E_core 5'/3' 28/20

- `chim_d7f10c39b2b97221` tandem_duplication_or_insertion: rules `U2_bit_drop;U7_junction_outside_hit;U8_domain_split;U5_indel_in_ss`; missing `-`; disagreements `passes_GA_but_indel_in_helix_or_strand`; aligned indel +458 (GDHLKICPQGSTCCSQEMEEKYSLQSKDDFKSVVSEQCNHLQAVFASRYKKFDEFFKELLENAEKSLNDMFVKTYGHLYMQNSELFKDLFVELKRYYVVGNVNLEEMLNDFWARLLERMFRLVNSQYHFTDEYLECVSKYTEQLKPFGDVPRKLKLQVTRAFVAARTFAQGLAVAGDVVSKVSVVNPTAQCTHALLKMIYCSHCRGLVTVKPCYNYCSNIMRGCLANQGDLDFEWNNFIDAMLMVAERLEGPFNIESVMDPIDVKISDAIMNMQDNSVQVSQKVFQGCGPPKPLPAGRISRSISESAFSARFRPHHPEERPTTAAGTSLDRLVTDVKEKLKQAKKFWSSLPSNVCNDERMAAGNGNEDDCWNGKGKSRYLFAVTGNGLANQGNNPEVQVDTSKPDILILRQIMALRVMTSKMKNAYNGNDVDFFDISDESSGEGSGSGCEYQQCPSEF)/-4 cols; junction residues in HMM inserts 0/460; delta bits 5'/3' -24.56/-202.72; HMMvar 5'/3' 717.19/539.08; junction DSSP 5' `-TTS-SS` 3' `EE-TT-S`; E_core 5'/3' 0/0

**LPO::MPO** (2 record(s))

| # | fusion | family | class | tier | bits vs GA (parents) | cov | E 5'/3' | indel SS | key residues | max split reads / samples | cohorts | artefact flags (rows flagged) | shortlist |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|  | LPO::MPO | PF03098 An_peroxidase | clean reconstitution | likely_intact | 607.98 / GA 25.0 (5' 602.76, 3' 627.48) | 0.998 | 53/53 | no_indel | swapped 2; retained 20 | 2 / 1 | depmap_22q2 | paralog_sequence=1, genomic_similarity=1, sample_burden=1 |  |
|  | LPO::MPO | PF03098 An_peroxidase | clean reconstitution | likely_intact | 609.61 / GA 25.0 (5' 602.76, 3' 627.48) | 0.998 | 40/43 | no_indel | retained 22 | 3 / 8 | tcga_gtex | genomic_similarity=5, sample_burden=5 |  |

- `chim_574c2bfcae5b0cca` clean_reconstitution: rules `I1_all_pass`; missing `-`; disagreements `likely_intact_but_broken_core_contacts`; aligned indel +0 (-)/-0 cols; junction residues in HMM inserts 0/6; delta bits 5'/3' 5.22/-19.5; HMMvar 5'/3' 5.22/-19.51; junction DSSP 5' `HHHHHHI` 3' `HHHHIII`; E_core 5'/3' 36/38

- `chim_ba474d15b67899b5` clean_reconstitution: rules `I1_all_pass`; missing `-`; disagreements `likely_intact_but_broken_core_contacts`; aligned indel +0 (-)/-0 cols; junction residues in HMM inserts 0/6; delta bits 5'/3' 6.85/-17.87; HMMvar 5'/3' 6.86/-17.87; junction DSSP 5' `PEP-SST` 3' `P-SSGGG`; E_core 5'/3' 23/28

**PIM3::ULK1** (1 record(s))

| # | fusion | family | class | tier | bits vs GA (parents) | cov | E 5'/3' | indel SS | key residues | max split reads / samples | cohorts | artefact flags (rows flagged) | shortlist |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|  | PIM3::ULK1 | PF00069 Pkinase | tandem duplication or insertion | uncertain | 147.05 / GA 31.70 (5' 209.86, 3' 240.06) | 0.653 | 45/34 | helix|loop|strand | swapped 1; retained 6 | 2 / 2 | target_consolidated, target_starfusion | none flagged |  |

- `chim_2d330991737cf02a` tandem_duplication_or_insertion: rules `U2_bit_drop;U3_partial_coverage;U7_junction_outside_hit;U8_domain_split;U5_indel_in_ss`; missing `-`; disagreements `passes_GA_but_indel_in_helix_or_strand;passes_GA_with_junction_residues_in_insert_states`; aligned indel +77 (AMRTLSEDTIRLFLQQIAGAMRLLHSKGIIHRDLKPQNILLSNPAGRRANPNSIRVKIADFGFARYLQSNMMAATLC)/-0 cols; junction residues in HMM inserts 9/79; delta bits 5'/3' -62.81/-93.01; HMMvar 5'/3' 58.91/28.86; junction DSSP 5' `----S-G` 3' `HHHHSS-`; E_core 5'/3' 31/28

## Data coverage

- Records: 2948 from 50226 translated-class screen rows; rebuild problems: {'parent CDS has an internal stop codon (selenocysteine UGA?)': 113}.
- Pfam 39.0 (2026-10-01T00:00:00Z), 408 family HMMs, unavailable: 1.
- UniProt 2026_04: 1461 genes, 41 without a reviewed entry.
- AlphaFold DB model versions ['6']: 1420 accessions, 8 without a usable model, 0 fetch/parse problems.
