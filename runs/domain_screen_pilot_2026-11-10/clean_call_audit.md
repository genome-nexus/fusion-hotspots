| Measure | depmap_22q2 | target_consolidated | target_starfusion | tcga_gtex |
|---|---|---|---|---|
| clean rows | 9,567 | 819 | 602 | 24,728 |
| distinct junctions (pair + breakpoints + family) | 831 | 169 | 96 | 1,234 |
| distinct samples | 1,282 | 415 | 404 | 10,874 |
| distinct gene pairs | 541 | 89 | 79 | 840 |
| same gene (5' == 3') | 0 | 0 | 0 | 0 |
| read-through: same chrom + strand, <= 200 kb | 227 | 119 | 69 | 3,885 |
| any different-gene pair on one chrom <= 200 kb | 302 | 150 | 99 | 6,171 |
| Compara human paralogues | 9,484 | 815 | 598 | 22,707 |
| same Pfam family, not Compara paralogues | 83 | 4 | 4 | 2,021 |
| gene not in Ensembl 116 (retired symbol) | 0 | 0 | 0 | 0 |
| both ends at annotated splice sites (ONLY_REF_SPLICE) | 9,437 | 0 | 600 | 24,610 |
| includes a non-reference splice site | 130 | 0 | 2 | 118 |
| one side < 10 aa or < 10% of its domain | 1,722 | 108 | 62 | 3,428 |
| both sides contribute substantially | 7,845 | 711 | 540 | 21,300 |
| chimera >= 95% identical to one parent domain | 7,095 | 393 | 286 | 11,577 |
| paralog_sequence_flag = flagged | 8,646 | 562 | 428 | 17,877 |
| all three artefact flags pass | 0 | 0 | 0 | 11 |
| **strict**: not same gene, not read-through, both sides substantial, flags pass | 0 | 0 | 0 | 7 |
| **strict distinct junctions** | 0 | 0 | 0 | 5 |


Top clean gene pairs, depmap_22q2:

| Pair | Rows | Samples | Junctions | Relation | Read-through rows | Families |
|---|---|---|---|---|---|---|
| ACTG2::ACTG1 | 827 | 818 | 4 | compara_paralog | 0 | PF00022 |
| CALM3::CALM2 | 657 | 450 | 3 | compara_paralog | 0 | PF13499 |
| EEF1A2::EEF1A1 | 584 | 437 | 5 | compara_paralog | 0 | PF00009, PF22594 |
| HDGFL3::HDGF | 435 | 435 | 1 | compara_paralog | 0 | PF00855 |
| FLNB::FLNA | 368 | 269 | 22 | compara_paralog | 0 | PF00307, PF00630 |
| FLNA::FLNB | 356 | 231 | 23 | compara_paralog | 0 | PF00307, PF00630 |
| EEF1A1::EEF1A2 | 354 | 250 | 5 | compara_paralog | 0 | PF00009, PF22594 |
| YBX1::YBX3 | 236 | 190 | 3 | compara_paralog | 0 | PF00313 |
| HSP90AB1::HSP90AA1 | 233 | 232 | 2 | compara_paralog | 0 | PF00183, PF13589 |
| RHOA::RHOC | 193 | 193 | 1 | compara_paralog | 0 | PF00071 |

Top clean gene pairs, target_consolidated:

| Pair | Rows | Samples | Junctions | Relation | Read-through rows | Families |
|---|---|---|---|---|---|---|
| ACTG2::ACTG1 | 198 | 153 | 5 | compara_paralog | 0 | PF00022 |
| CMTM8::CMTM7 | 138 | 138 | 1 | compara_paralog | 0 | PF01284 |
| HBB::HBG2 | 71 | 29 | 8 | compara_paralog | 71 | PF00042 |
| HBB::HBG1 | 42 | 24 | 11 | compara_paralog | 42 | PF00042 |
| COL3A1::COL1A1 | 30 | 14 | 9 | compara_paralog | 0 | PF01391, PF01410 |
| SAA4::SAA1 | 30 | 30 | 1 | compara_paralog | 0 | PF00277 |
| CPSF7::CPSF6 | 27 | 27 | 1 | compara_paralog | 0 | PF00076 |
| FTL::FTH1 | 24 | 6 | 8 | compara_paralog | 0 | PF00210 |
| HBG1::HBA2 | 24 | 6 | 10 | compara_paralog | 0 | PF00042 |
| ACTG1::ACTA2 | 22 | 20 | 3 | compara_paralog | 0 | PF00022 |

Top clean gene pairs, target_starfusion:

| Pair | Rows | Samples | Junctions | Relation | Read-through rows | Families |
|---|---|---|---|---|---|---|
| ACTG2::ACTG1 | 153 | 153 | 2 | compara_paralog | 0 | PF00022 |
| CMTM8::CMTM7 | 138 | 138 | 1 | compara_paralog | 0 | PF01284 |
| HBB::HBG2 | 36 | 29 | 2 | compara_paralog | 36 | PF00042 |
| SAA4::SAA1 | 30 | 30 | 1 | compara_paralog | 0 | PF00277 |
| HBB::HBG1 | 28 | 24 | 2 | compara_paralog | 28 | PF00042 |
| CPSF7::CPSF6 | 27 | 27 | 1 | compara_paralog | 0 | PF00076 |
| ACTG1::ACTA2 | 21 | 20 | 2 | compara_paralog | 0 | PF00022 |
| SLC25A4::SLC25A6 | 10 | 10 | 1 | compara_paralog | 0 | PF00153 |
| TAF15::FUS | 9 | 9 | 1 | compara_paralog | 0 | PF00076 |
| EEF1A2::EEF1A1 | 8 | 8 | 1 | compara_paralog | 0 | PF00009 |

Top clean gene pairs, tcga_gtex:

| Pair | Rows | Samples | Junctions | Relation | Read-through rows | Families |
|---|---|---|---|---|---|---|
| ACTG2::ACTG1 | 2643 | 2589 | 5 | compara_paralog | 0 | PF00022 |
| SAA4::SAA1 | 1786 | 1786 | 1 | compara_paralog | 0 | PF00277 |
| CMTM8::CMTM7 | 1133 | 1133 | 1 | compara_paralog | 0 | PF01284 |
| ACTG1::ACTA2 | 1080 | 926 | 4 | compara_paralog | 0 | PF00022 |
| FLNA::FLNB | 731 | 707 | 14 | compara_paralog | 0 | PF00307, PF00630 |
| PRSS27::TPSAB1 | 671 | 671 | 1 | not_paralog_same_family | 0 | PF00089 |
| ACTA2::ACTG1 | 637 | 605 | 5 | compara_paralog | 0 | PF00022 |
| ACTG1::ACTG2 | 529 | 460 | 5 | compara_paralog | 0 | PF00022 |
| CELA2A::CTRC | 401 | 203 | 6 | not_paralog_same_family | 401 | PF00089 |
| KRT2::KRT1 | 394 | 347 | 5 | compara_paralog | 394 | PF00038 |
