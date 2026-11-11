# Domain-completing fusion pilot: open tables (GRCh38) and MSK-IMPACT 50k (GRCh37)

**Question.** Are there fusions where a protein domain is truncated in both partners and the in-frame junction rebuilds a complete, chimeric domain? (Paralogous ETS fusions are the model case.)

**Short answer.**

- No convincing domain-completing fusion was found in any of the five cohorts.
- MSK-IMPACT 50k reproduces #153 exactly: 0 clean, 0 ambiguous.
- The RNA open tables give tens of thousands of `clean_reconstitution` rows: 35,716 rows, 2,330 distinct junctions. These are overwhelmingly a known RNA-seq artefact.
  - Highly expressed paralogue pairs (actins, calmodulins, EEF1A, SAA, CMTM, haemoglobins, collagens, keratins, filamins) are joined at their homologous position, at annotated splice sites.
  - 92–99.5% of clean rows are Ensembl Compara paralogue pairs.
  - The paralog-sequence test flags 69–90% of them.
- The classifier applies its definition correctly. That definition, however, is not specific for this question, and it over-calls in two ways (see [Is "clean" over-calling?](#is-clean-over-calling)).
- After the three artefact flags, a handful of low-support calls remain: 1–8 junction reads, FFPM ≤ 1. None involves an oncogenic family.

## Files

| File | Contents |
|---|---|
| `<cohort>.domain_screen.tsv.gz` | Every shared-family row of the cohort (all classes except `no_shared_family`) with the three artefact flags and their sub-checks. Rows are sorted by class, gene pair, breakpoints, family and sample. |
| `combined_candidates.tsv.gz` | All five cohorts. Only rows of class clean, ambiguous, tandem, gap, out-of-frame, junction-stop and splice-site-unresolved (56,216 rows). |
| `flag_summary.json` | Screen summaries and flag counts per cohort. |
| `clean_call_audit.{json,md}` | The clean-call breakdown below (`scripts/domain_screen_pilot_clean_audit.py`). |
| `tier2_candidates.json.gz` | 609 distinct junctions that pass every check except sample burden (`--tier2`). |
| `tier2_exclusions.json` | Rows left out of tier 2 because a breakpoint-window measurement is missing or malformed, by reason (empty in this run). |
| `report_tables.md` | Generated: every non-trivial call grouped by distinct junction (with every sample and its reads), the FLI1::ERG rows and full MSK flag reasons (`scripts/domain_screen_pilot_report.py`). |
| `run_manifest.json` | Commit SHAs, builds, service URLs, cache dirs, SHA-256 of every input, exact command lines, attempt timings. |

Raw downloads, normalized tables and caches are not committed; they live in the git-ignored `runs/.domain_screen_pilot_work/`.

## Class counts per cohort (shared-family rows)

| Class | MSK-IMPACT 50k (GRCh37) | TARGET STAR-Fusion | TARGET consolidated | DepMap 22Q2 | TCGA+GTEx |
|---|---|---|---|---|---|
| clean_reconstitution | 0 | 602 | 819 | 9,567 | 24,728 |
| ambiguous_junction_register | 0 | 65 | 95 | 1,802 | 9,449 |
| tandem_duplication_or_insertion | 2 | 94 | 211 | 317 | 1,286 |
| gap_or_deletion | 0 | 21 | 116 | 243 | 809 |
| out_of_frame | 2 | 60 | 103 | 585 | 4,010 |
| premature_stop_at_junction | 0 | 0 | 0 | 0 | 4 |
| splice_site_unresolved | 3 | 0 | 26 | 1,177 | 20 |
| not_split_on_both | 66 | 2,008 | 6,436 | 20,904 | 46,911 |
| no_transcript_model | 17 | 6,004 | 6,543 | 6,715 | 47,751 |
| no_shared_family (events) | 5,148 | 226,514 | 247,194 | 331,518 | 642,375 |
| *screened fusion events* | 5,201 | 231,188 | 254,143 | 356,103 | 716,420 |
| *events sharing a family* | 53 | 4,396 | 4,396 | 16,897 | 64,118 |
| *shared-family rows* | 90 | 8,854 | 14,349 | 41,310 | 134,968 |

The 1,177 DepMap `splice_site_unresolved` rows (an intronic break inside a splice site) probably come from the non-reference-splice calls that DepMap's unfiltered table keeps. They were not investigated further.

## MSK-IMPACT 50k vs #153

| Class | #153 | This run | Match |
|---|---|---|---|
| clean_reconstitution | 0 | 0 | yes |
| ambiguous_junction_register | 0 | 0 | yes |
| tandem_duplication_or_insertion | 2 | 2 (FLI1::ERG ETS `YYYDKNI`; STAT5B::STAT5A PF02864 `KQPPQVLKTQTKFAAT`) | yes |
| out_of_frame | 2 | 2 (HCK::LYN, BIRC2::BIRC3) | yes |
| splice_site_unresolved | 3 | 3 | yes |
| not_split_on_both | 66 | 66 | yes |
| no_transcript_model | 17 | 17 | yes |
| no_shared_family | 5,148 | 5,148 | yes |
| family rows / events sharing a family | 90 / 53 | 90 / 53 | yes |

The MSK screen ran on a new, empty GRCh37 cache (www.genomenexus.org, grch37.rest.ensembl.org). Its output at the final screen commit is byte-identical to the first run of this pilot (`run_manifest.json`).

MSK non-trivial calls with flags (full reasons in `report_tables.md`):

| Fusion | Sample | Family | Class | Junction peptide | Reads (split/pairs) | Frame (translated / alt) | Paralog | Genomic | Normal |
|---|---|---|---|---|---|---|---|---|---|
| FLI1::ERG | P-0053737-T01-IM6 | PF00178 ETS | tandem (7 aa `YYYDKNI`) | `YDKLSRALRYYYDKNI\|YYYDKNIMTKVHGKR` | 0/6 IMPRECISE | in / out | pass | **FLAG**: breakpoint windows 83.3% over 275 bp; 3 shared-family calls in sample; 0 split reads, imprecise | pass |
| STAT5B::STAT5A | P-0029972-T01-IM6 | PF02864 | tandem (16 aa) | `KQPPQVLKTQTKFAAT\|KQPPQVLKTQTKFAA` | 0/5 IMPRECISE | in / out | pass | **FLAG**: windows 95.3% over 253 bp; both breakpoints in a segdup (fracMatch 0.974); 0 split reads | pass |
| HCK::LYN | P-0053737-T01-IM6 | PF07714 | out_of_frame | `NKHTKVAVKTMKPGS\|ILPGRSQPHEDPAA` | 0/6 IMPRECISE | out / in | pass | **FLAG**: windows 80.9% over 183 bp; burden 3; 0 split reads | pass |
| BIRC2::BIRC3 | P-0013189-T02-IM6 | PF00653 | out_of_frame | `ESGDDPWVEHAKWFP\|R` | 10/12 PRECISE | out / out | pass | **FLAG**: 5' in L1MB3 (LINE/L1) | pass |
| KMT2B::KMT2E | P-0037777-T01-IM6 | PF00628, PF00856 | splice_site_unresolved (2 rows) | - | 0/8 IMPRECISE | unknown | n/a (no CDS) | **FLAG**: 3' in FLAM_C (Alu); 0 split reads | pass |
| ADAMTSL1::PTPRD | P-0063201-T01-IM7 | PF07679 | splice_site_unresolved | - | 10/11 PRECISE | unknown | n/a (no CDS) | pass | pass |

The MSK FLI1::ERG sample (P-0053737) carries 3 distinct shared-family calls: EP300::CREBBP (0 split / 40 pairs; `not_split_on_both` for all seven of its families), FLI1::ERG (0/6) and HCK::LYN (0/6). All three are paired-end-only and imprecise. Two of them have similar breakpoint regions: FLI1::ERG at 83.3% over 275 bp (±150 bp) and HCK::LYN at 80.9% over 183 bp. EP300::CREBBP's do not pass the threshold: 86.9% over only 61 bp (±150 bp) and 76.0% over 612 bp (±500 bp). This is the sample-level artefact signature the brief described.

## The TARGET FLI1::ERG row

Both TARGET tables carry the same call, which is FusionInspector-validated in the consolidated table:

| Field | Value |
|---|---|
| Sample | TARGET-10-PAPSPG-04A (B-ALL, tumour) |
| Breakpoints (GRCh38) | chr11:128782023 (FLI1 exon 5, c.655) → chr21:38445621 (ERG exon 2, c.19) |
| Transcripts | NM_002017 / NM_182918 |
| Reads | 2 junction / 0 spanning, FFPM 0.0206, FI-validated TRUE (consolidated) |
| Junction peptide | `TTSHTDQSSRLSVKE\|GSLISCE` |
| Frame | **out of frame** under the STAR-Fusion (translated) convention; in frame only under the off-by-one alternative |
| PF00178 ETS | `not_split_on_both` (FLI1's ETS is lost, being 3' of the break; ERG's is intact) |
| PF02198 SAM/PNT | `not_split_on_both` (both intact) |
| paralog_sequence_flag | pass: best single sequence ERG (ENSG00000157554) 100% over 53/100 nt; below the 90% coverage threshold |
| genomic_similarity_flag | **flagged**: sample has 3 shared-family calls. Breakpoint windows are not similar (+/-150: 72% over 25 bp). No segdup/repeat. |
| normal_recurrence_flag | pass: junction absent from TCGA/GTEx normals; gene pair absent from TCGA+GTEx altogether |

So the TARGET FLI1::ERG call joins FLI1's N-terminus to almost all of ERG at exon boundaries. It does not split the ETS domain, and it is frameshifted under the caller's own convention. It is not a domain-completing fusion. The MSK FLI1::ERG rebuilds the ETS domain only through a DNA break inside the ETS exons. An RNA caller anchored on splice sites cannot report that kind of junction (see [Limits](#limits)).

## Is "clean" over-calling?

A clean call requires all of the following:

- a copy of the family is truncated in both partners (`junction_domain_classifier.py:391-404`);
- in frame;
- no junction stop;
- the chimeric domain has exactly one parent's length with no indel, under every co-optimal alignment.

So it really does need a domain cut in both partners. By construction, though, every clean call is a stitch at the homologous position of two homologous domains. That is exactly what template switching or multi-mapping between paralogous transcripts produces, as does read-through between neighbouring family members.

The breakdown below (`clean_call_audit.md`, Ensembl 116 Compara paralogues; read-through = same chromosome, same transcribed strand, ≤200 kb):

| Measure | TARGET STAR-Fusion | TARGET consolidated | DepMap 22Q2 | TCGA+GTEx |
|---|---|---|---|---|
| clean rows | 602 | 819 | 9,567 | 24,728 |
| distinct junctions (pair + breakpoints + family) | 96 | 169 | 831 | 1,234 |
| distinct samples | 404 | 415 | 1,282 | 10,874 |
| distinct gene pairs | 79 | 89 | 541 | 840 |
| same gene (5' == 3') | 0 | 0 | 0 | 0 |
| read-through (same chrom + strand, ≤200 kb) | 69 | 119 | 227 | 3,885 |
| any different-gene pair on one chromosome ≤200 kb | 99 | 150 | 302 | 6,171 |
| Compara human paralogues | 598 (99.3%) | 815 (99.5%) | 9,484 (99.1%) | 22,707 (91.8%) |
| same Pfam family, not Compara paralogues | 4 | 4 | 83 | 2,021 |
| both ends at annotated splice sites (ONLY_REF_SPLICE) | 600 | n/r* | 9,437 | 24,610 |
| one side contributes < 10 aa or < 10% of its domain | 62 | 108 | 1,722 | 3,428 |
| both sides contribute substantially | 540 | 711 | 7,845 | 21,300 |
| chimera ≥ 95% identical to one parent domain | 286 | 393 | 7,095 | 11,577 |
| paralog_sequence_flag = flagged | 428 | 562 | 8,646 | 17,877 |
| all three artefact flags pass | 0 | 0 | 0 | 11 |
| strict: not read-through, both sides substantial, all flags pass | 0 | 0 | 0 | 7 rows / 5 junctions |

\* The consolidated TARGET table has no splice-type column. MSK has no clean rows.

The top clean pairs are the same in every RNA cohort. Rows / samples:

- **TCGA+GTEx:** ACTG2::ACTG1 2,643 / 2,589; SAA4::SAA1 1,786; CMTM8::CMTM7 1,133; ACTG1::ACTA2 1,080; FLNA::FLNB 731; PRSS27::TPSAB1 671; KRT2::KRT1 394 (all read-through).
- **DepMap:** ACTG2::ACTG1 827 / 818; CALM3::CALM2 657; EEF1A2::EEF1A1 584; HDGFL3::HDGF 435; FLNB::FLNA 368; HSP90AB1::HSP90AA1 233.
- **TARGET:** ACTG2::ACTG1 153; CMTM8::CMTM7 138; HBB::HBG2 / HBG1 (read-through).

ACTG2::ACTG1 joins c.987 to c.985 of two actins, and the chimera is 99.5% identical to ACTG2. Kinase "clean" calls follow the same pattern. ERBB3::EGFR at `YTHQSDVWSYG|VTVWELMTFG` (76 TCGA+GTEx rows, 28 DepMap), YES1::FYN, EPHB3::EPHB1 and MAP4K4::MINK1 are paralogue template switches at conserved kinase motifs, and nearly all are flagged.

**Verdict.**

1. **Over-calls by construction, not by bug.** The classifier does what its definition says, but on RNA fusion tables "clean" is not specific for a real domain-completing fusion. It cannot tell a rebuilt paralogous domain from mis-mapping between paralogues. 69–90% of clean rows are flagged by the paralog-sequence test, and 47–74% have a chimera ≥95% identical to a parent domain.
2. **Over-calls through weak truncation.** "Truncated" means any cut inside a domain span, however small. 10–18% of clean rows take < 10 residues or < 10% of the domain from one side, so the domain is in effect intact in the other partner. Example: PTN::MDK keeps 1 residue of PTN's domain (`96-158 kept 96-96`).

**Follow-ups (not in this PR):**

- Require a minimum contribution from each side.
- Report or gate the chimera's identity to each parent.
- Exclude or flag adjacent same-strand pairs inside the classifier.
- Treat known paralogue-artefact pairs (STAR-Fusion's own paralogue lists) as a separate class.

## Artefact flags (definitions and thresholds)

Each flag is one of:

- `flagged`;
- `pass`;
- `not_assessed`: should have run but could not, with the reason;
- `not_applicable`: a `no_transcript_model` row that never reached junction translation.

None passes silently, and every flag has a reason column in the TSVs.

- **`paralog_sequence_flag`.** Takes 50 nt of the fused CDS either side of the predicted junction and aligns it locally to each single candidate source: the 5' and 3' parent CDS, and the CDS of every Ensembl 116 Compara paralogue of either partner. **Flagged** when one sequence covers ≥90% of the window at ≥90% identity, i.e. one unfused transcript would explain the "junction". `not_assessed` when there is no fused CDS (`splice_site_unresolved`), or when the comparison set is incomplete (a missing parent CDS, a failed or unresolved paralogue lookup, or a Compara member without a CDS), unless an available sequence already explains the junction (then `flagged`). A verified "no paralogues" is complete and noted in the reason. Before full traceback, targets are screened by local-alignment score with a cutoff derived from the window length (`explaining_min_score`), so windows clipped at the CDS ends (< 100 nt) are tested as strictly as full ones.
- **`genomic_similarity_flag`.** Flagged when any of these sub-checks fires (each has its own column):
  - **Breakpoint similarity.** Reference ±150 bp and ±500 bp around the two breakpoints, aligned locally in both orientations. Flagged at **≥80% identity over ≥100 aligned bp**. Percent identity and aligned length are reported for both windows.
  - **`breakpoint_repeat_flag`.** The breakpoint base lies in a UCSC `genomicSuperDups` or `rmsk` element of the row's own build (hg19 for MSK, hg38 for the open tables). Both tracks were downloaded and used, so this check is assessed for every translated row.
  - **`sample_burden_flag`.** The sample has ≥2 distinct shared-family calls in its cohort.
  - **`read_evidence_flag`.** 0 split/junction reads, or an `IMPRECISE` breakpoint (MSK).
- **`normal_recurrence_flag`.** Flagged when the exact junction (gene pair and both breakpoint coordinates on GRCh38) appears in another TCGA or GTEx **normal** sample of the TCGA+GTEx table. MSK GRCh37 breakpoints are lifted base by base with Ensembl `map`. A TCGA+GTEx call that is itself from a normal sample is also flagged. GTEx cell-line samples (EBV-transformed lymphocytes, cultured fibroblasts) are classed `cell_line` by `normalize_table`. They are counted separately (`cell_line_exact_junction_samples`) and never flag. The reason also gives the gene pair's counts by tissue class.

Flag counts per cohort, over shared-family rows:

| Cohort | paralog_sequence | genomic_similarity | breakpoint_repeat | sample_burden | read_evidence | normal_recurrence |
|---|---|---|---|---|---|---|
| MSK 50k | 0 flagged / 70 pass / 3 n/a | 57 / 16 | 36 / 37 | 13 / 60 | 29 / 44 | 0 / 73 |
| TARGET STAR-Fusion | 826 / 2,024 | 2,750 / 100 | 418 / 2,432 | 2,722 / 128 | 0 / 2,850 | 1,414 / 1,436 |
| TARGET consolidated | 1,139 / 6,641 / 26 n/a | 7,717 / 89 | 623 / 7,183 | 7,689 / 117 | 0 / 7,806 | 2,332 / 5,474 |
| DepMap 22Q2 | 20,061 / 13,357 / 1,177 n/a | 34,567 / 28 | 6,208 / 28,387 | 34,558 / 37 | 0 / 34,595 | 19,913 / 14,682 |
| TCGA+GTEx | 38,412 / 48,785 / 20 n/a | 85,281 / 1,936 | 21,449 / 65,768 | 84,501 / 2,716 | 0 / 87,217 | 78,792 / 8,425 |

(flagged / pass [/ not_assessed]; `no_transcript_model` rows are `not_applicable`.)

In the RNA tables the sample-burden sub-check fires on almost every row: one RNA sample typically carries many shared-family paralogue calls. Requiring all three flags to pass is therefore dominated by burden < 2, so the candidates below come in two tiers.

## Strongest candidates after the artefact flags

### Tier 1: clean / ambiguous / indel calls passing all three flags

25 rows; 18 distinct junctions (the two TARGET junctions appear in both TARGET tables). Full flag reasons are in `report_tables.md`.

| Fusion | Cohort | Sample (tumour type) | Family | Class | Junction peptide | Reads split/pairs, FFPM | Frame (translated / alt) | Note |
|---|---|---|---|---|---|---|---|---|
| GPC6::GPC5 | TCGA | LIHC-CC-A3M9-TP | PF01153 | gap_or_deletion | `AKGFSLADIPYQEIAG\|MPDDMNFSDVKQIHQ` | 8/34, 0.60 | in / out | same-strand neighbours 361 kb apart on chr13; 32 aa of GPC6 + 47 aa of GPC5, ~470 aa deleted |
| GPC4::GPC3 | TCGA | TGCT-2G-AAKD-TP | PF01153 | gap_or_deletion | `SKGFNKNDAPLHEINE\|LAYDLDVDDAPGNSQ` | 8/25, 0.97 | in / out | **neighbouring genes** (same strand, 121 kb apart on Xq26, inside the 200 kb cut-off, so read-through by this report's definition; the 5'→3' order runs against GPC3/GPC4 transcription); 41 aa of GPC4 + 52 aa of GPC3, ~470 aa deleted |
| GPC3::GPC4 | TCGA | LIHC-2Y-A9GU-TP | PF01153 | tandem | `GDGMIKVKNQLRFLAG\|DHLKICPQGSTCCSQ` | 5/23, 0.52 | in / out | same pair in the other order, same strand, 257 kb apart (outside the 200 kb cut-off) |
| LPO::MPO | TCGA | 3 LAML (AB-2862, AB-2982, AB-3012) | PF03098 | clean | `CEFINTTARVPCFLAG\|DTRSSEMPELTSMHT` | 1–2/4–6, 0.07–0.11 | in / out | neighbours 23 kb apart on chr17, opposite strands; MPO is very highly expressed in AML |
| PIM3::ULK1 | TARGET (both tables) | TARGET-10-PARWRJ-09B | PF00069 | tandem (76 aa insertion) | `GSGALLKDTVYTDFDA\|MRTLSEDTIRLFLQQ` | 2/0, 0.019, FI TRUE | in / out | kinase domains of unrelated kinases |
| COL1A1::COL18A1 | TARGET (both tables) | TARGET-40-PAMRHD-01A | PF01391 | tandem | `GPPGPRGRTGDAGPV\|GPPGPPGREGPPGRT` | 2/1, 0.034, FI FALSE | in / out | collagen repeats |
| GNAL::GNAS | TCGA | THCA-KS-A41I-TP | PF00503 | clean | `PKVTRAKFFIRDLFL\|RISTASGDGRHYCYP` | 1/2, 0.045 | in / out | chimera 99.4% identical to GNAL |
| KRT82::KRT7 | TCGA | BLCA-CF-A5UA-TP | PF00038 | clean | `SLRPCVENEFVALKK\|DVDAAYMSKVELEAK` | 1/2, 0.050 | in / out | type II keratin cluster |
| VIM::KRT18 | TCGA | LUAD-05-4415-TP | PF00038 | clean | `SLQEEIAFLKKLHEE\|EVKGLQAQIASSGLT` | 1/5, 0.10 | in / out | intermediate filaments |
| YWHAB::YWHAZ | TCGA | LUSC-18-3419-TP | PF00244 | clean | `EKIEAELQDICNDVL\|SLLEKFLIPNASQAE` | 1/2, 0.055 | in / out | 14-3-3 paralogues |
| COL1A1::COL12A1 | TCGA | 3 (BRCA, BLCA, COAD) | PF01391 | gap_or_deletion | `YDEKSTGGISVPGPM\|GPPGPRGDIGPPGPQ` | 1/1 each | in / out | 3 aa of COL1A1's domain |
| AFM::ALB | TCGA | LIHC-UB-A7MB-TP | PF00273 | clean | `DGCCEGDVVQCIRDT\|ADLAKYICENQDSIS` | 2/0, 0.031 | in / out | read-through (same strand, 75 kb) |
| CRIP1::CRIP2 | TCGA | BRCA-E2-A14W-TP | PF00412 | clean | `CEKCGKTLTSGGHAE\|HDGKPFCHKPCYATL` | 1/1, 0.036 | in / out | read-through (same strand, 10 kb) |
| PTN::MDK | TCGA | BRCA-AC-A62X-TP | PF01091 | clean | `TQRCKIPCNWKKQFGA\|DCKYKFENWGACDGG` | 1/1, 0.032 | in / out | 1 residue from PTN |
| KLK3::PRSS8 | TCGA | PRAD-G9-6354-TP | PF00089 | ambiguous | `MLCAGRWTGGKSTCS\|GDSGGPLSCPVEGLW` | 1/1, 0.035 | in / out | KLK3 is highly expressed in prostate |
| CMTM7::CMTM8 | TCGA | BLCA-XF-A9SK-TP | PF01284 | gap_or_deletion | `SYPRTHAALLKVAQM\|GLCFNGSAFVLYLSA` | 1/6, 0.13 | in / out | |
| COL3A1::COL1A1 | TCGA | BLCA-XF-A9SZ-TP | PF01391 | gap_or_deletion | `GQQGAIGSPGPAGPR\|GPQGPRGDKGETGEQ` | 1/5, 0.13 | in / out | |
| SPNS2::SPNS3 | TCGA | BLCA-DK-A1A3-TP | PF07690 | tandem | `TLLFSNWAITADILM\|VMPCLEAVALILLIL` | 1/4, 0.097 | in / out | |

All of these are in frame under the STAR-Fusion convention, and all three flags pass for every one. Every candidate has thin support, and nearly all are adjacent or co-expressed paralogues. The best-supported are the glypican pairs (8 junction reads; FFPM 0.5–1.0, in liver cancer and a germ-cell tumour). They are same-strand tandem paralogues: GPC4::GPC3 is 121 kb apart, inside the 200 kb cut-off, so it counts as neighbouring-gene read-through by this report's definition. GPC3::GPC4 (257 kb) and GPC6::GPC5 (361 kb) fall just outside it. The two gap calls keep only ~40–50 aa from each side. A read-through or local rearrangement is more likely than a domain-completing fusion. Validation would need the reads.

### Tier 2: every check except sample burden

609 distinct junctions (`tier2_candidates.json.gz`). Criteria:

- clean, ambiguous, tandem or gap;
- in frame (translated convention);
- paralog-sequence pass and normal-recurrence pass;
- breakpoint-window similarity measured on both windows (numeric identity and aligned length) and not high; no segdup/repeat; split reads present (rows without both window measurements are excluded and counted in `tier2_exclusions.json`: 0 in this run);
- not read-through;
- both sides contribute ≥10 aa and ≥10% of the domain.

Of these, 554 are Compara paralogues. The top 15 by number of samples:

| # | Fusion | Cohort | Class | Family | Samples | Max split reads | Junction peptide | Frame | Chimera id. to parent |
|---|---|---|---|---|---|---|---|---|---|
| 1 | COL7A1::COL1A1 | TCGA+GTEx | ambiguous | PF01391 | 63 | 4 | `GDKGPRGDNGDPGDK\|GESGPSGPAGPTGAR` | in / out | 1.0 |
| 2 | KRT8::KRT5 | DepMap | clean | PF00038 | 38 | 8 | `NKRTEMENEFVLIKK\|DVDAAYMNKVELEAK` | in / out | 0.83 |
| 3 | FLNA::FLNB | DepMap | clean | PF00630 | 36 | 2 | `RKDGSCGVAYVVQEPG\|NYEVSIKFNDEHIPE` | in / out | 0.89 |
| 4 | PTMA::PTMS | DepMap | gap | PF03247 | 33 | 3 | `EAESATGKRAAEDDE\|DEADPKRQKTENGAS` | in / out | 0.96 |
| 5 | GFAP::VIM | TCGA+GTEx | clean | PF00038 | 18 | 3 | `ASSNMHEAEEWYRSK\|FADLSEAANRNNDAL` | in / out | 0.89 |
| 6 | EIF4A1::EIF4A2 | DepMap | clean | PF00270 | 17 | 2 | `TQALVLAPTRELAQQ\|IQKVILALGDYMGAT` | in / out | 0.96 |
| 7 | FLNB::FLNA | DepMap | clean | PF00630 | 15 | 2 | `DVMVFVEDPEGNKEE\|AKVTANNDKNRTFSV` | in / out | 0.88 |
| 8 | KLK4::KLK3 | TCGA+GTEx | clean | PF00089 | 14 | 1 | `MFCAGGGQDQKDSCN\|GDSGGPLVCNGVLQG` | in / out | 0.90 |
| 9 | COL12A1::COL1A1 | DepMap | tandem | PF01391 | 11 | 8 | `GERGISGAIGPPGPRG\|HRGFSGLQGPPGPPG` | in / out | 0.89 |
| 10 | COL1A1::COL5A1 | DepMap | gap | PF01391 | 11 | 2 | `GKPGRPGERGPPGPQ\|GPRGVQGPPGPAGKP` | in / out | 0.91 |
| 11 | FLNB::FLNA | DepMap | clean | PF00630 | 11 | 1 | `RTFEMSDFIVDTRDAG\|YGGLSLSIEGPSKVD` | in / out | 0.87 |
| 12 | FLNA::FLNB | DepMap | clean | PF00630 | 10 | 2 | `DAKDAGEGLLAVQIT\|DQEGKPKRAIVHDNK` | in / out | 0.81 |
| 13 | COL3A1::COL1A1 | TCGA+GTEx | tandem | PF01391 | 10 | 1 | `GMPGPRGSPGPQGVK\|GESGPSGPAGPTGAR` | in / out | 0.98 |
| 14 | KLK3::PRSS8 | TCGA+GTEx | ambiguous | PF00089 | 10 | 1 | `MLCAGRWTGGKSTCS\|GDSGGPLSCPVEGLW` | in / out | 0.92 |
| 15 | COL1A1::COL3A1 | DepMap | gap | PF01391 | 9 | 13 | `GLPGERGRPGAPGPA\|GPKGNDGAPGKNGER` | in / out | 0.89 |

Tier 2 is the same artefact class with the paralogues diverged below 90% nucleotide identity. It is dominated by fibroblast collagens, epithelial keratins, filamins and prostate KLKs: highly expressed, co-expressed paralogues that recur across samples.

## Review re-run (832a267)

The #154 review found four ways a check could pass on missing data. Each now has a regression test that fails on 854bb1d and passes at 832a267:

1. A missing parent CDS, or a failed paralogue lookup, let `paralog_sequence_flag` pass.
2. A missing split-read count did not propagate to `genomic_similarity_flag`.
3. The traceback cutoff was fixed for 100 nt windows.
4. Tier 2 read blank window measurements as zero.

Only the flag pass, tables, audit and packaging were re-run. The screen outputs are unchanged.

| Measure (shared-family rows) | MSK 50k | TARGET SF | TARGET cons. | DepMap | TCGA+GTEx |
|---|---|---|---|---|---|
| any flag value changed | 0 | 0 | 0 | 0 | 0 |
| rows whose paralogue comparison is incomplete (now `not_assessed`) | 0 | 0 | 0 | 0 | 0 |
| `read_evidence_flag` not_assessed | 0 | 0 | 0 | 0 | 0 |
| best paralogue target changed (reason or target only) | 2 | 341 | 376 | 275 | 761 |
| best-target identity changed (clipped windows < 100 nt) | 0 | 1 | 1 | 37 | 56 |
| all three flags pass, before → after | 15 → 15 | 54 → 54 | 50 → 50 | 0 → 0 | 86 → 86 |
| clean rows passing all three flags, before → after | 0 → 0 | 0 → 0 | 0 → 0 | 0 → 0 | 11 → 11 |

Other measures:

- **Tier 1:** 25 → 25 rows; 18 → 18 junctions.
- **Tier 2:** 609 → 609 junctions; 0 rows excluded for missing window measurements.
- **Top-10 list:** unchanged.

The review's bugs were real but latent in these data:

- Every translated row's parent CDS and every Compara member's CDS is in the Ensembl 116 release files.
- All MSK GRCh37 homology lookups, re-fetched under a new cache key that records completeness, succeeded.
- Every row has a split-read count and both window measurements.

The cutoff fix changed the best single-sequence target on short windows. One example is PTMA::PTMS: before, the PTMA parent at 91.0% over 82% of the window; after, the PTMS parent at 82.3% over 100%. None of these changes reaches the 90%/90% "explains" threshold.


### Review round 2 (fail-closed parsing; outputs unchanged)

Two more fail-closed rules followed:

1. **REST paralogue lookups.** A REST paralogue lookup is now complete only when the response resolves the gene: a non-empty `data` list whose entries carry an `id` and an explicit `homologies` list, with an id and a sequence for every member. Responses `{}`, `{"data": []}`, `{"data": [{}]}`, an entry without `homologies`, or an unresolved gene are all incomplete, so the paralog flag is `not_assessed`. A resolved gene with `homologies: []` is the only verified empty. The cache key is now `paralog_lookup_v2_<gene>`, so entries stored under the earlier rules are re-fetched.
2. **Tier-2 window measurements.** Tier 2 now also rejects non-finite or out-of-range identities and negative or fractional aligned lengths. Such rows are `unassessed` and counted in `tier2_exclusions.json`.

Neither rule changes these outputs, so the flag pass was not re-run:

- **REST lookups.** Only the MSK GRCh37 rows use REST paralogue lookups; the GRCh38 cohorts read the Ensembl 116 release files. All 58 MSK lookups were re-fetched under the new rule. All 58 are complete and identical, in sequences and completeness, to what the committed flags used.
- **Window measurements.** Across all five committed TSVs, all 530,164 non-blank window measurements are valid (finite, identity in [0, 1], integral non-negative length). 0 are non-finite, out of range, fractional or non-numeric. Blanks occur only outside the tier-2 classes.
- **Tier 2.** Re-run under the new rule, it reproduces `tier2_candidates.json` and `tier2_exclusions.json` byte for byte: 609 junctions, 0 exclusions.

## How the pilot was run

- **Inputs.**
  - The four GRCh38 open tables were fetched with `fetch_table`, which writes each download to a new directory and checks the SHA-256 pinned in `open_fusion_manifest.py`. They were then normalized with `normalize_table`.
  - MSK-IMPACT 50k (GRCh37) is the archive used in #152/#153 (SHA-256 `c9522b3e…9e7d38`).
  - Each download went to its own new directory, and every Python process that read one ran with `-I`.
- **Screen.**
  - `cfh domain-screen --mode any-domain`.
  - Open tables use `--input-format open-fusion-table` (`src/cfh/ingestion/open_fusion_sv_adapter.py`). This maps each normalized row to a structural-variant row unfiltered. 5'→3' order is explicit, positions and build are copied verbatim, junction reads become split reads and spanning fragments become paired-end reads.
  - `--omit-unshared-rows` counts `no_shared_family` events without writing them.
  - GRCh38 uses grch38.genomenexus.org and rest.ensembl.org; GRCh37 uses www.genomenexus.org and grch37.rest.ensembl.org.
- **Caches.**
  - Each build has its own cache directory, created empty for this pilot and stamped with its build. One GRCh37 directory serves MSK; one GRCh38 directory is shared by the four open tables, so genes are not re-fetched.
  - Failed attempts were re-run unchanged on the same cache. The failures were Genome Nexus/Ensembl timeouts and a `/`-containing DepMap symbol, fixed in 13002e3; `run_manifest.json` lists every attempt.
  - The GRCh38 CDS cache was seeded from Ensembl release 116 `cds.all.fa.gz`, the release REST was serving. All 1,042 REST-fetched entries were identical to the file.
- **Determinism.** Re-running MSK and both TARGET tables at the final screen commit gave byte-identical TSVs and summaries.
- **Commits.**
  - Screens: a8a33b1 (MSK, TARGET ×2, TCGA+GTEx) and 13002e3 (DepMap). The only change between them percent-encodes gene symbols in single Genome Nexus lookups, which affects the DepMap symbol `THRA1/BTR`.
  - Flags: 832a267 (review fixes, after #154 review; the 13002e3 flag outputs are superseded, see [Review re-run](#review-re-run-832a267)).
- **Flag data.**
  - Reference windows come from UCSC `hg38.2bit` / `hg19.2bit`, matched base-for-base against 112 Ensembl REST windows. The Ensembl sequence endpoint was taking 15–55 s per 1 kb region during the pilot.
  - Paralogue and parent CDS come from Ensembl 116 release files (`cds.all`, `pep.all`, Compara `protein_default` homologies; BSD checksums match Ensembl's `CHECKSUMS`). MSK GRCh37 rows use grch37.rest.ensembl.org.
  - Repeat tracks are UCSC `genomicSuperDups` and `rmsk` for hg19 and hg38.

### Coordinate conventions and frame

- `translated_frame_status` uses the screen's default convention, `both_bases_retained`: the reported 5' base is the last base kept, and the reported 3' base is the first base kept.
- `alternative_convention_frame_status` drops the reported 3' base.
- **Open tables.** STAR-Fusion's `LeftBreakpoint` / `RightBreakpoint` are the last donor-exon base and the first acceptor-exon base, i.e. `both_bases_retained`. The data agree. On distinct shared-family junctions with a frame, "in frame (translated) / out (alt)" outnumbers "out / in" 482:121 in TARGET STAR-Fusion, 4,369:344 in TCGA+GTEx and 3,066:163 in DepMap. **For the open tables, read the translated frame as the frame**; the alternative is an off-by-one control.
- **MSK (DNA).** The SV caller's convention is undocumented, so a call whose two frames disagree has an unresolved frame, as in #153.
- **Mechanism label.** The translator labels any join between two coding positions `mid_exon_dna`, exon edges included. For STAR-Fusion junctions at annotated splice sites, that join is exactly the spliced exon–exon product. The label describes how the product was assembled; it is not evidence of a DNA break.

## Limits

- **No reads.** The open tables have read counts but no read or junction sequences, so there is no read-level paralogue or multi-mapping test. The paralog-sequence flag tests only whether a ±50 nt window of the *predicted* fused CDS could come from one unfused transcript. It cannot see soft-clipped, multi-mapped or low-complexity reads, and it passes paralogues diverged below 90% nucleotide identity, which is why tier 2 is still mostly artefacts.
- **RNA calls sit on splice sites.** STAR-Fusion junctions are spliced exon–exon joins. A domain-completing fusion that needs a DNA break *inside* an exon cannot be seen in these tables: the MSK FLI1::ERG joins inside the ERG/FLI1 ETS-encoding exons and depends on that. In RNA it would appear only as an exon–exon join, which typically keeps or drops the ETS domain whole.
- **"Clean" is a candidate, not proof.** Pfam envelopes from the parent transcripts stand in for the family; there is no HMM. A clean call says the protein sequence is consistent with a rebuilt domain. It is not evidence of expression, folding or function. As shown above, on RNA tables it is dominated by paralogue artefacts.
- **Normal recurrence uses one table.** Only the TCGA+GTEx table is searched, by exact GRCh38 coordinates. The same junction under a different symbol is counted and named in the reason, but it does not flag. MSK calls are lifted from GRCh37 base by base; all 73 translated MSK rows lifted.
- **The thresholds are judgment calls.** Paralog sequence ≥90%/≥90%; breakpoint windows ≥80% over ≥100 bp; burden ≥2; read-through ≤200 kb, same strand; minor side < 10 aa or < 10%. All are recorded per row so they can be re-cut.
- **Caller conventions.** Frame for MSK depends on an undocumented breakpoint convention. The open tables follow STAR-Fusion's.
