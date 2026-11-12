# Sequence-level functional-domain checks: domain-screen pilot, 2026-11-11

PR 1 of 3 (sequence metrics). PR 2 adds structure prediction and PR 3 adds
stability and protein-language-model scores; both read `chimera_records.json.gz`.
Methods, every metric and every tier rule are documented in
[`docs/domain_function_seq.md`](../../docs/domain_function_seq.md).

**The tiers are HEURISTIC.** No calibrated intact/broken cutoff exists for these
metrics; calibrating one would need assayed chimeras. A tier summarises the
reported metrics with documented rules. It does not predict activity. There is
no combined score.

This run replaces the first one (commit 334c050) after review. The section
[Changes after review](#changes-after-review) gives the before/after numbers.

## Run

| | |
|---|---|
| input | `runs/domain_screen_pilot_2026-11-10/combined_candidates.tsv.gz` (sha256 in `run_manifest.json`) |
| code | commit `00680e0` (`code_dirty: false`). Records were rebuilt from scratch at `16caf07`; this assessment reused them because the records-cache fingerprint matched (`records_cache_reused: true`; record-building code, screen and transcript caches unchanged) |
| classes | clean_reconstitution, tandem_duplication_or_insertion, gap_or_deletion, ambiguous_junction_register |
| screen rows in those classes | 50,226 (all five cohorts) |
| rows used | 50,113. The other 113 have a parent CDS with an internal UGA (selenoproteins such as GPX1/2/3/4/6), so the translated fusion stops early and no full domain exists; they are excluded and counted in the manifest |
| records | **3,034** chimeric domains (sequence + origin) over **2,948** distinct sequences |
| Pfam | 39.0 (2026-10-01), 408 per-family HMMs fetched from the InterPro API. PF12265 returns HTTP 410 (retired), so 1 record has HMM = not_assessed |
| UniProtKB | 2026_04, 1,461 genes, 41 without a reviewed human entry |
| AlphaFold DB | model v6, 1,420 accessions, 1,412 with a usable model, 0 fetch failures. 8 very large proteins have no single model (FBN1, HSPG2, COL7A1, PLEC, AHNAK, MUC5B, TTN, …) |
| tools | pyhmmer 0.12.3, gemmi 0.7.5, biopython 1.88, freesasa 2.2.1 |

Files:
- `chimera_records.json.gz`: one record per chimeric domain and origin, schema `domain_function_seq/2`, with every screen row kept as an observation.
- `domain_function_metrics.tsv.gz`: every metric, tier, fired rule, missing input and metric disagreement, per record.
- `tier_counts.json`: tier counts by record, by sample and by row, per cohort and class.
- `shortlist_for_structure.tsv`: input to PR 2.
- `run_manifest.json`: releases, parameters, rules, required inputs, the records fingerprint and the read-support rule.
- `report_tables.md`: the generated tables, repeated below.

## Summary

| tier (HEURISTIC) | records | share |
|---|---|---|
| likely_intact | 2,178 | 72 % |
| uncertain | 584 | 19 % |
| not_assessed | 125 | 4 % |
| likely_broken | 147 | 5 % |

Records per cohort (a record seen in several cohorts counts in each):

| cohort | likely_intact | uncertain | not_assessed | likely_broken | total |
|---|---|---|---|---|---|
| depmap_22q2 | 963 | 224 | 52 | 53 | 1,292 |
| msk_impact_50k_2026 (GRCh37) | 0 | 2 | 0 | 0 | 2 |
| target_consolidated | 168 | 108 | 31 | 29 | 336 |
| target_starfusion | 89 | 34 | 12 | 17 | 152 |
| tcga_gtex | 1,614 | 347 | 78 | 90 | 2,129 |

* **The junction class does not decide the tier.** Of 1,679 clean
  reconstitutions, 1,492 are likely_intact. So are 589 of 718
  ambiguous-register, 56 of 303 gap/deletion and 41 of 334 tandem/insertion
  records. Conversely, 19 clean reconstitutions are likely_broken: the chimera
  falls below GA, or an in-domain point site is lost. Most tandem and gap
  records are `uncertain` (270 of 334 and 165 of 303). The usual reasons are an
  indel residue in a parent helix or strand (U5, 199 records) and HMMER
  splitting the domain around the indel (U8, 244).
* **A clean junction can still lose a site.** For example, S100A10::S100A11
  takes S100A10's N-terminus, which lacks the loop residues that bind Ca²⁺ in
  S100A11 (the known reason S100A10 is calcium-insensitive). So S100A11's Ca²⁺
  site 33 is truly absent. AFM::ALB and AFP::ALB behave the same way. All 81 B3
  records rest on sites inside the Pfam domain.
* **Metrics often disagree**, and the disagreements are listed per record
  rather than resolved:
  - 864 likely_intact records break buried-core contacts in the SCHEMA-like count;
  - 231 pass GA with an indel residue in a parent helix or strand;
  - 222 pass GA with junction residues placed in HMM insert states;
  - 71 pass GA but lose a point site;
  - 95 fall below GA while every point key site is retained.
* **What likely_intact means here:**
  - the domain still scores as its Pfam family (≥ GA, ≥ 80 % match-state
    coverage, within 10 bits of the weaker parent);
  - no point key site and no motif-scale feature (≤ 10 residues) *inside the
    Pfam domain* is lost, altered or unresolved;
  - no indel touches parent helix/strand or a key site;
  - **both parent AlphaFold DB structures are present and validated**, and both
    parents' UniProt annotations are assessable.

  It does *not* mean:
  - longer annotated regions are whole: 267 likely_intact records list a
    longer region as lost or altered, reported in `key_span_lost_or_altered`;
  - flank sites survive: 182 records list a flank site as lost or altered, in
    `key_flank_features`;
  - there was anything to lose: in 417 records neither parent has an annotated
    feature inside the domain (`no_annotations_in_domain` on both sides). They
    were assessed, but "nothing lost" is vacuous there.

## Read support (corrected)

Each biological sample counts once. The two TARGET tables are two call sets
over the same TARGET RNA-seq with identical barcodes, so they share one
`TARGET:` sample namespace. A sample's read support is the **largest**
split-read count among its rows (tables, callers or event rows), never their
sum. `max_split_reads_one_sample` is the best single sample, and
`total_split_reads_over_samples` sums the per-sample maxima. The first report
treated each table's row as its own sample and summed reads over rows. Its
phrase "48 split reads across 19 samples" for GPC6::GPC5 was wrong: 48 is one
DepMap cell line (ACH-001552), and the record has 17 distinct samples.

## Named cases

Read support is given as best sample / sum over samples (samples, rows).

| fusion | records | outcome |
|---|---|---|
| **FLI1::ERG** | 1 (MSK, 1 sample, 0 split reads; read_evidence flagged) | **uncertain** (U2 bit drop −14.6/−16.4; U5 indel in helix). Tandem +7 (`IYYYDKN`). Ets domain 106.0 bits vs GA 26.9 (parents 120.6 / 122.4), coverage 0.99. 7 of the 9 junction residues sit in HMM insert states. The duplicated residues come from a 3' helix-adjacent segment (3' DSSP `HHHTTTT`). FLI1's ETS DNA-binding region is conserved (residues past the crossover are replaced by identical ERG residues); ERG's is swapped. SCHEMA E 0/0. Disagreement: it passes GA, but the insertion falls in helix and insert states. |
| **STAT5B::STAT5A** | 1 (MSK, 1 sample, 0 split reads; read_evidence and breakpoint_repeat flagged) | **uncertain** (U5, U7, U8). Tandem +16. STAT_bind domain 174.6 bits vs GA 27.0 (parents 175.9 / 174.8), coverage 1.0. HMMER splits the domain into more hits than either parent, and the best hit does not span the junction. Indel residues are loop/strand. No UniProt feature inside the domain; STAT5B's NMI-interaction region in the flank is retained. E 0/0. |
| **GPC4::GPC3** | 5 | **likely_broken** (gap; B1, B2, B3): 8.1 bits vs GA 26.4, 9 % coverage, and GPC3's annotated cleavage site 358–359 is lost (newly detected; 334c050 missed it). This is the best-supported record (25 / 57, 3 samples), shortlist #17. **likely_broken** (gap, 1 sample; B3): 345.4 bits, coverage 0.41, and loses the same cleavage site. **uncertain** (gap): 490.7 bits, split domain; 4 / 11, 6 samples, 8 rows (two TARGET samples appear in both tables). **likely_intact** ×2 (ambiguous register, 1 sample each): 614–618 bits. |
| **GPC6::GPC5** | 4 | **likely_broken** (gap): 10.2 bits vs GA 26.4, 8 % HMM coverage, i.e. only a glypican fragment. This is the best-supported record (48 / 84, 17 samples, 19 rows; the 48 is DepMap ACH-001552), shortlist #16. **likely_broken** (gap, 2 / 4, 3 samples). **uncertain** (gap). **likely_intact** (ambiguous register; DepMap, 1 sample). |
| **GPC3::GPC4** | 3 | **likely_intact** (ambiguous register): 781.9 bits; 1 / 1, **1 sample** in both TARGET tables (the first report said 2 samples). **uncertain** (tandem +458): 5 / 7, 2 samples; shortlist #8. **uncertain** (gap): 302.4 bits, coverage 0.44. |
| **LPO::MPO** | 2 (both clean_reconstitution) | **likely_intact** ×2. An_peroxidase 608.0 / 609.6 bits vs GA 25.0 (parents 602.8 / 627.5), coverage 0.998. All 22 in-domain heme, Ca²⁺ and active-site features are retained or conserved in the TCGA/GTEx record (3 / 12, 8 samples), shortlist #7. In the DepMap record (2 / 2, 1 sample), MPO's covalent-heme residue 409 and its autolysis site 409–410 carry LPO's aligned residue (`swapped`, allowed), so MPO's third heme link may not form. Disagreement: E 40–53 with E_core 23–38. |
| **PIM3::ULK1** | 1 (TARGET; no artefact flag fired) | **uncertain** (U2, U3, U5, U7, U8). Read support 2 / 2: **1 sample** seen in both TARGET tables (the first report said 2 samples). Tandem +77. Pkinase 147.1 vs GA 31.7 (parents 209.9 / 240.1), coverage 0.65, split domain. The indel spans helix, strand and loop. E 45/34, E_core 31/28. ULK1's ATP-binding site 22–30 is swapped for PIM3's. Shortlist #12. |

## Shortlist for PR 2

27 rows in `shortlist_for_structure.tsv`, one record per fusion and family:
- **18 `ranked`:** at least one supporting row passes all six pilot artefact
  flags.
- **7 `fill`:** a row passes every flag except `genomic_similarity` and
  `sample_burden`. Those two flags fire on 3,015 and 3,008 of the 3,034
  records, because shared-domain partners are paralogues, usually in segmental
  duplications, so they cannot separate candidates here. The fired flags are
  named on each row. All 7 are keratin or myosin pairs; treat them as
  artefact-suspect.
- **2 `named_case`:** FLI1::ERG and STAT5B::STAT5A.

Order: tier, then fewest broken contacts (max E over the two parent
structures), then GA margin, then best single-sample split reads and number of
samples. Likely_broken records with strong read support (GPC6::GPC5,
GPC4::GPC3) stay in the list on purpose: a structure prediction of a broken
call is as informative as one of an intact call. PR 2 can group shortlisted
records by `sequence_id` to predict each sequence once.

## Changes after review

Six blocking issues from the cross-vendor review were fixed, each with a
regression test that fails on 334c050:
1. Missing or rejected AlphaFold DB models gave likely_intact.
2. Absent key residues were reported as retained.
3. UniProt features with unusable coordinates were silently dropped.
4. Samples and reads were double-counted.
5. The record key ignored the junction's origin.
6. A stale records cache could be reused.

Two refinements followed from inspecting the new calls:
- Lost or altered motif-scale span features now block likely_intact (rule U9;
  the non-blocking review item).
- Key residues are tiered only inside the Pfam domain, not the ±30 flanks.

| | before (334c050) | after (16caf07) | change |
|---|---|---|---|
| records | 2,948 | 3,034 | +86 (same sequence from different origins, now separate) |
| distinct sequences | 2,948 | 2,948 | 0 |
| likely_intact | 2,162 | 2,178 | +16 |
| uncertain | 577 | 584 | +7 |
| not_assessed | 103 | 125 | +22 |
| likely_broken | 106 | 147 | +41 |

| class | intact | uncertain | not_assessed | broken | records |
|---|---|---|---|---|---|
| clean_reconstitution | 1449 → 1492 | 63 → 69 | 81 → 99 | 10 → 19 | 1603 → 1679 |
| tandem_duplication_or_insertion | 44 → 41 | 270 → 270 | 2 → 4 | 16 → 19 | 332 → 334 |
| gap_or_deletion | 61 → 56 | 168 → 165 | 2 → 2 | 70 → 80 | 301 → 303 |
| ambiguous_junction_register | 608 → 589 | 76 → 80 | 18 → 20 | 10 → 29 | 712 → 718 |

Where records moved (matched on fusion, family, sequence and junction
residues):

| before → after | records | cause |
|---|---|---|
| likely_intact → likely_broken | 25 | B3: an in-domain point site absent from the chimera with nothing at its aligned column (formerly reported as retained) |
| uncertain → likely_broken | 14 | B3, as above |
| not_assessed → likely_broken | 2 | B3, as above |
| likely_intact → uncertain | 15 | U9: a motif-scale feature lost or altered |
| not_assessed → uncertain | 2 | U9 |
| likely_intact → not_assessed | 10 | a missing or rejected parent structure. 10 of the 11 records the round-1 review cited go here (e.g. FBN1::FBN3: no model for either parent; PLEC::DSP: no model / sequence mismatch). The 11th, COL12A1::COL1A1, is now uncertain via U9, and its rejected structure is still listed in its missing inputs |
| uncertain → not_assessed | 1 | missing structure and indel secondary structure |
| (new origins) | 86 | 66 likely_intact, 15 not_assessed, 5 uncertain |

**Round 2.** One more blocking issue: a feature that was both lost and unresolved
counted only as lost, and for long regions loss is reported only, so the
unresolved residue could not block likely_intact. Unresolved in-domain positions
are now counted per record (`key_unresolved_positions`), independently of each
feature's worst status. Measured on this run, no record has an unresolved
in-domain position, so no tier, rule or shortlist row changed.

No likely_intact record now has a parent structure status other than `ok`,
or an unassessed UniProt annotation.

**Shortlist: unchanged.** The same 27 fusion+family groups appear in the same
order, with the same records; only their read-support columns changed (e.g.
PIM3::ULK1 1 sample, not 2). Every record that changed tier was outside the
shortlist or ineligible for it.

## Limits

* The tiers are heuristic and uncalibrated. The 50 % / 80 % coverage, 10-bit,
  5-residue and 10-residue motif thresholds are documented choices, not fitted
  values.
* Passing GA supports membership of the Pfam family. It does not show a correct
  fold or activity. HMMER absorbs inserted residues in insert states, so the
  insert count and the indel's secondary structure are reported alongside.
* SCHEMA E uses each parent's AlphaFold monomer model, with contacts restricted
  to pLDDT ≥ 70. E is not normalised and has no universal threshold, so it is
  not used for tiering. For divergent paralogues (LPO/MPO, glypicans) even a
  clean junction gives a large E.
* Key residues come only from UniProtKB/Swiss-Prot. They are tiered inside the
  Pfam domain span, so a site just outside an approximate Pfam boundary is
  only reported, in `key_flank_features`. Longer span features are reported,
  not tiered.
* The parent alignment that places absent residues (BLOSUM62, global, free end
  gaps) can misplace columns in low-complexity repeats such as collagen
  Gly-X-Y; such "lost" or "swapped" calls inherit that uncertainty. Survival of
  a parent's own residues does not depend on the alignment.
* No structure of the chimera itself is modelled here; that is PR 2. Stability
  and language-model scores (PROVEAN, ESM, INDELi, FoldX) are PR 3.
* HMMvar here is the difference in HMMER full-sequence bits between equal-flank
  windows. Windows of an indel chimera differ in length from the parent's.
* Only 2 MSK records are translated shared-domain chimeras (FLI1::ERG,
  STAT5B::STAT5A), and both have 0 split reads in the source table.
* Artefact flags come unchanged from the pilot. Flag counts are per row
  (observation), not per sample.

---

# Generated tables

## Tier counts (HEURISTIC)

Distinct chimeric-domain records (sequence + origin), all cohorts; 2948 distinct sequences:

| class | likely_intact | uncertain | not_assessed | likely_broken | total |
|---|---|---|---|---|---|
| clean_reconstitution | 1492 | 69 | 99 | 19 | 1679 |
| tandem_duplication_or_insertion | 41 | 270 | 4 | 19 | 334 |
| gap_or_deletion | 56 | 165 | 2 | 80 | 303 |
| ambiguous_junction_register | 589 | 80 | 20 | 29 | 718 |
| **all** | 2178 | 584 | 125 | 147 | 3034 |

Records per cohort and class (a record seen in several cohorts counts in each):

| cohort | class | likely_intact | uncertain | not_assessed | likely_broken | total |
|---|---|---|---|---|---|---|
| depmap_22q2 | ambiguous_junction_register | 197 | 38 | 5 | 15 | 255 |
| depmap_22q2 | clean_reconstitution | 732 | 36 | 44 | 7 | 819 |
| depmap_22q2 | gap_or_deletion | 22 | 62 | 0 | 24 | 108 |
| depmap_22q2 | tandem_duplication_or_insertion | 12 | 88 | 3 | 7 | 110 |
| msk_impact_50k_2026 | tandem_duplication_or_insertion | 0 | 2 | 0 | 0 | 2 |
| target_consolidated | ambiguous_junction_register | 18 | 15 | 2 | 2 | 37 |
| target_consolidated | clean_reconstitution | 123 | 11 | 27 | 6 | 167 |
| target_consolidated | gap_or_deletion | 11 | 36 | 2 | 17 | 66 |
| target_consolidated | tandem_duplication_or_insertion | 16 | 46 | 0 | 4 | 66 |
| target_starfusion | ambiguous_junction_register | 14 | 9 | 1 | 2 | 26 |
| target_starfusion | clean_reconstitution | 71 | 7 | 11 | 5 | 94 |
| target_starfusion | gap_or_deletion | 3 | 5 | 0 | 9 | 17 |
| target_starfusion | tandem_duplication_or_insertion | 1 | 13 | 0 | 1 | 15 |
| tcga_gtex | ambiguous_junction_register | 469 | 46 | 13 | 18 | 546 |
| tcga_gtex | clean_reconstitution | 1089 | 53 | 64 | 12 | 1218 |
| tcga_gtex | gap_or_deletion | 34 | 88 | 0 | 49 | 171 |
| tcga_gtex | tandem_duplication_or_insertion | 22 | 160 | 1 | 11 | 194 |

Distinct biological samples per cohort and class, by their record's tier (a TARGET sample in both TARGET tables appears in both rows here but once in record totals):

| cohort | class | likely_intact | uncertain | not_assessed | likely_broken | total |
|---|---|---|---|---|---|---|
| depmap_22q2 | ambiguous_junction_register | 1092 | 682 | 5 | 22 | 1801 |
| depmap_22q2 | clean_reconstitution | 8553 | 728 | 250 | 10 | 9541 |
| depmap_22q2 | gap_or_deletion | 76 | 124 | 0 | 42 | 242 |
| depmap_22q2 | tandem_duplication_or_insertion | 17 | 261 | 26 | 13 | 317 |
| msk_impact_50k_2026 | tandem_duplication_or_insertion | 0 | 2 | 0 | 0 | 2 |
| target_consolidated | ambiguous_junction_register | 62 | 27 | 2 | 4 | 95 |
| target_consolidated | clean_reconstitution | 741 | 16 | 52 | 9 | 818 |
| target_consolidated | gap_or_deletion | 17 | 69 | 4 | 26 | 116 |
| target_consolidated | tandem_duplication_or_insertion | 29 | 176 | 0 | 5 | 210 |
| target_starfusion | ambiguous_junction_register | 46 | 14 | 1 | 4 | 65 |
| target_starfusion | clean_reconstitution | 559 | 11 | 23 | 8 | 601 |
| target_starfusion | gap_or_deletion | 3 | 6 | 0 | 12 | 21 |
| target_starfusion | tandem_duplication_or_insertion | 1 | 92 | 0 | 1 | 94 |
| tcga_gtex | ambiguous_junction_register | 7540 | 567 | 240 | 391 | 8738 |
| tcga_gtex | clean_reconstitution | 22543 | 967 | 902 | 139 | 24551 |
| tcga_gtex | gap_or_deletion | 347 | 169 | 0 | 293 | 809 |
| tcga_gtex | tandem_duplication_or_insertion | 28 | 1243 | 1 | 13 | 1285 |

Screen rows (observations) per cohort and class, by their record's tier:

| cohort | class | likely_intact | uncertain | not_assessed | likely_broken | total |
|---|---|---|---|---|---|---|
| depmap_22q2 | ambiguous_junction_register | 1092 | 682 | 5 | 23 | 1802 |
| depmap_22q2 | clean_reconstitution | 8573 | 728 | 250 | 10 | 9561 |
| depmap_22q2 | gap_or_deletion | 76 | 125 | 0 | 42 | 243 |
| depmap_22q2 | tandem_duplication_or_insertion | 17 | 261 | 26 | 13 | 317 |
| msk_impact_50k_2026 | tandem_duplication_or_insertion | 0 | 2 | 0 | 0 | 2 |
| target_consolidated | ambiguous_junction_register | 62 | 27 | 2 | 4 | 95 |
| target_consolidated | clean_reconstitution | 741 | 16 | 52 | 9 | 818 |
| target_consolidated | gap_or_deletion | 17 | 69 | 4 | 26 | 116 |
| target_consolidated | tandem_duplication_or_insertion | 29 | 177 | 0 | 5 | 211 |
| target_starfusion | ambiguous_junction_register | 46 | 14 | 1 | 4 | 65 |
| target_starfusion | clean_reconstitution | 559 | 11 | 23 | 8 | 601 |
| target_starfusion | gap_or_deletion | 3 | 6 | 0 | 12 | 21 |
| target_starfusion | tandem_duplication_or_insertion | 1 | 92 | 0 | 1 | 94 |
| tcga_gtex | ambiguous_junction_register | 8031 | 785 | 240 | 391 | 9447 |
| tcga_gtex | clean_reconstitution | 22617 | 967 | 902 | 139 | 24625 |
| tcga_gtex | gap_or_deletion | 347 | 169 | 0 | 293 | 809 |
| tcga_gtex | tandem_duplication_or_insertion | 28 | 1244 | 1 | 13 | 1286 |

## Rules fired, missing inputs, metric disagreements (records)

| rule | records |
|---|---|
| B1_below_GA | 57 |
| B2_low_coverage | 52 |
| B3_point_feature_lost | 81 |
| I1_all_pass | 2178 |
| N1_missing | 125 |
| U1_parent_below_GA | 225 |
| U2_bit_drop | 147 |
| U3_partial_coverage | 163 |
| U5_indel_in_ss | 199 |
| U6_indel_near_key | 67 |
| U7_junction_outside_hit | 140 |
| U8_domain_split | 244 |
| U9_motif_lost_or_altered | 28 |

| missing input | records |
|---|---|
| hmm | 1 |
| indel_secondary_structure | 26 |
| key_annotation_3p:low_uniprot_identity | 13 |
| key_annotation_3p:no_uniprot_entry | 79 |
| key_annotation_5p:low_uniprot_identity | 7 |
| key_annotation_5p:no_uniprot_entry | 66 |
| structure_3p:afdb_sequence_mismatch | 8 |
| structure_3p:low_uniprot_identity | 13 |
| structure_3p:no_afdb_model | 16 |
| structure_3p:no_uniprot_entry | 79 |
| structure_5p:afdb_sequence_mismatch | 9 |
| structure_5p:low_uniprot_identity | 7 |
| structure_5p:no_afdb_model | 14 |
| structure_5p:no_uniprot_entry | 66 |

| disagreement | records |
|---|---|
| below_GA_but_key_points_retained | 95 |
| below_GA_but_no_broken_contacts | 82 |
| likely_intact_but_broken_core_contacts | 864 |
| passes_GA_but_indel_in_helix_or_strand | 231 |
| passes_GA_but_key_point_lost | 71 |
| passes_GA_with_junction_residues_in_insert_states | 222 |

| parent | structure status (AlphaFold DB) | records |
|---|---|---|
| 3p | afdb_sequence_mismatch | 13 |
| 3p | low_uniprot_identity | 14 |
| 3p | no_afdb_model | 19 |
| 3p | no_uniprot_entry | 81 |
| 3p | ok | 2907 |
| 5p | afdb_sequence_mismatch | 9 |
| 5p | low_uniprot_identity | 7 |
| 5p | no_afdb_model | 16 |
| 5p | no_uniprot_entry | 69 |
| 5p | ok | 2933 |

| parent | key-annotation status (UniProt) | records |
|---|---|---|
| 3p | annotated | 2265 |
| 3p | low_uniprot_identity | 14 |
| 3p | no_annotations_in_domain | 674 |
| 3p | no_uniprot_entry | 81 |
| 5p | annotated | 2237 |
| 5p | low_uniprot_identity | 7 |
| 5p | no_annotations_in_domain | 721 |
| 5p | no_uniprot_entry | 69 |

## Shortlist for PR 2 structure prediction

One record per fusion+family, ordered by tier, then lowest SCHEMA E (max of the two parent structures), then GA margin, then read support. `ranked`: a row passes all six artefact flags; `fill`: a row passes all but the near-universal genomic_similarity and sample_burden flags (fired flags named); `named_case`: appended named fusions.

| # | fusion | family | class | tier | bits vs GA (parents) | cov | E 5'/3' | indel SS | key residues | split reads: best sample / sum over samples (samples) | cohorts | artefact flags (rows flagged) | shortlist |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | KRT82::KRT7 | PF00038 Filament | clean reconstitution | likely_intact | 388.31 / GA 40.0 (5' 343.76, 3' 403.96) | 1.0 | 0/0 | no_indel | lost 1; swapped 8; retained 7; span lost/altered KRT7:Region:127-144(Linker 1) | 1 / 2 (2) | tcga_gtex | genomic_similarity=1, sample_burden=1 | ranked;flagged:genomic_similarity,sample_burden |
| 2 | YWHAB::YWHAZ | PF00244 14-3-3 | clean reconstitution | likely_intact | 365.32 / GA 33.20 (5' 366.91, 3' 365.88) | 1.0 | 0/0 | no_indel | conserved 2; retained 2 | 1 / 6 (6) | depmap_22q2, tcga_gtex | genomic_similarity=5, sample_burden=5 | ranked;flagged:genomic_similarity,sample_burden |
| 3 | PTN::MDK | PF01091 PTN_MK_C | clean reconstitution | likely_intact | 96.97 / GA 25.29 (5' 98.99, 3' 96.97) | 1.0 | 0/0 | no_indel | swapped 3; retained 1 | 3 / 22 (12) | depmap_22q2, tcga_gtex | genomic_similarity=11, sample_burden=11 | ranked;flagged:genomic_similarity,sample_burden |
| 4 | CRIP1::CRIP2 | PF00412 LIM | clean reconstitution | likely_intact | 54.54 / GA 22.89 (5' 57.47, 3' 50.12) | 0.983 | 2/3 | no_indel | no annotated feature in window | 5 / 10 (5) | depmap_22q2, tcga_gtex | genomic_similarity=5, sample_burden=5 | ranked;flagged:genomic_similarity,sample_burden |
| 5 | AFM::ALB | PF00273 Serum_albumin | clean reconstitution | likely_intact | 181.12 / GA 34.40 (5' 177.4, 3' 195.4) | 1.0 | 25/26 | no_indel | swapped 9; conserved 3; retained 5 | 2 / 3 (2) | tcga_gtex | genomic_similarity=1, sample_burden=1 | ranked;flagged:genomic_similarity,sample_burden |
| 6 | KLK3::PRSS8 | PF00089 Trypsin | ambiguous junction register | likely_intact | 221.94 / GA 20.60 (5' 227.5, 3' 247.06) | 1.0 | 34/42 | no_indel | conserved 3; retained 3 | 1 / 10 (10) | tcga_gtex | genomic_similarity=9, sample_burden=9 | ranked;flagged:genomic_similarity,sample_burden |
| 7 | LPO::MPO | PF03098 An_peroxidase | clean reconstitution | likely_intact | 609.61 / GA 25.0 (5' 602.76, 3' 627.48) | 0.998 | 40/43 | no_indel | conserved 10; retained 12 | 3 / 12 (8) | tcga_gtex | genomic_similarity=5, sample_burden=5 | ranked;flagged:genomic_similarity,sample_burden |
| 8 | GPC3::GPC4 | PF01153 Glypican | tandem duplication or insertion | uncertain | 585.44 / GA 26.39 (5' 610.0, 3' 788.15) | 0.895 | 0/0 | helix|loop|strand | retained 1 | 5 / 7 (2) | depmap_22q2, tcga_gtex | genomic_similarity=1, sample_burden=1 | ranked;flagged:genomic_similarity,sample_burden |
| 9 | VIM::KRT18 | PF00038 Filament | clean reconstitution | uncertain | 354.58 / GA 40.0 (5' 382.38, 3' 370.36) | 0.994 | 0/0 | no_indel | lost 5; swapped 6; conserved 1; retained 9; span lost/altered KRT18:Region:70-373(Necessary for interaction with PNN);KRT18:Region:77-128(Interaction with TRADD);KRT18:Region:80-115(Coil 1A);KRT18:Region:116-132(Linker 1);KRT18:Region:133-224(Coil 1B) | 1 / 1 (1) | tcga_gtex | none flagged | ranked |
| 10 | COL3A1::COL1A1 | PF01391 Collagen | gap or deletion | uncertain | 47.53 / GA 27.20 (5' 35.91, 3' 45.04) | 0.983 | 0/0 | loop | lost 2; retained 1; span lost/altered COL3A1:Region:168-1196(Triple-helical region);COL1A1:Region:179-1192(Triple-helical region) | 1 / 1 (1) | tcga_gtex | none flagged | ranked |
| 11 | SPNS2::SPNS3 | PF07690 MFS_1 | tandem duplication or insertion | uncertain | 128.47 / GA 33.40 (5' 134.31, 3' 137.17) | 0.914 | 10/7 | helix|loop | no annotated feature in window | 1 / 1 (1) | tcga_gtex | none flagged | ranked |
| 12 | PIM3::ULK1 | PF00069 Pkinase | tandem duplication or insertion | uncertain | 147.05 / GA 31.70 (5' 209.86, 3' 240.06) | 0.653 | 45/34 | helix|loop|strand | swapped 1; conserved 1; retained 4 | 2 / 2 (1) | target_consolidated, target_starfusion | none flagged | ranked |
| 13 | GNAL::GNAS | PF00503 G-alpha | clean reconstitution | not_assessed | 398.71 / GA 48.79 (5' 398.62, 3' 393.33) | 0.997 | 0/not_assessed | no_indel | not assessed (3p: low_uniprot_identity) | 1 / 1 (1) | tcga_gtex | none flagged | ranked |
| 14 | COL1A1::COL18A1 | PF01391 Collagen | tandem duplication or insertion | likely_broken | 26.18 / GA 27.20 (5' 36.38, 3' 38.22) | 0.867 | 0/0 | loop | lost 2; span lost/altered COL1A1:Region:179-1192(Triple-helical region);COL18A1:Region:900-1021(Triple-helical region 3 (COL3)) | 2 / 4 (3) | target_consolidated, target_starfusion, tcga_gtex | genomic_similarity=2, sample_burden=2 | ranked;flagged:genomic_similarity,sample_burden |
| 15 | COL1A1::COL12A1 | PF01391 Collagen | gap or deletion | likely_broken | 25.61 / GA 27.20 (5' 39.54, 3' 35.57) | 0.967 | 0/not_assessed | not_assessed | lost 2; retained 1; span lost/altered COL1A1:Region:179-1192(Triple-helical region);COL12A1:Region:2747-2898(Triple-helical region (COL2) with 1 imperfection) | 3 / 63 (58) | target_consolidated, target_starfusion, tcga_gtex | genomic_similarity=56, sample_burden=56 | ranked;flagged:genomic_similarity,sample_burden |
| 16 | GPC6::GPC5 | PF01153 Glypican | gap or deletion | likely_broken | 10.16 / GA 26.39 (5' 840.44, 3' 618.29) | 0.083 | 0/0 | no_indel | no annotated feature in window | 48 / 84 (17) | depmap_22q2, target_consolidated, target_starfusion, tcga_gtex | genomic_similarity=18, sample_burden=18 | ranked;flagged:genomic_similarity,sample_burden |
| 17 | GPC4::GPC3 | PF01153 Glypican | gap or deletion | likely_broken | 8.07 / GA 26.39 (5' 788.15, 3' 610.0) | 0.086 | 0/0 | no_indel | lost 1; POINT LOST GPC3:Site:358-359(Cleavage) | 25 / 57 (3) | depmap_22q2, tcga_gtex | genomic_similarity=2, sample_burden=2 | ranked;flagged:genomic_similarity,sample_burden |
| 18 | CMTM7::CMTM8 | PF01284 MARVEL | gap or deletion | likely_broken | 29.89 / GA 32.79 (5' 75.67, 3' 81.35) | 0.444 | 8/7 | helix|loop | no annotated feature in window | 4 / 5 (2) | tcga_gtex | genomic_similarity=1, sample_burden=1 | ranked;flagged:genomic_similarity,sample_burden |
| 19 | MYH10::MYH14 | PF00063 Myosin_head | clean reconstitution | likely_intact | 1013.82 / GA 33.29 (5' 1013.61, 3' 945.32) | 1.0 | 0/0 | no_indel | swapped 1; conserved 1; retained 2 | 1 / 1 (1) | tcga_gtex | genomic_similarity=1, sample_burden=1 | fill;flagged:genomic_similarity,sample_burden |
| 20 | KRT13::KRT5 | PF00038 Filament | ambiguous junction register | likely_intact | 415.35 / GA 40.0 (5' 409.22, 3' 411.45) | 1.0 | 0/0 | no_indel | lost 3; swapped 6; retained 6; span lost/altered KRT5:Region:168-203(Coil 1A);KRT5:Region:204-222(Linker 1);KRT5:Region:316-338(Linker 12) | 7 / 7 (1) | depmap_22q2 | genomic_similarity=1, sample_burden=1 | fill;flagged:genomic_similarity,sample_burden |
| 21 | KRT13::KRT7 | PF00038 Filament | clean reconstitution | likely_intact | 413.58 / GA 40.0 (5' 409.22, 3' 403.96) | 1.0 | 0/0 | no_indel | lost 2; swapped 8; retained 6; span lost/altered KRT7:Region:127-144(Linker 1);KRT7:Region:237-260(Linker 12) | 4 / 4 (1) | depmap_22q2 | genomic_similarity=1, sample_burden=1 | fill;flagged:genomic_similarity,sample_burden |
| 22 | KRT8::KRT5 | PF00038 Filament | clean reconstitution | likely_intact | 412.91 / GA 40.0 (5' 408.94, 3' 411.45) | 1.0 | 0/0 | no_indel | lost 1; swapped 8; conserved 1; retained 7; span lost/altered KRT5:Region:204-222(Linker 1) | 8 / 91 (46) | depmap_22q2, tcga_gtex | genomic_similarity=46, sample_burden=46 | fill;flagged:genomic_similarity,sample_burden |
| 23 | KRT19::KRT7 | PF00038 Filament | clean reconstitution | likely_intact | 408.09 / GA 40.0 (5' 407.15, 3' 403.96) | 1.0 | 0/0 | no_indel | swapped 11; conserved 1; retained 7 | 6 / 13 (3) | depmap_22q2 | genomic_similarity=3, sample_burden=3 | fill;flagged:genomic_similarity,sample_burden |
| 24 | KRT5::KRT8 | PF00038 Filament | clean reconstitution | likely_intact | 407.07 / GA 40.0 (5' 411.45, 3' 408.94) | 1.0 | 0/0 | no_indel | swapped 8; conserved 1; retained 8 | 5 / 8 (2) | depmap_22q2 | genomic_similarity=2, sample_burden=2 | fill;flagged:genomic_similarity,sample_burden |
| 25 | KRT5::KRT13 | PF00038 Filament | ambiguous junction register | likely_intact | 405.32 / GA 40.0 (5' 411.45, 3' 409.22) | 0.997 | 0/0 | no_indel | lost 1; swapped 7; retained 7; span lost/altered KRT13:Region:104-139(Coil 1A) | 3 / 3 (1) | depmap_22q2 | genomic_similarity=1, sample_burden=1 | fill;flagged:genomic_similarity,sample_burden |
| 26 | FLI1::ERG | PF00178 Ets | tandem duplication or insertion | uncertain | 106.01 / GA 26.89 (5' 120.58, 3' 122.37) | 0.988 | 0/0 | helix|loop | swapped 1; conserved 1 | 0 / 0 (1) | msk_impact_50k_2026 | genomic_similarity=1, sample_burden=1, read_evidence=1 | named_case;flagged:genomic_similarity,sample_burden,read_evidence |
| 27 | STAT5B::STAT5A | PF02864 STAT_bind | tandem duplication or insertion | uncertain | 174.64 / GA 27.0 (5' 175.88, 3' 174.81) | 1.0 | 0/0 | loop|strand | no annotated feature in window | 0 / 0 (1) | msk_impact_50k_2026 | genomic_similarity=1, breakpoint_repeat=1, read_evidence=1 | named_case;flagged:genomic_similarity,breakpoint_repeat,read_evidence |

## Named cases (every record)


**FLI1::ERG** (1 record(s))

| # | fusion | family | class | tier | bits vs GA (parents) | cov | E 5'/3' | indel SS | key residues | split reads: best sample / sum over samples (samples) | cohorts | artefact flags (rows flagged) | shortlist |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|  | FLI1::ERG | PF00178 Ets | tandem duplication or insertion | uncertain | 106.01 / GA 26.89 (5' 120.58, 3' 122.37) | 0.988 | 0/0 | helix|loop | swapped 1; conserved 1 | 0 / 0 (1) | msk_impact_50k_2026 | genomic_similarity=1, sample_burden=1, read_evidence=1 |  |

- `chim_e37f8e6307348429` tandem_duplication_or_insertion: rules `U2_bit_drop;U5_indel_in_ss`; missing `-`; disagreements `passes_GA_but_indel_in_helix_or_strand;passes_GA_with_junction_residues_in_insert_states`; aligned indel +7 (IYYYDKN)/-0 cols; junction residues in HMM inserts 7/9; delta bits 5'/3' -14.57/-16.36; HMMvar 5'/3' -7.57/-9.22; junction DSSP 5' `TTTTSEE` 3' `HHHTTTT`; E_core 5'/3' 0/0

**STAT5B::STAT5A** (1 record(s))

| # | fusion | family | class | tier | bits vs GA (parents) | cov | E 5'/3' | indel SS | key residues | split reads: best sample / sum over samples (samples) | cohorts | artefact flags (rows flagged) | shortlist |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|  | STAT5B::STAT5A | PF02864 STAT_bind | tandem duplication or insertion | uncertain | 174.64 / GA 27.0 (5' 175.88, 3' 174.81) | 1.0 | 0/0 | loop|strand | no annotated feature in window | 0 / 0 (1) | msk_impact_50k_2026 | genomic_similarity=1, breakpoint_repeat=1, read_evidence=1 |  |

- `chim_a44830a515dbd281` tandem_duplication_or_insertion: rules `U7_junction_outside_hit;U8_domain_split;U5_indel_in_ss`; missing `-`; disagreements `passes_GA_but_indel_in_helix_or_strand`; aligned indel +16 (TKQPPQVLKTQTKFAA)/-0 cols; junction residues in HMM inserts 0/18; delta bits 5'/3' -1.23/-0.16; HMMvar 5'/3' 9.17/10.24; junction DSSP 5' `EEEEEEE` 3' `EEEEPPP`; E_core 5'/3' 0/0

**GPC4::GPC3** (5 record(s))

| # | fusion | family | class | tier | bits vs GA (parents) | cov | E 5'/3' | indel SS | key residues | split reads: best sample / sum over samples (samples) | cohorts | artefact flags (rows flagged) | shortlist |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|  | GPC4::GPC3 | PF01153 Glypican | gap or deletion | likely_broken | 345.41 / GA 26.39 (5' 788.15, 3' 610.0) | 0.407 | 8/0 | helix|loop|strand | lost 1; POINT LOST GPC3:Site:358-359(Cleavage) | 1 / 1 (1) | depmap_22q2 | genomic_similarity=1, sample_burden=1 |  |
|  | GPC4::GPC3 | PF01153 Glypican | ambiguous junction register | likely_intact | 613.73 / GA 26.39 (5' 788.15, 3' 610.0) | 0.993 | 12/9 | no_indel | retained 1 | 2 / 2 (1) | tcga_gtex | genomic_similarity=1, sample_burden=1, normal_recurrence=1 |  |
|  | GPC4::GPC3 | PF01153 Glypican | gap or deletion | likely_broken | 8.07 / GA 26.39 (5' 788.15, 3' 610.0) | 0.086 | 0/0 | no_indel | lost 1; POINT LOST GPC3:Site:358-359(Cleavage) | 25 / 57 (3) | depmap_22q2, tcga_gtex | genomic_similarity=2, sample_burden=2 |  |
|  | GPC4::GPC3 | PF01153 Glypican | ambiguous junction register | likely_intact | 618.41 / GA 26.39 (5' 788.15, 3' 610.0) | 0.993 | 48/54 | no_indel | retained 1 | 1 / 1 (1) | tcga_gtex | genomic_similarity=1, sample_burden=1 |  |
|  | GPC4::GPC3 | PF01153 Glypican | gap or deletion | uncertain | 490.71 / GA 26.39 (5' 788.15, 3' 610.0) | 0.831 | 3/4 | no_indel | retained 1 | 4 / 11 (6) | depmap_22q2, target_consolidated, target_starfusion, tcga_gtex | genomic_similarity=8, sample_burden=8 |  |

- `chim_426f4dd887cffb77` gap_or_deletion: rules `B3_point_feature_lost`; missing `-`; disagreements `passes_GA_but_key_point_lost;passes_GA_but_indel_in_helix_or_strand`; aligned indel +0 (-)/-242 cols; junction residues in HMM inserts 0/6; delta bits 5'/3' -442.74/-264.59; HMMvar 5'/3' -416.08/-237.97; junction DSSP 5' `HHTTPPP` 3' `HHHHHHH`; E_core 5'/3' 5/0

- `chim_4274ed455d6f850a` ambiguous_junction_register: rules `I1_all_pass`; missing `-`; disagreements `likely_intact_but_broken_core_contacts`; aligned indel +0 (-)/-0 cols; junction residues in HMM inserts 0/6; delta bits 5'/3' -174.43/3.73; HMMvar 5'/3' -174.42/3.69; junction DSSP 5' `E-EE-TT` 3' `E--TT-S`; E_core 5'/3' 0/3

- `chim_5cb7d12fa77dce03` gap_or_deletion: rules `B1_below_GA;B2_low_coverage;B3_point_feature_lost`; missing `-`; disagreements `below_GA_but_no_broken_contacts`; aligned indel +0 (-)/-0 cols; junction residues in HMM inserts 0/6; delta bits 5'/3' -780.08/-601.93; HMMvar 5'/3' -756.01/-577.9; junction DSSP 5' `E-EE-TT` 3' `TS-SSS-`; E_core 5'/3' 0/0

- `chim_de2f675a496bc89f` ambiguous_junction_register: rules `I1_all_pass`; missing `-`; disagreements `likely_intact_but_broken_core_contacts`; aligned indel +0 (-)/-0 cols; junction residues in HMM inserts 0/6; delta bits 5'/3' -169.74/8.41; HMMvar 5'/3' -169.73/8.38; junction DSSP 5' `HHHHHHH` 3' `HHHHHHH`; E_core 5'/3' 20/28

- `chim_f1b9207cfcabc3b5` gap_or_deletion: rules `U2_bit_drop;U8_domain_split`; missing `-`; disagreements `-`; aligned indel +0 (-)/-0 cols; junction residues in HMM inserts 0/6; delta bits 5'/3' -297.45/-119.29; HMMvar 5'/3' -278.95/-100.84; junction DSSP 5' `E-EE-TT` 3' `HHHHHHH`; E_core 5'/3' 0/2

**GPC6::GPC5** (4 record(s))

| # | fusion | family | class | tier | bits vs GA (parents) | cov | E 5'/3' | indel SS | key residues | split reads: best sample / sum over samples (samples) | cohorts | artefact flags (rows flagged) | shortlist |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|  | GPC6::GPC5 | PF01153 Glypican | ambiguous junction register | likely_intact | 630.94 / GA 26.39 (5' 840.44, 3' 618.29) | 0.975 | 17/10 | no_indel | no annotated feature in window | 1 / 1 (1) | depmap_22q2 | genomic_similarity=1, sample_burden=1 |  |
|  | GPC6::GPC5 | PF01153 Glypican | gap or deletion | likely_broken | 122.1 / GA 26.39 (5' 840.44, 3' 618.29) | 0.157 | 0/0 | no_indel | no annotated feature in window | 2 / 4 (3) | depmap_22q2, target_consolidated | genomic_similarity=3, sample_burden=3 |  |
|  | GPC6::GPC5 | PF01153 Glypican | gap or deletion | likely_broken | 10.16 / GA 26.39 (5' 840.44, 3' 618.29) | 0.083 | 0/0 | no_indel | no annotated feature in window | 48 / 84 (17) | depmap_22q2, target_consolidated, target_starfusion, tcga_gtex | genomic_similarity=18, sample_burden=18 |  |
|  | GPC6::GPC5 | PF01153 Glypican | gap or deletion | uncertain | 490.85 / GA 26.39 (5' 840.44, 3' 618.29) | 0.488 | 0/0 | no_indel | no annotated feature in window | 1 / 1 (1) | tcga_gtex | genomic_similarity=1, sample_burden=1 |  |

- `chim_36fcb16b2b0b1744` ambiguous_junction_register: rules `I1_all_pass`; missing `-`; disagreements `likely_intact_but_broken_core_contacts`; aligned indel +0 (-)/-0 cols; junction residues in HMM inserts 0/6; delta bits 5'/3' -209.5/12.65; HMMvar 5'/3' -209.52/12.68; junction DSSP 5' `E-EE-TT` 3' `E--TT-S`; E_core 5'/3' 2/3

- `chim_4ef2912a09cba631` gap_or_deletion: rules `B2_low_coverage`; missing `-`; disagreements `-`; aligned indel +0 (-)/-0 cols; junction residues in HMM inserts 0/6; delta bits 5'/3' -718.34/-496.2; HMMvar 5'/3' -709.58/-487.38; junction DSSP 5' `HHHHHHH` 3' `-S---S-`; E_core 5'/3' 0/0

- `chim_cad4f6a750d94835` gap_or_deletion: rules `B1_below_GA;B2_low_coverage`; missing `-`; disagreements `below_GA_but_key_points_retained;below_GA_but_no_broken_contacts`; aligned indel +0 (-)/-0 cols; junction residues in HMM inserts 0/6; delta bits 5'/3' -830.29/-608.14; HMMvar 5'/3' -801.6/-579.4; junction DSSP 5' `E-EE-TT` 3' `-S---S-`; E_core 5'/3' 0/0

- `chim_ddc27a1899220c5f` gap_or_deletion: rules `U2_bit_drop;U3_partial_coverage;U7_junction_outside_hit;U8_domain_split`; missing `-`; disagreements `-`; aligned indel +0 (-)/-0 cols; junction residues in HMM inserts 0/6; delta bits 5'/3' -349.59/-127.44; HMMvar 5'/3' -343.08/-120.88; junction DSSP 5' `HHHHHHH` 3' `-S---S-`; E_core 5'/3' 0/0

**GPC3::GPC4** (3 record(s))

| # | fusion | family | class | tier | bits vs GA (parents) | cov | E 5'/3' | indel SS | key residues | split reads: best sample / sum over samples (samples) | cohorts | artefact flags (rows flagged) | shortlist |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|  | GPC3::GPC4 | PF01153 Glypican | tandem duplication or insertion | uncertain | 585.44 / GA 26.39 (5' 610.0, 3' 788.15) | 0.895 | 0/0 | helix|loop|strand | retained 1 | 5 / 7 (2) | depmap_22q2, tcga_gtex | genomic_similarity=1, sample_burden=1 |  |
|  | GPC3::GPC4 | PF01153 Glypican | gap or deletion | uncertain | 302.39 / GA 26.39 (5' 610.0, 3' 788.15) | 0.44 | 15/14 | no_indel | conserved 1 | 1 / 1 (1) | depmap_22q2 | genomic_similarity=1, sample_burden=1 |  |
|  | GPC3::GPC4 | PF01153 Glypican | ambiguous junction register | likely_intact | 781.9 / GA 26.39 (5' 610.0, 3' 788.15) | 0.95 | 56/48 | no_indel | conserved 1 | 1 / 1 (1) | target_consolidated, target_starfusion | genomic_similarity=2, sample_burden=2 |  |

- `chim_3df47c1f84878b99` tandem_duplication_or_insertion: rules `U2_bit_drop;U7_junction_outside_hit;U8_domain_split;U5_indel_in_ss`; missing `-`; disagreements `passes_GA_but_indel_in_helix_or_strand`; aligned indel +458 (GDHLKICPQGSTCCSQEMEEKYSLQSKDDFKSVVSEQCNHLQAVFASRYKKFDEFFKELLENAEKSLNDMFVKTYGHLYMQNSELFKDLFVELKRYYVVGNVNLEEMLNDFWARLLERMFRLVNSQYHFTDEYLECVSKYTEQLKPFGDVPRKLKLQVTRAFVAARTFAQGLAVAGDVVSKVSVVNPTAQCTHALLKMIYCSHCRGLVTVKPCYNYCSNIMRGCLANQGDLDFEWNNFIDAMLMVAERLEGPFNIESVMDPIDVKISDAIMNMQDNSVQVSQKVFQGCGPPKPLPAGRISRSISESAFSARFRPHHPEERPTTAAGTSLDRLVTDVKEKLKQAKKFWSSLPSNVCNDERMAAGNGNEDDCWNGKGKSRYLFAVTGNGLANQGNNPEVQVDTSKPDILILRQIMALRVMTSKMKNAYNGNDVDFFDISDESSGEGSGSGCEYQQCPSEF)/-4 cols; junction residues in HMM inserts 0/460; delta bits 5'/3' -24.56/-202.72; HMMvar 5'/3' 717.19/539.08; junction DSSP 5' `-TTS-SS` 3' `EE-TT-S`; E_core 5'/3' 0/0

- `chim_5d0e757b69d5e4e6` gap_or_deletion: rules `U2_bit_drop;U3_partial_coverage;U8_domain_split`; missing `-`; disagreements `-`; aligned indel +0 (-)/-0 cols; junction residues in HMM inserts 0/6; delta bits 5'/3' -307.61/-485.76; HMMvar 5'/3' -207.84/-385.95; junction DSSP 5' `HHHHHHH` 3' `HHHHHHH`; E_core 5'/3' 7/7

- `chim_6014defa0203ac26` ambiguous_junction_register: rules `I1_all_pass`; missing `-`; disagreements `likely_intact_but_broken_core_contacts`; aligned indel +0 (-)/-0 cols; junction residues in HMM inserts 0/6; delta bits 5'/3' 171.9/-6.26; HMMvar 5'/3' 171.89/-6.22; junction DSSP 5' `HHHHHHH` 3' `HHHHHHH`; E_core 5'/3' 28/20

**LPO::MPO** (2 record(s))

| # | fusion | family | class | tier | bits vs GA (parents) | cov | E 5'/3' | indel SS | key residues | split reads: best sample / sum over samples (samples) | cohorts | artefact flags (rows flagged) | shortlist |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|  | LPO::MPO | PF03098 An_peroxidase | clean reconstitution | likely_intact | 609.61 / GA 25.0 (5' 602.76, 3' 627.48) | 0.998 | 40/43 | no_indel | conserved 10; retained 12 | 3 / 12 (8) | tcga_gtex | genomic_similarity=5, sample_burden=5 |  |
|  | LPO::MPO | PF03098 An_peroxidase | clean reconstitution | likely_intact | 607.98 / GA 25.0 (5' 602.76, 3' 627.48) | 0.998 | 53/53 | no_indel | swapped 2; conserved 10; retained 10 | 2 / 2 (1) | depmap_22q2 | paralog_sequence=1, genomic_similarity=1, sample_burden=1 |  |

- `chim_4277e22644ac73b3` clean_reconstitution: rules `I1_all_pass`; missing `-`; disagreements `likely_intact_but_broken_core_contacts`; aligned indel +0 (-)/-0 cols; junction residues in HMM inserts 0/6; delta bits 5'/3' 6.85/-17.87; HMMvar 5'/3' 6.86/-17.87; junction DSSP 5' `PEP-SST` 3' `P-SSGGG`; E_core 5'/3' 23/28

- `chim_d562ab789b6093f7` clean_reconstitution: rules `I1_all_pass`; missing `-`; disagreements `likely_intact_but_broken_core_contacts`; aligned indel +0 (-)/-0 cols; junction residues in HMM inserts 0/6; delta bits 5'/3' 5.22/-19.5; HMMvar 5'/3' 5.22/-19.51; junction DSSP 5' `HHHHHHI` 3' `HHHHIII`; E_core 5'/3' 36/38

**PIM3::ULK1** (1 record(s))

| # | fusion | family | class | tier | bits vs GA (parents) | cov | E 5'/3' | indel SS | key residues | split reads: best sample / sum over samples (samples) | cohorts | artefact flags (rows flagged) | shortlist |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|  | PIM3::ULK1 | PF00069 Pkinase | tandem duplication or insertion | uncertain | 147.05 / GA 31.70 (5' 209.86, 3' 240.06) | 0.653 | 45/34 | helix|loop|strand | swapped 1; conserved 1; retained 4 | 2 / 2 (1) | target_consolidated, target_starfusion | none flagged |  |

- `chim_97417cac3653327e` tandem_duplication_or_insertion: rules `U2_bit_drop;U3_partial_coverage;U7_junction_outside_hit;U8_domain_split;U5_indel_in_ss`; missing `-`; disagreements `passes_GA_but_indel_in_helix_or_strand;passes_GA_with_junction_residues_in_insert_states`; aligned indel +77 (AMRTLSEDTIRLFLQQIAGAMRLLHSKGIIHRDLKPQNILLSNPAGRRANPNSIRVKIADFGFARYLQSNMMAATLC)/-0 cols; junction residues in HMM inserts 9/79; delta bits 5'/3' -62.81/-93.01; HMMvar 5'/3' 58.91/28.86; junction DSSP 5' `----S-G` 3' `HHHHSS-`; E_core 5'/3' 31/28

## Data coverage

- Records: 3034 from 50226 translated-class screen rows; rebuild problems: {'parent CDS has an internal stop codon (selenocysteine UGA?)': 113}.
- Pfam 39.0 (2026-10-01T00:00:00Z), 408 family HMMs, unavailable: 1.
- UniProt 2026_04: 1461 genes, 41 without a reviewed entry.
- AlphaFold DB model versions ['6']: 1420 accessions, 8 without a usable model, 0 fetch/parse problems.
