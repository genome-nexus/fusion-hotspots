# Structure prediction of the shortlisted chimeric domains (2026-11-13)

PR 2 of the domain-function stack. Input: the 27-record structure shortlist
and records of `runs/domain_function_seq_2026-11-11/` (schema
`domain_function_seq/2`; sha256 of both inputs in `run_manifest.json`).
Method, definitions and thresholds: `docs/domain_structure.md`.

**Read this first.** Most of these chimeras join near-identical paralogues.
An AlphaFold-family model fed essentially the parents' MSA will usually
predict a confident, parent-like fold whatever the junction does. A
`fold_supported` flag is therefore **weak evidence** that the chimera folds;
it is a screen-out flag, and only `fold_doubtful` carries real information.
Nothing here is a calibrated predictor of folding or activity, and PR 1's
tiers are unchanged (`pr1_tier` is copied; it matches PR 1's
`domain_function_metrics.tsv.gz` for all 27 records).

## Run

* ColabFold 1.6.3 (alphafold-colabfold 2.3.20, JAX 0.10.2, CPU) on an Apple
  M3, 16 GB; AlphaFold2 weights `alphafold_params_2021-07-14`
  (`params_model_1_ptm`, sha256 in the manifest). Install and weights live in
  `~/.cache/cfh_colabfold` (outside the repository, ~2.4 GB).
* Settings, identical for chimeras and parents (`settings_id`
  `set_f6b2876d27c8`, which includes the MSA source): `alphafold2_ptm`, **1 model** (model_1_ptm), **1
  recycle**, **templates off**, **seed 0**, 1 seed, MSA depth `--max-msa
  256:512`, no Amber relax. MSAs: `mmseqs2_uniref_env` on the public server
  `https://api.colabfold.com`, 8 sequences per ticket with 20 s pauses, each
  a3m cached by sequence and never re-queried.
* Jobs: 71 distinct sequences (27 chimeric windows, 44 parent windows; 22,812
  residues). Each chimera window is its Pfam domain ±30 residues, exactly the
  sequence-step window; each parent window is the parent Pfam domain ±30.
* Wall time: MSAs ~53 min (9 tickets); predictions **11.8 h** in total (23 s
  for 113 residues to 82 min for the 1,027-residue GPC3::GPC4 window); the
  `predict` invocation took 12.7 h end to end. `analyze` takes ~30 s.
* Outcome: **71/71 predictions validated, 0 failed**, 0 MSA errors. The
  only line of `predict.log` matching `error|fail` is the closing summary
  `predict finished: 71 ok, 0 failed`.
* Checkpoint migration (no re-prediction): the predictions were made under
  key `set_e20d3cb4e3e2`, which did not include the MSA server. After the
  #160 review the key includes it. The 71 checkpoints were moved to
  `set_f6b2876d27c8` only after the provenance check passed. The work
  directory holds exactly one MSA cache, `msa_c3444ece` = sha256
  (`mmseqs2_uniref_env|https://api.colabfold.com`)[:8], and the predictor
  reads MSAs only from its own config's cache. Every job re-validated and has
  its a3m there. Each checkpoint records `migrated_from`, and the manifest
  counts 71 migrated.
* Committed: the metrics TSV, flag counts, manifest, this report, and the
  rank-1 model plus pLDDT/PAE (PAE rounded to whole Å; metrics use full
  precision) for the 19 distinct sequences of the named cases below (~2.5
  MB). The other 52 models stay in the work cache; MSAs are not committed.

## Flags

| flag | records |
|---|---|
| fold_supported | 14 |
| fold_doubtful | 9 |
| not_assessed | 4 |

Every `not_assessed` has a real missing input:

* COL3A1::COL1A1, COL1A1::COL18A1 and CMTM7::CMTM8 fire
  `N1_parents_not_confident`. Neither parent's domain pLDDT reaches 70 under
  these settings (53/60, 47/47, 58/64), so there is no confident reference
  fold to compare against.
* GNAL::GNAS fires `N3_afdb_comparison_missing`. The 3' AFDB comparison was
  refused because the UniProt entry PR 1 mapped for GNAS (P84996) matches the
  parent window at 0.3 % identity; see below.

No record was `not_assessed` for a failed prediction, MSA or mapping.

Rules fired (a record can fire several): D3 junction very low 8, D4 junction
pLDDT drop 7, D5 halves not packed 5, D1 domain pLDDT low 3, D2 domain pLDDT
drop 3, D6 fold differs 3, N1 3, N3 1, S1 all pass 14.

By PR 1 tier (unchanged): likely_intact 13 supported / 1 doubtful
(KRT5::KRT13); uncertain 1 supported (VIM::KRT18) / 5 doubtful / 1 not
assessed; likely_broken 3 doubtful / 2 not assessed; PR 1 not_assessed 1
not assessed (GNAL::GNAS).

By cohort (a record counts in each cohort it was seen in): tcga_gtex 10 / 5 /
4 (supported / doubtful / not assessed), depmap_22q2 8 / 4 / 0,
target_consolidated 0 / 3 / 1, target_starfusion 0 / 3 / 1,
msk_impact_50k_2026 0 / 2 / 0.

## Per-record measures

C = chimera, 5'/3' = the parent predictions under the same settings. pLDDT
means over the Pfam domain ("dom") and the junction window ("jx", junction
±10 residues, widened past any inserted residues); Δ = chimera minus the
better parent. "half PAE" = mean PAE between the 5'- and 3'-derived halves
(parents: their own residues before/after their junction). TM = mapped
TM-score of the chimera domain against each parent prediction and against the
AFDB parent model, through PR 1's parent alignment, normalised by the parent
domain; a pair counts only when both residues are inside their own Pfam
domains (flank residues never count). jx RMSD = local RMSD (Å) of mapped
residues within 8 Å of the junction. Full columns (RMSD, chimera-normalised TM, sequence-free TM-align,
pTM of parents, AFDB pLDDT, wall time) are in
`domain_structure_metrics.tsv.gz`.

| fusion | Pfam | PR1 tier | flag | rules | dom pLDDT C | 5' | 3' | Δdom | jx mean C | jx min C | Δjx | half PAE C | 5' | 3' | pTM C | TM vs 5' | TM vs 3' | jx RMSD 5' | jx RMSD 3' | TM AFDB 5' | TM AFDB 3' |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| KRT82::KRT7 | PF00038 | likely_intact | fold_supported | S1_all_pass | 87.19 | 88.62 | 86.0 | -1.43 | 95.57 | 91.81 | -0.22 | 25.18 | 24.75 | 24.82 | 0.32 | 0.5136 | 0.6306 | 0.035 | 0.051 | 0.471 | 0.4966 |
| YWHAB::YWHAZ | PF00244 | likely_intact | fold_supported | S1_all_pass | 94.97 | 95.33 | 95.16 | -0.36 | 97.39 | 94.38 | -0.13 | 4.19 | 4.11 | 4.12 | 0.87 | 0.9965 | 0.9993 | 0.038 | 0.026 | 0.9944 | 0.9909 |
| PTN::MDK | PF01091 | likely_intact | fold_supported | S1_all_pass | 71.98 | 72.54 | 70.14 | -0.56 | 70.09 | 52.12 | -1.87 | 13.45 | 13.26 | 12.97 | 0.41 | 0.7228 | 0.8089 | 0.325 | 0.273 | 0.7165 | 0.7696 |
| CRIP1::CRIP2 | PF00412 | likely_intact | fold_supported | S1_all_pass | 83.54 | 83.12 | 81.48 | 0.42 | 81.24 | 71.06 | 0.64 | 5.8 | 5.98 | 6.23 | 0.5 | 0.993 | 0.9915 | 0.116 | 0.1 | 0.975 | 0.9723 |
| AFM::ALB | PF00273 | likely_intact | fold_supported | S1_all_pass | 91.86 | 92.26 | 95.48 | -3.62 | 87.86 | 82.62 | -7.46 | 4.31 | 4.18 | 3.09 | 0.85 | 0.9642 | 0.9741 | 0.592 | 0.599 | 0.9667 | 0.9651 |
| KLK3::PRSS8 | PF00089 | likely_intact | fold_supported | S1_all_pass | 92.27 | 94.36 | 93.43 | -2.09 | 92.03 | 85.56 | -3.19 | 3.52 | 2.83 | 3.15 | 0.8 | 0.9824 | 0.8611 | 0.381 | 3.731 | 0.9799 | 0.8586 |
| LPO::MPO | PF03098 | likely_intact | fold_supported | S1_all_pass | 92.93 | 92.43 | 95.09 | -2.16 | 93.32 | 86.94 | -2.58 | 4.88 | 5.06 | 4.19 | 0.93 | 0.9855 | 0.9749 | 0.162 | 0.218 | 0.9794 | 0.9745 |
| GPC3::GPC4 | PF01153 | uncertain | fold_doubtful | D3_junction_very_low;D5_halves_not_packed | 77.49 | 73.86 | 83.15 | -5.66 | 79.45 | 22.39 | -13.22 | 27.74 | 27.05 | 15.62 | 0.47 | 0.7981 | 0.6651 | 22.292 | 27.438 | 0.8063 | 0.6722 |
| VIM::KRT18 | PF00038 | uncertain | fold_supported | S1_all_pass | 84.0 | 85.91 | 86.35 | -2.35 | 91.17 | 78.75 | -0.47 | 25.1 | 25.25 | 25.65 | 0.32 | 0.5584 | 0.5418 | 0.113 | 0.081 | 0.5591 | 0.5175 |
| COL3A1::COL1A1 | PF01391 | uncertain | not_assessed | N1_parents_not_confident | 63.37 | 53.23 | 60.3 | 3.07 | 67.55 | 56.56 | 10.79 | 24.47 | 25.28 | 26.02 | 0.12 | 0.1794 | 0.2325 | 12.913 | 0.706 | 0.2896 | 0.2383 |
| SPNS2::SPNS3 | PF07690 | uncertain | fold_doubtful | D2_domain_plddt_drop;D3_junction_very_low;D4_junction_plddt_drop;D5_halves_not_packed | 74.14 | 87.34 | 88.13 | -13.99 | 70.56 | 24.19 | -21.33 | 19.14 | 7.71 | 6.7 | 0.61 | 0.5856 | 0.5683 | 8.122 | 8.146 | 0.5844 | 0.5809 |
| PIM3::ULK1 | PF00069 | uncertain | fold_doubtful | D3_junction_very_low;D4_junction_plddt_drop;D5_halves_not_packed | 88.08 | 95.05 | 90.19 | -6.97 | 78.93 | 39.12 | -16.37 | 24.56 | 3.61 | 5.87 | 0.5 | 0.6599 | 0.5463 | 10.906 | 11.045 | 0.6573 | 0.5425 |
| GNAL::GNAS | PF00503 | not_assessed | not_assessed | N3_afdb_comparison_missing | 92.99 | 92.16 | 89.42 | 0.83 | 88.27 | 58.59 | -1.39 | 5.63 | 5.98 | 6.57 | 0.86 | 0.9992 | 0.9513 | 0.045 | 0.085 | 0.9969 |  |
| COL1A1::COL18A1 | PF01391 | likely_broken | not_assessed | N1_parents_not_confident | 54.89 | 46.55 | 47.23 | 7.66 | 60.34 | 40.56 | 7.38 | 27.02 | 22.02 | 24.88 | 0.13 | 0.2138 | 0.1047 | 0.523 |  | 0.1751 | 0.0979 |
| COL1A1::COL12A1 | PF01391 | likely_broken | fold_doubtful | D1_domain_plddt_low;D3_junction_very_low;D4_junction_plddt_drop;D6_fold_differs | 67.38 | 68.92 | 77.32 | -9.94 | 61.08 | 40.5 | -17.6 | 14.25 | 24.93 | 24.41 | 0.23 | 0.255 | 0.2625 | 12.057 | 15.336 | 0.1878 | 0.2714 |
| GPC6::GPC5 | PF01153 | likely_broken | fold_doubtful | D1_domain_plddt_low;D2_domain_plddt_drop;D3_junction_very_low;D4_junction_plddt_drop;D5_halves_not_packed;D6_fold_differs | 54.62 | 82.42 | 76.81 | -27.8 | 52.93 | 45.41 | -37.89 | 24.88 | 13.51 | 26.47 | 0.23 | 0.0407 | 0.0491 | 23.975 | 22.235 | 0.0366 | 0.042 |
| GPC4::GPC3 | PF01153 | likely_broken | fold_doubtful | D1_domain_plddt_low;D2_domain_plddt_drop;D3_junction_very_low;D4_junction_plddt_drop;D5_halves_not_packed;D6_fold_differs | 59.02 | 83.15 | 73.86 | -24.13 | 52.98 | 45.16 | -39.8 | 23.16 | 15.7 | 27.13 | 0.28 | 0.0528 | 0.0425 |  | 28.1 | 0.0478 | 0.0403 |
| CMTM7::CMTM8 | PF01284 | likely_broken | not_assessed | N1_parents_not_confident | 67.98 | 58.48 | 64.33 | 3.65 | 71.02 | 64.75 | 5.55 | 14.17 | 14.22 | 15.0 | 0.37 | 0.2069 | 0.2375 | 4.408 | 5.349 | 0.1931 | 0.2421 |
| MYH10::MYH14 | PF00063 | likely_intact | fold_supported | S1_all_pass | 79.92 | 79.76 | 82.94 | -3.02 | 87.15 | 79.0 | -4.08 | 15.57 | 16.44 | 15.41 | 0.76 | 0.9293 | 0.9831 | 0.077 | 0.092 | 0.9592 | 0.9587 |
| KRT13::KRT5 | PF00038 | likely_intact | fold_supported | S1_all_pass | 84.85 | 84.68 | 85.03 | -0.18 | 82.67 | 55.19 | -2.88 | 23.14 | 23.21 | 22.22 | 0.3 | 0.6215 | 0.5395 | 0.107 | 0.086 | 0.5144 | 0.4246 |
| KRT13::KRT7 | PF00038 | likely_intact | fold_supported | S1_all_pass | 85.13 | 84.68 | 86.0 | -0.87 | 87.5 | 69.25 | -4.54 | 22.98 | 22.93 | 22.15 | 0.3 | 0.7278 | 0.498 | 0.114 | 0.118 | 0.4832 | 0.4534 |
| KRT8::KRT5 | PF00038 | likely_intact | fold_supported | S1_all_pass | 84.97 | 85.04 | 85.03 | -0.07 | 91.87 | 83.5 | -0.9 | 24.45 | 24.64 | 24.71 | 0.31 | 0.5171 | 0.7395 | 0.082 | 0.051 | 0.5141 | 0.4469 |
| KRT19::KRT7 | PF00038 | likely_intact | fold_supported | S1_all_pass | 83.39 | 85.68 | 86.0 | -2.61 | 79.33 | 67.75 | -5.5 | 24.63 | 24.76 | 24.02 | 0.31 | 0.497 | 0.5926 | 0.076 | 0.054 | 0.5752 | 0.4672 |
| KRT5::KRT8 | PF00038 | likely_intact | fold_supported | S1_all_pass | 85.6 | 85.03 | 85.04 | 0.56 | 94.62 | 92.12 | -0.88 | 24.91 | 24.96 | 24.86 | 0.31 | 0.7884 | 0.5592 | 0.036 | 0.11 | 0.4776 | 0.5233 |
| KRT5::KRT13 | PF00038 | likely_intact | fold_doubtful | D3_junction_very_low | 85.94 | 85.03 | 84.68 | 0.91 | 81.97 | 48.59 | -5.31 | 22.83 | 22.19 | 23.21 | 0.32 | 0.6794 | 0.5176 | 0.088 | 0.119 | 0.4267 | 0.5335 |
| FLI1::ERG | PF00178 | uncertain | fold_doubtful | D4_junction_plddt_drop | 89.6 | 94.14 | 93.64 | -4.54 | 77.65 | 50.59 | -16.69 | 5.76 | 2.48 | 2.68 | 0.72 | 0.9507 | 0.951 | 3.245 | 3.24 | 0.9482 | 0.9439 |
| STAT5B::STAT5A | PF02864 | uncertain | fold_doubtful | D3_junction_very_low;D4_junction_plddt_drop | 88.61 | 92.23 | 91.5 | -3.62 | 69.9 | 33.28 | -24.93 | 8.04 | 3.77 | 3.95 | 0.85 | 0.9747 | 0.9724 | 8.53 | 8.53 | 0.968 | 0.9681 |

## Named cases

### FLI1::ERG (ETS domain, 7-residue tandem repeat IYYYDKN) — `fold_doubtful` (D4)

The repeat duplicates the C-terminal end of the DNA-recognition helix H3
(`…KLSRALRYYYDKN` + `IYYYDKN` + `IMTKV…`). In the chimera model H3 is kept
as a helix from Y332 through K345 (mean pLDDT 87 vs 93 for the same
residues in FLI1), including both H3 arginines and the first YYY. The helix
then stops one residue early (N346 is helical in both parents), and the
7-residue copy is modelled as a **non-helical loop/bulge between H3 and the
β3 strand** with low confidence (pLDDT 51–70, mean 57). The following β-strands
resume in the parental register with pLDDT ~90. Domain-level TM is 0.95
against both parents, the domain pLDDT drop is −4.5, and inter-half PAE is
5.8 Å vs 2.5/2.7 Å. The junction window drops by −16.7 pLDDT (min 50.6),
which fires D4. Local RMSD around the junction is 3.2 Å over 8 residues.

Verdict: the insertion is **accommodated as an extruded loop after H3**
rather than lengthening or breaking H3, but the last turn of H3 is frayed
and lower-confidence. Whether DNA binding survives cannot be judged from a
monomer without DNA. **AlphaFold 3 with the ETS target DNA would be the
better test and is out of scope here.**

### STAT5B::STAT5A (DNA-binding domain, 16-residue tandem repeat) — `fold_doubtful` (D3, D4)

The repeat TKQPPQVLKTQTKFAA duplicates a β-strand segment of the STAT
DNA-binding domain. The model rebuilds the strand from the **second** copy
(KFAA 363–366, pLDDT 88–97, parental position). The end of the first copy
plus most of the insert (347–362) is pushed out as a ~16-residue loop of
very low confidence (mean pLDDT 45, min 33). The rest of the domain is
parent-like: domain pLDDT 88.6 vs 92.2/91.5, mapped TM 0.97 against both
parents and AFDB. Inter-half PAE is 8.0 Å vs 3.8/4.0 Å, below the D5
margin. The local RMSD of 8.5 Å around the junction comes from the
displaced first-copy residues.

Verdict: the core fold is kept, with a disordered excursion out of a strand
in the DNA-binding domain. Its effect on DNA binding (a dimer on DNA) is not
modelled.

### Glypicans (PF01153)

* **GPC6::GPC5** (17 samples, 4 cohorts) and **GPC4::GPC3** (3 samples) —
  `fold_doubtful` (D1–D6 all fire). Both are 466-residue deletions: the
  first ~30–40 residues of the glypican core are joined directly to the
  parents' C-terminal tail, which is disordered in the parents themselves
  (pLDDT ~25–35 there). The resulting ~80–95-residue "domain" is predicted
  at pLDDT 55/59 against 82/77 and 83/74 for the parents. Parent-normalised
  TM is 0.04–0.05: the glypican core is absent. This is a clear screen-out,
  consistent with PR 1's likely_broken. Note that PR 1's
  `aligned_deleted_columns` reads 0 for these records (see Limits); the
  screen's `pilot_indel_residues` is −466.
* **GPC3::GPC4** (2 samples) — `fold_doubtful` (D3, D5). This is the
  opposite case, a near-whole-domain tandem duplication (+457 residues). A
  complete GPC3 core (residues 14–524, mapped TM 0.80 against GPC3) is
  followed by an almost complete second core from GPC4 (inserted 525–982,
  mean pLDDT 81). The junction sits in GPC3's C-terminal stalk, which is
  low-confidence in GPC3 itself (pLDDT ~30–40). D5 fires because the two
  "halves" are two separate domains (inter-half PAE 27.7 Å) while the GPC4
  parent baseline is one domain (15.6 Å). Read this flag as "two domains
  in tandem", not as misfolding of either copy.

### LPO::MPO (peroxidase, PF03098) — `fold_supported`

A clean reconstitution, with the junction at G369 in a loop between a strand
and a helix. Secondary structure around the junction is identical to both
parents. Domain pLDDT is 92.9 vs 92.4/95.1, junction-window min 86.9,
inter-half PAE 4.9 vs 5.1/4.2 Å, and mapped TM is 0.99/0.97 (0.98/0.97
against AFDB). LPO and MPO are close paralogues, so this is the "confident
fold, weak evidence" situation: nothing argues against folding. Heme
binding, its covalent links and MPO's proteolytic processing are not
modelled.

### PIM3::ULK1 (protein kinase, PF00069) — `fold_doubtful` (D3, D4, D5)

PIM3 and ULK1 are not close paralogues, so this flag is more informative
than most. The junction falls in PIM3's activation segment (after
`…YTDFD`). The ULK1 part restarts 77 alignment columns earlier, so the
chimera carries a **second copy of the catalytic segment**: ULK1's αE
region, HRD catalytic loop (HRDLKPQN), DFG motif and activation loop.

The ULK1 copy is fairly confident (mean pLDDT 81), but the PIM3 residues
just before the junction drop from a mean of 94 to 70 (198–205: pLDDT 39–68, below 50 from residue 202).
**Inter-half PAE is 24.6 Å against 3.6/5.9 Å in the parents**: the model
does not place the PIM3 N-terminal part and the ULK1 C-terminal part
confidently relative to each other. Mapped TM is 0.66 vs PIM3 and 0.55 vs
ULK1. Verdict: a doubtful chimeric kinase domain; the halves do not pack
as one kinase fold in this model.

## Other observations

* **Keratin / vimentin rods (PF00038, 8 records).** 7 are `fold_supported`
  and KRT5::KRT13 is `fold_doubtful` (D3: one junction-window residue at
  pLDDT 48.6). For a monomeric coiled-coil rod, whole-domain TM-scores sit
  around 0.5 even between our parent model and AFDB (0.45–0.61). Inter-half
  PAE is ~25 Å in parents and chimeras alike. The informative values are the
  junction ones: local RMSD 0.04–0.12 Å and junction pLDDT within 5.5 of the
  better parent.
* **COL1A1::COL12A1** is `fold_doubtful` rather than `not_assessed` only
  because one parent window (COL12A1) reaches pLDDT 77. Our collagen parent
  models agree with AFDB at TM ~0.3, so D6 is not meaningful for collagen.
  Treat this flag as weak; PR 1 already has it likely_broken.
* **SPNS2::SPNS3** (+186-residue duplication in an MFS transporter) is
  `fold_doubtful` (D2–D5): domain pLDDT −14 vs the better parent,
  inter-half PAE 19.1 vs 6.7–7.7 Å.
* **GNAL::GNAS** is `not_assessed` (N3). Its ColabFold measures look
  parent-like (TM 0.999/0.951 against the two parent predictions, domain
  pLDDT 93.0 vs 92.2/89.4). However, the UniProt entry PR 1 mapped for GNAS
  is P84996 (ALEX, an alternative-reading-frame product), which matches the
  parent window at 0.3 % identity. The AFDB 3' comparison was therefore
  refused, and a record lacking a required comparison cannot be
  `fold_supported`. It was `fold_supported` in the first version of this
  run.
* As a check of the reduced CPU settings, our parent-window models agree
  with AFDB at mapped TM ≥ 0.84 for every globular parent; the exceptions
  are the rods (above), collagen, CMTM7 (0.72) and PTN/MDK (0.77/0.76).

## Limits

* One model (model_1_ptm), one recycle, reduced MSA depth, no templates,
  one seed. On the FLI1::ERG benchmark the reduced MSA depth reproduced
  the default-depth model (TM 0.998, 0.16 Å), and recycles beyond the first
  raised pLDDT by < 1. They are
  still not the 5-model consensus a GPU run would give. The settings are one
  config file (`--config`) for a later GPU rerun, and the cached MSAs are
  reused.
* Every prediction is a monomeric domain window ±30 residues. Partners
  (keratin heterodimers, STAT dimers, DNA, heme, the rest of each fusion
  protein) are absent.
* The flag thresholds are documented heuristics, not calibrated. A
  confident fold for a near-identical paralogue chimera is weak evidence;
  a doubtful one is a reason to deprioritise, not proof of misfolding.
* Upstream (PR 1, not changed here): `place_chimera` reports 0 deleted
  alignment columns when the junction has a hybrid-codon residue, because
  the hybrid takes the 3' column next to the 3' junction residue. This
  affects the reported `aligned_deleted_columns` of e.g. the 466-residue
  glypican deletions. The structure measures here do not depend on that
  count (halves, mapping and TM normalisation use positions directly), but
  PR 1 consumers should use `pilot_indel_residues` for indel size.
