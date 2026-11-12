# Structure prediction of shortlisted chimeric domains

`cfh.domain_structure` (driver: `scripts/domain_structure.py`) is PR 2 of the
domain-function stack. It takes the 27-record structure shortlist written by
the sequence step (`docs/domain_function_seq.md`) and asks whether a
structure predictor builds the chimeric domain as well as it builds the two
parent domains, under identical settings. PR 1's tiers are not changed; this
step adds separate measures and a separate flag.

**The flag is a screen-out, not a confirmation.** For near-identical paralogue
chimeras (most of the shortlist) an AlphaFold-family model sees essentially
the parents' MSA and will usually predict a confident, parent-like fold
whatever the junction does. A confident chimera fold is therefore *weak*
evidence that the chimera folds; a failed or clearly worse fold is
informative. Nothing here is a calibrated predictor of folding or activity.

## Inputs and jobs

* `shortlist_for_structure.tsv` and `chimera_records.json.gz`
  (`domain_function_seq/2`) from `runs/domain_function_seq_2026-11-11/`.
  Each shortlisted `record_id` must exist in the records file with the same
  `sequence_id` (else the run stops).
* Per record three sequences are predicted: the chimeric domain with the same
  ±30-residue window the sequence step used (`chimeric_domain_flanked`) and
  each parent's Pfam domain ±30 (`ParentDomain.window()`).
* Jobs are keyed by sequence (`job_<sha256[:16]>`), so a parent window shared
  by several fusions (e.g. KRT5, KRT7, COL1A1) or a chimera shared by several
  records is predicted once. The 2026-11-11 shortlist gives 71 distinct
  sequences (27 chimeras, 44 parent windows; 22,812 residues).

## Predictor and settings

ColabFold 1.6.3 (`colabfold_batch`, AlphaFold2 weights
`alphafold_params_2021-07-14`, `alphafold2_ptm`), CPU JAX on an Apple M3 (no
GPU). Settings are one `PredictorConfig` (JSON override with `--config`); the
settings that change a prediction are hashed into a `settings_id`, and cached
predictions live under it, so changing a setting never reuses old results.

| setting | value | why |
|---|---|---|
| models | 1 (`model_1_ptm`) | pTM head needed for pTM; 5 models would be 5x the CPU time |
| recycles | 1 (2 passes) | on the FLI1::ERG benchmark pLDDT rose 73.9 -> 79.9 at recycle 1 and only to 80.6 by recycle 3 |
| MSA depth (`--max-msa`) | 256:512 | default 512:5120 took 10.5 min for 150 aa; 256:512 gave the same domain (TM 0.998, RMSD 0.16 A vs the default) in 6.3 min with 3 recycles |
| templates | off | ColabFold default; avoids copying a parent's PDB structure into the chimera |
| seed | 0 (1 seed) | |
| relax | off | pLDDT/PAE/CA geometry do not need Amber |
| MSA | `mmseqs2_uniref_env` on the public server `https://api.colabfold.com` | user's choice |

The same settings are used for chimeras and parents, so every comparison is
like-for-like. Wall time per prediction and the total are in
`run_manifest.json` (`wall_time`).

**Moving to a GPU host.** Install ColabFold there and pass a config such as
`{"colabfold_command": ["/path/to/colabfold_batch"], "msa_python":
"/path/to/python", "data_dir": "/path/to/weights", "num_models": 5,
"num_recycle": 3, "max_msa": "512:5120"}`. The MSAs are cached by sequence
and reused (they do not depend on the prediction settings).

**Install (outside the git tree).** `$CFH_COLABFOLD_HOME` (default
`~/.cache/cfh_colabfold`): `env/` (uv venv, `pip install
"colabfold[alphafold]"`, ~0.5 GB) and `weights/params/` (only the five
`params_model_*_ptm.npz` streamed out of the 3.7 GB tar, ~1.9 GB, plus an
empty `download_finished.txt` marker). The work directory
`runs/.domain_structure_work/` (MSAs, per-job checkpoints, logs) is
gitignored. The run stops cleanly (`ABORT`, exit status 2) if free disk falls
below `min_free_disk_gb` (5 GB), checked before each MSA ticket and each
prediction.

### MSA server etiquette

* All uncached sequences are sent in tickets of `msa_batch_size` (8)
  sequences, one ticket at a time, with a 20 s pause between tickets; the
  server's `RATELIMIT` / queue answers are handled by ColabFold's own
  back-off (`run_mmseqs2`).
* The user agent is `cfh-domain-structure/1` (no personal contact details).
* Every returned a3m is checked to start with its query and cached as
  `msa_<server hash>/<sequence hash>.a3m`; a cached a3m is never re-queried.
  The MSA worker (`msa_worker.py`) runs inside the ColabFold env.

### Fail closed and checkpoints

A prediction is cached (`predictions/<settings_id>/<job>/result.json`, model
`model.pdb.gz`, scores `scores.json.gz`) only after validation: exit status 0,
exactly one rank-1 model and score file, pLDDT and PAE of the query's length,
finite pLDDT/PAE/pTM, and a single-chain model whose sequence is the query.
Cached results are re-validated on load; a corrupt checkpoint is discarded.
Any failure -- MSA server error, missing MSA, ColabFold exit status, timeout
(`job_timeout_seconds`), missing model file, parse error -- is written to
`failure.json` with its reason, never read back as a result, and retried on
the next run. A record with any failed prediction, or a failed mapping, is
`not_assessed` with the reason in `structure_flag_missing`.

An interrupted `predict` resumes at job granularity: finished jobs are
skipped, the interrupted job is rerun.

## Definitions

Positions are each protein's own 1-based residue numbers.

* **Domain** -- the Pfam span (chimera: `chimeric_domain_span`; parent: its
  own Pfam span), not the ±30 flank.
* **Junction bond** -- between chimera residues `J` and `J+1`
  (`J = junction_index`, the last 5'-derived residue, hybrid included); in
  the 5' parent between `F` and `F+1`, in the 3' parent between `T-1` and `T`.
* **Junction window** -- 10 residues each side of the bond (`J-9..J+10`),
  widened to 10 residues beyond any junction-inserted residue
  (`junction_window` column); parents use the same ±10 around their own bond.
* **Halves** -- chimera domain residues from the 5' parent vs from the 3'
  parent (hybrid residue in neither). A parent's halves are its own domain
  residues before / after its junction bond, so each parent gives the
  like-for-like baseline (`*_half_sizes`).
* **Residue mapping** -- chimera residue -> parent residue through PR 1's
  parent-window alignment (`align_parents` / `place_chimera`: BLOSUM62,
  global, gap -11/-1, free end gaps). Junction-inserted chimera residues and
  columns where the parent has a gap are unmapped. Comparisons use mapped
  domain residues only (`vs_*_n_mapped`); fewer than 10 is a mapping failure.

## Measures (`domain_structure_metrics.tsv.gz`, one row per record)

Each measure is reported on its own; there is no combined score.

Per prediction (`chimera_*`, `parent5_*`, `parent3_*`):

* `domain_plddt` -- mean pLDDT over the domain;
* `junction_plddt_mean`, `junction_plddt_min` -- over the junction window;
* `domain_pae` -- mean PAE over ordered domain residue pairs `i != j`;
* `interhalf_pae` -- mean PAE between the two halves, symmetrised
  (`(mean PAE[A,B] + mean PAE[B,A]) / 2`). This is the most direct measure of
  whether the model places the 5'- and 3'-derived parts confidently relative
  to each other;
* `ptm` -- pTM of the whole predicted window (ColabFold reports pTM for the
  full sequence, flanks included);
* `wall_seconds`, `prediction_status`, `job`.

Deltas: `delta_domain_plddt_vs_better_parent`,
`delta_junction_plddt_mean_vs_better_parent` = chimera minus the higher of
the two parents.

Chimera vs each ColabFold parent model (`vs_5p_*`, `vs_3p_*`), over mapped
domain residues:

* `tm` -- TM-score with the fixed mapping, normalised by the chimera domain
  length (own implementation of the TM-score search: seed fragments L, L/2,
  ... >= 4, iterative re-superposition; Zhang & Skolnick 2004);
* `rmsd` -- RMSD of all mapped CA pairs after least-squares superposition;
* `local_rmsd`, `n_local` -- RMSD after local superposition of the mapped
  residues whose CA is within 8 A of a junction residue (`J`, `J+1`, and any
  junction-inserted residue) in the chimera model;
* `tmalign_seqfree` -- sequence-independent TM-align (`tmtools`) of the two
  domain regions, normalised by the chimera domain; secondary, register-free.

Chimera vs the AlphaFold DB parent models used in PR 1 (`vs_afdb_5p_*`,
`vs_afdb_3p_*`): same measures, after mapping the parent window onto the AFDB
sequence (identity >= 0.9 required, `afdb_*_identity`); `afdb_*_model` names
the entry/version. As a check of the CPU settings, `parent5_vs_afdb_tm` /
`parent3_vs_afdb_tm` compare our parent-window model with AFDB over the
parent domain, and `afdb_*_domain_plddt` gives AFDB's own domain pLDDT.

## HEURISTIC structure flag

`structure_flag`, with `structure_flag_rules` (rules fired) and
`structure_flag_missing` (reasons). Rules are also in `run_manifest.json`
(`flag_rules`).

| rule | flag | condition |
|---|---|---|
| N0_prediction_missing | not_assessed | a chimera or parent prediction failed / is missing (reason recorded) |
| N1_parents_not_confident | not_assessed | neither parent's domain pLDDT >= 70: no confident reference fold |
| D1_domain_plddt_low | fold_doubtful | chimera domain pLDDT < 70 (a parent is >= 70) |
| D2_domain_plddt_drop | fold_doubtful | chimera minus better parent domain pLDDT <= -10 |
| D3_junction_very_low | fold_doubtful | chimera junction-window min pLDDT < 50 while the better parent's junction-window mean >= 70 |
| D4_junction_plddt_drop | fold_doubtful | chimera minus better parent junction-window mean pLDDT <= -15 |
| D5_halves_not_packed | fold_doubtful | chimera inter-half PAE >= 10 A and >= 5 A above the lower parent inter-half PAE |
| D6_fold_differs | fold_doubtful | mapped TM-score < 0.5 against both parent predictions |
| N2_metric_missing | not_assessed | no doubt rule fired, but a measure needed for fold_supported is missing (named) |
| S1_all_pass | fold_supported | everything assessed, no doubt rule fired |

Doubt rules fire on whatever is available; `fold_supported` needs every
input (all three predictions, both mappings, the pLDDT deltas, the junction
minimum, and inter-half PAE for the chimera and at least one parent).

Thresholds, justified:

* **pLDDT 70 / 50.** AlphaFold's own confidence bands: >= 70 "confident"
  backbone, < 50 "very low", typically disordered (Jumper et al. 2021,
  Nature 596:583; Tunyasuvunakool et al. 2021, Nature 596:590).
* **TM-score 0.5.** Above 0.5 two structures share the same fold with high
  probability; below it they mostly do not (Xu & Zhang 2010, Bioinformatics
  26:889).
* **-10 domain / -15 junction pLDDT.** Roughly one confidence band; the
  junction window is short and noisier, so its threshold is wider.
* **Inter-half PAE >= 10 A and +5 A over the parent.** PAE below ~5 A means
  the relative placement is confidently predicted; 10 A is a position error
  of about two helix turns. The +5 A margin over the parents' own baseline
  keeps multi-lobed domains, whose parents already have high inter-half PAE,
  from firing.

These are documented heuristics, not calibrated values. Collagen
triple-helix domains (PF01391) are not foldable as monomers, so their parents
are typically low-confidence and the records end up `not_assessed`
(N1) -- by design, not as a failure.

## Outputs (`runs/domain_structure_<date>/`)

* `domain_structure_metrics.tsv.gz` -- one row per shortlisted record (ids,
  PR 1 tier unchanged as `pr1_tier`, all measures, flag);
* `models/<job>.pdb.gz` -- the top-ranked (only) model per distinct sequence;
  `models/<job>.scores.json.gz` -- pLDDT (2 decimals), PAE (1 decimal; the
  metrics use full precision), pTM. MSAs and other models are not committed;
* `flag_counts.json` -- flags overall, per cohort, per PR 1 tier, rules fired;
* `run_manifest.json` -- code commit, ColabFold/AlphaFold/JAX versions, weight
  source and sha256, settings and `settings_id`, MSA server, input sha256,
  per-job and total wall time, job index (which records use each model);
* `REPORT.md` -- results, named cases and caveats.

## Out of scope

* For ETS-domain chimeras (FLI1::ERG) the relevant test is DNA binding; an
  AlphaFold 3 model with the target DNA would be the better test and is not
  done here.
* Multimers (e.g. STAT dimers, keratin heterodimers) are not modelled; every
  prediction is a monomeric domain window.
