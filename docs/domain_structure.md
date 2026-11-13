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
settings that change a prediction -- including the MSA source (`msa_mode`,
`host_url`) -- are hashed into a `settings_id`, and cached predictions live
under it, so changing a setting or MSA server never reuses old results.

*Checkpoint provenance.* Each checkpoint's `result.json` records
`msa_source` (`msa_mode`, `host_url`, MSA cache name), `a3m_sha256` and
`msa_provenance`:

* `recorded` -- written by `predict` at prediction time (sha256 of the exact
  a3m handed to ColabFold);
* `inferred` -- written afterwards by the one-off
  `scripts/domain_structure_inferred_provenance.py`
  (`cfh.domain_structure.provenance`), together with
  `msa_provenance_evidence`; its `a3m_sha256` is the a3m in the cache when
  the record was written, not at prediction time.

A checkpoint is used only if it has such a record, the recorded or inferred
`msa_source` equals the current config's, and an inferred record carries its
evidence. Anything else is kept on disk but never used, and there is no
automatic migration between settings ids (checkpoints under another id are
listed and ignored). `msa_provenance` (and `*_msa_provenance` per prediction)
is a column of the metrics TSV, and the manifest counts both kinds.

*Inferred provenance of the 2026-11-13 run.* Its 71 predictions were made at
3da242a under `set_e20d3cb4e3e2`, before the key included the MSA server and
before checkpoints recorded their source, so no per-job record exists. The
user chose to keep them with inferred provenance rather than re-predict
(~12 h CPU). The evidence (`EVIDENCE_2026_11_12`) is run-level:

1. a single `predict` launch without `--config` (session transcript), so the
   default config applied, whose `host_url` at 3da242a was
   `https://api.colabfold.com`;
2. a single invocation in `predict.log` (2026-11-12 16:23:38 to 2026-11-13
   05:03:03);
3. a single MSA cache, `msa_c3444ece` = sha256(`mmseqs2_uniref_env|
   https://api.colabfold.com`)[:8], filled 16:23:37-17:16:51, before the
   first prediction (17:17);
4. the FLI1::ERG chimera a3m was seeded from the benchmark fetch, whose
   request (`bench/req.json`) names `https://api.colabfold.com`.

The step is never run by `predict` or `analyze`. It refuses, and changes
nothing, unless the config's MSA source is exactly the evidence's (it cannot
relabel to another server), the work directory holds only that one MSA
cache, the legacy settings id equals the config's settings without
`host_url`, and every job re-validates with its query's a3m in that cache.
It then writes the inferred record and moves each checkpoint to the config's
settings id.

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

Committed run (`runs/domain_structure_2026-11-13/`): 71 sequences, 9 MSA
tickets in ~53 min, then 11.8 h of prediction wall time on the M3 CPU
(23 s for 113 residues to 82 min for the 1,027-residue GPC3::GPC4 window;
time grows roughly with length squared), all 71 validated, none failed.

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
Loaded checkpoints go through the same check (`check_scores`: numeric and
finite pLDDT[L], PAE[L][L] and pTM, model sequence = query); a corrupt
checkpoint is discarded and re-predicted, never used.
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
  columns where the parent has a gap are unmapped. A pair is kept only when
  the chimera residue is in the chimeric domain **and** the parent residue is
  in that parent's Pfam domain (flank residues never count), for both the
  ColabFold-parent and the AFDB comparisons (`vs_*_n_mapped`). Fewer than 10
  pairs, more pairs than the parent or chimera domain length, a mapping that
  is not one-to-one, or a TM-score outside [0, 1] is a mapping failure.

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

* `tm` -- TM-score with the fixed mapping, normalised by the **parent**
  domain length (own implementation of the TM-score search: seed fragments
  L, L/2, ... >= 4, iterative re-superposition; Zhang & Skolnick 2004). The
  parent is the reference, as a native structure is when a model is scored:
  chimera residues with no parent counterpart (a junction insertion) do not
  lower `tm`, parent domain residues missing from the chimera (a deletion) do.
  A first version normalised by the chimera domain, which capped a
  whole-domain duplication such as GPC3::GPC4 (985-residue chimeric domain,
  ~520 mapped) near 0.5 however well the mapped residues matched;
* `tm_chim_norm` -- the same superposition search normalised by the chimera
  domain length (what fraction of the chimeric domain sits like the parent);
* `rmsd` -- RMSD of all mapped CA pairs after least-squares superposition;
* `local_rmsd`, `n_local` -- RMSD after local superposition of the mapped
  residues whose CA is within 8 A of a junction residue (`J`, `J+1`, and any
  junction-inserted residue) in the chimera model;
* `tmalign_seqfree` -- sequence-independent TM-align (`tmtools`) of the two
  domain regions, normalised by the parent domain; secondary, register-free.

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
| D6_fold_differs | fold_doubtful | mapped TM-score (normalised by the parent domain) < 0.5 against both parent predictions |
| N2_metric_missing | not_assessed | no doubt rule fired, but a measure needed for fold_supported is missing (named); a NaN/inf value counts as missing |
| N3_afdb_comparison_missing | not_assessed | no doubt rule fired, but the comparison with either parent's AlphaFold DB model was refused or failed (no model, sequence mismatch, identity < 0.9, mapping failure) |
| S1_all_pass | fold_supported | everything assessed, no doubt rule fired |

Doubt rules fire on whatever is available (finite values only);
`fold_supported` needs every input: all three predictions, both mapped
ColabFold-parent comparisons, **both AFDB parent comparisons**, the pLDDT
deltas, the junction minimum, and inter-half PAE for the chimera and at least
one parent. Every rule input must be a finite number; NaN or inf is missing.

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

Known blind spots, kept visible rather than tuned away after seeing the data:

* **Rod domains (PF00038, keratin/vimentin).** A monomeric ~310-residue
  coiled-coil rod is one long, gently curved helix; its curvature is
  arbitrary without the partner chain, so whole-domain TM-scores sit around
  0.5 even between our parent model and AFDB (`parent*_vs_afdb_tm`), and
  inter-half PAE is ~25 A in parents and chimeras alike. For these records
  the junction measures (pLDDT, local RMSD) carry the information; D6 is
  close to its threshold by construction.
* **Whole-domain tandem duplications** (GPC3::GPC4). The "3' half" is a
  complete second domain, so D5 compares two separate units with a parent
  baseline that is a single domain; D5 firing there means "two domains, not
  one", not necessarily misfolding.
* **Collagen** (below).

These are documented heuristics, not calibrated values. Collagen
triple-helix domains (PF01391) are not foldable as monomers, so their parents
are typically low-confidence and the records end up `not_assessed`
(N1) -- by design, not as a failure.

## Outputs (`runs/domain_structure_<date>/`)

* `domain_structure_metrics.tsv.gz` -- one row per shortlisted record (ids,
  PR 1 tier unchanged as `pr1_tier`, the screen's `pilot_indel_residues` and
  `hybrid_residue`, all measures, flag). `aligned_deleted_columns` /
  `n_deleted_columns` come from PR 1's `place_chimera`, which counts 0
  deleted columns when the junction has a hybrid-codon residue (the hybrid
  takes the 3' column next to the 3' junction residue), e.g. 0 for the
  466-residue glypican deletions; use `pilot_indel_residues` for indel size;
* `models/<job>.pdb.gz` -- the top-ranked (only) model per distinct sequence;
  `models/<job>.scores.json.gz` -- pLDDT (2 decimals), PAE
  (`--pae-decimals`, default 1; the metrics always use the full-precision
  cached scores), pTM. `--export-fusions A::B,C::D` writes only those
  records' chimera and parent models (the committed 2026-11-13 run exports
  the named cases, PAE in whole A, ~2.5 MB); every model stays in the work
  cache and is listed in the manifest with `exported`. MSAs are never
  committed;
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
