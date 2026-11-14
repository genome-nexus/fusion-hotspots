# Domain stability and protein-language-model scores

Run date: 2026-11-13. Stacked on PR #160 (`polly/domain-structure`, base `09207fc`).

ESM-2 assessed **224 of 3,039 records (7.37%)**: all 45 priority records and 179 of the 2,994 random-order records. The other 2,815 records are **not assessed**, with blank PLM metrics. These are domain/origin records, not independent biological samples. No folding-energy or ΔΔG estimate was produced; these are sequence-plausibility rankings, not proof of stability or function.

## Coverage

| group | records | assessed | not assessed |
|---|---:|---:|---:|
| PR #155 shortlist | 27 | 27 | 0 |
| Additional named priority records | 13 | 13 | 0 |
| Additional PCAWG DNA events from #161 | 5 | 5 | 0 |
| Random-order remainder | 2,994 | 179 | 2,815 |
| **Total** | **3,039** | **224** | **2,815** |

The #155 collection contributes 219 assessed / 3,034 records; the five additional #161 DNA records are all assessed. All 27 shortlisted records retain their #160 flags: 14 `fold_supported`, nine `fold_doubtful`, and four `not_assessed`. The rank reference contains the 219 assessed #155 records; DNA records do not shape that reference. One assessed record (GPC3::GPC4) used windowed context; four planned records require it.

The not-assessed rows lack required cached positions: 2,809 first fail at the chimera, five at parent5 and one at parent3. No missing result is assigned a neutral or favourable score. Shared cached sequences can complete records outside the visited prefixes, so assessed coverage exceeds the number of visits made by these two invocations.

## Scoring completion and provenance

Both inherited drivers (PIDs 26648, 37339) and model workers (26694, 37408) were confirmed absent with `ps` at 21:30 local, before final analysis. No additional scoring was started and no scoring process was killed. The existing PLM environment was reused.

* MPS: the log records the deadline at 21:24:46 after 171 record visits. The old 60-second shutdown wait then raised `TimeoutExpired`; the worker subsequently exited. Valid per-sequence cache files were already durable. Its stored `score_status.json` is from 17:09 and is **not** this invocation’s final status. The driver now waits for graceful model teardown; regression coverage exercises delayed shutdown.
* CPU: deadline after 27 reverse-order visits at 21:29:09; completed at 21:29:10, zero failures, 13,216.9 seconds including the final in-flight record.

The [archived scoring logs](scoring_logs/) preserve both outcomes and are hashed in [run_manifest.json](run_manifest.json). Traceback source excerpts can differ from the loaded code because this worktree was updated while the inherited driver was running. The manifest records analysis source commit `2e96f9452487d171d53e525ddbfbfdc3904bdf70`, model/weight identity, input hashes and the final stored status files. Record visits include cache hits and do not equal new independent predictions.

## Named cases

All named cases and all five PCAWG DNA events are assessed. Values below are nats per residue; Δ is chimera mean domain PLL minus the specified parent mean. Junction dip is the chimera junction-window mean minus its own domain mean. These are separate measures; there is no combined score. Ranges span distinct domain/origin records. Full record IDs, parent PLLs, paired-context shifts, parent-specific dips and disagreements appear in [report_tables.md](report_tables.md) and [plm_metrics.tsv.gz](plm_metrics.tsv.gz).

| fusion | source | records | Δ vs 5′ | Δ vs 3′ | junction dip | excess dip vs 5′ | excess dip vs 3′ | PLM rank |
|---|---|---:|---:|---:|---:|---:|---:|---|
| FLI1::ERG | #155 | 1 | +0.0796 | +0.0113 | +0.0011 | +0.0740 | -0.0768 | not_bottom_decile |
| STAT5B::STAT5A | #155 | 1 | +0.0418 | -0.0020 | +0.1669 | -0.1139 | -0.0065 | not_bottom_decile |
| GPC6::GPC5 | #155 | 4 | -0.0363 to +0.0587 | -0.0835 to +0.0116 | -0.1758 to +0.0051 | -0.2413 to +0.0760 | -0.0962 to +0.1961 | not_bottom_decile |
| GPC4::GPC3 | #155 | 5 | -0.0054 to +0.0823 | +0.0088 to +0.0965 | -0.0425 to +0.1183 | -0.0848 to +0.1874 | -0.0934 to +0.0915 | not_bottom_decile |
| GPC3::GPC4 | #155 | 3 | -0.0253 to -0.0011 | -0.0395 to -0.0153 | +0.0163 to +0.0882 | -0.1141 to +0.0488 | +0.0090 to +0.1862 | not_bottom_decile |
| LPO::MPO | #155 | 2 | +0.0194 to +0.0345 | +0.0200 to +0.0351 | +0.0329 to +0.0366 | -0.0034 to +0.0037 | +0.0086 to +0.1276 | not_bottom_decile |
| PIM3::ULK1 | #155 | 1 | -0.0198 | +0.0458 | -0.1572 | +0.0835 | -0.2618 | not_bottom_decile |
| GNA14::GNAQ | #155 | 1 | -0.0403 | -0.1099 | -0.0226 | +0.1144 | -0.0883 | not_bottom_decile |
| ABCA6::ABCA10 | #155 | 2 | -0.0280 to -0.0113 | -0.0409 to +0.0107 | -0.1895 to -0.1791 | +0.0052 to +0.0072 | -0.0720 to -0.0026 | not_bottom_decile |
| MAP3K3::STRADA | PCAWG DNA | 1 | +0.0581 | +0.0146 | -0.0335 | +0.1732 | +0.0052 | not_bottom_decile |
| GNA14::GNAQ | PCAWG DNA | 1 | -0.0331 | -0.0961 | -0.0342 | +0.1071 | -0.1021 | not_bottom_decile |
| ABCA6::ABCA10 | PCAWG DNA | 1 | -0.1138 | -0.1697 | +0.1022 | +0.3400 | -0.1801 | not_bottom_decile |
| AKR1C2::AKR1C3 | PCAWG DNA | 1 | -0.0124 | +0.0069 | -0.0130 | -0.0157 | +0.0035 | not_bottom_decile |
| KIAA1919::SLC16A10 | PCAWG DNA | 1 | -0.2372 | -0.2328 | -0.2919 | -0.5930 | +0.0608 | bottom_decile |

**FLI1::ERG and STAT5B::STAT5A:** domain likelihood is near or above the parents, and neither ranks in the bottom decile. Both remain #155 `uncertain` and #160 `fold_doubtful`; PLM plausibility does not remove the structural concern.

**Glypican pairs:** all 12 records across GPC6::GPC5, GPC4::GPC3 and GPC3::GPC4 are outside the PLM bottom decile. This includes deletion records called `likely_broken` by #155 and three shortlisted glypican records with #160 `fold_doubtful`. Near-paralog sequence plausibility can coexist with an unfavourable sequence tier or structural flag.

**LPO::MPO:** both clean-reconstitution records have slightly higher domain PLL than either parent and are outside the bottom decile. Both are #155 `likely_intact`; the shortlisted record is #160 `fold_supported`, while the other has no #160 assessment.

**PIM3::ULK1:** Δ vs 5′ = −0.0198, Δ vs 3′ = +0.0458; junction dip = −0.1572 and excess vs 3′ = −0.2618. It remains #155 `uncertain` / #160 `fold_doubtful`, despite not ranking in the PLM bottom decile.

**PCAWG DNA events:** MAP3K3::STRADA is slightly above both parents and #155-method `uncertain`. ABCA6::ABCA10, AKR1C2::AKR1C3 and GNA14::GNAQ have `likely_broken` sequence tiers but are outside the PLM bottom decile. KIAA1919::SLC16A10 is `likely_broken` and PLM bottom decile: Δ = −0.2372/−0.2328, junction dip = −0.2919, worst parent-adjusted dip excess = −0.5930. None of these five additional DNA records has a #160 structure assessment. Their `pr155_tier` column copies the same sequence-tier method’s results from the #161 metrics input; these records were not in the original #155 collection. DNA candidates’ expression and biological function are not established by PLM scores.

## Disagreements

**39 records** have at least one descriptive disagreement. Category counts overlap; the record-level table is in [report_tables.md](report_tables.md#disagreements). No not-assessed PLM result takes part in a PLM disagreement.

| comparison | records |
|---|---:|
| `plm_low_but_fold_supported` | 2 |
| `plm_low_but_tier_likely_intact` | 20 |
| `plm_not_low_but_fold_doubtful` | 9 |
| `plm_not_low_but_tier_likely_broken` | 13 |
| `tier_likely_intact_but_fold_doubtful` | 1 |

Among the #155 records, 20 `likely_intact` cases are PLM bottom decile and ten `likely_broken` cases are outside it. The latter category rises to 13 when the three additional DNA cases are included. Two #160 `fold_supported` records are PLM bottom decile; all nine `fold_doubtful` records are outside it. This disagreement is a reason to retain separate measures, not to re-tier the candidates.

## Method

**Model.** ESM-2 `esm2_t12_35M_UR50D` (35 M parameters, 12 layers; fair-esm 2.0.0, torch
2.5.1, Apple-silicon MPS, fp32). Weights sha256 `7f21e80e61d1…cc55e` (in
`run_manifest.json`). Chosen over `esm2_t30_150M_UR50D` because on this 16 GB M3 the 150M model
was ~4× slower per masked pseudo-log-likelihood (26 s vs 6.5 s for a 243-residue PLL on an idle
machine) and the full set is ~927 k masked positions; fp16 on MPS was slower than fp32. The existing PLM env pins torch 2.5.1 for MPS compatibility on macOS 14.5. Model, env and torch
cache together use ~1.0 GB under `~/.cache/cfh_plm` (well under the 5 GB optional-tool budget). The `cfh`
package never imports torch: `src/cfh/domain_plm/worker.py` runs under the PLM env's Python and
talks JSON lines to the driver, so CI needs no ML dependencies.

**Score.** Masked pseudo-log-likelihood (PLL): for each residue *i*, mask *i* alone and take
`log p(true residue | rest)` (natural log, nats). Three sequences per record, the same windows
PR #160 predicted: chimeric domain ±30 residues, and each parent's Pfam domain ±30. Flanks are
context only. Sequences over the chosen 1,022-residue context limit (4 records: GPC3::GPC4 and three
MYH fusions) are scored per position in the 1,022-residue window centred on that position; those
rows carry `plm_note = <role>_windowed_context_<n>aa`.

**Length-robust comparisons.** Indels change the domain length, so summed PLLs are never
compared. All measures are per residue:

| column | definition |
|---|---|
| `chimera_domain_pll`, `parent5_domain_pll`, `parent3_domain_pll` | mean PLL over the Pfam domain residues |
| `delta_vs_5p`, `delta_vs_3p` | chimera mean minus parent mean |
| `delta_vs_better_parent` | chimera mean minus the higher parent mean |
| `context_shift_5p`, `context_shift_3p` | **paired** mean of (chimera PLL − parent PLL) over the *same* residues (each chimera residue copied from that parent, scored in both contexts). Inserted/hybrid residues have no partner and deleted residues are absent, so length drops out entirely. `n_paired_*` gives the number of pairs |
| `chimera_junction_pll`, `chimera_junction_min` | mean and minimum PLL over 20 residues around the crossover (10 through the last 5'-derived/hybrid residue + 10 first 3'-derived, clipped to available context) |
| `chimera_junction_dip` | junction-window mean minus the same sequence's domain mean (self-normalised) |
| `parent5_junction_dip`, `parent3_junction_dip` | the same at each parent's own crossover residues |
| `dip_excess_vs_5p`, `dip_excess_vs_3p`, `dip_excess_worst` | chimera dip minus parent dip (and the lower of the two), so a crossover in a naturally low-scoring loop is not counted against the chimera |

Per-residue junction-window profiles are in `junction_profiles.tsv.gz`.

**Rank label, used only to list disagreements.** There is no validated way to combine these
measures or to call a chimeric domain functional from them, so no combined score and no cutoff
is defined. To compare with PR #155's tier and PR #160's fold flag, the PLM needs a categorical
reading, and it is a within-run rank, not a calibrated threshold: `plm_rank = bottom_decile` when
`delta_vs_better_parent` **or** `dip_excess_worst` is in the lowest 10 % of the scored PR #155
records by empirical percentile (fraction of reference values ≤ the record value), otherwise `not_bottom_decile`. PCAWG/TCGA DNA records are ranked
against the PR #155 distribution but do not shape it. Because the reference is the *scored*
subset (priority records + a seeded random sample of the rest), the decile is a sample estimate. Tied values share a percentile; a tie spanning the 10 % boundary is not labelled bottom decile.

**Fail-safe (every rule has a regression test in `tests/test_domain_plm.py`).**

* Missing PLM env or weights, worker crash, worker error, out-of-order or wrong-sequence reply →
  `not_assessed` for every affected record.
* Any NaN, ±inf or positive log-probability, or a reply with missing/extra score positions →
  the whole reply fails, nothing from it is cached, and its records are `not_assessed`. Truncated domain sequences or out-of-bounds domain spans also fail closed.
* Sequence mismatch: each record's domain span, parent domains, junction index and domain length
  are re-checked against PR #155's `domain_function_metrics.tsv.gz`; any difference →
  `not_assessed`. Non-canonical residues → `not_assessed`.
* Cache: `<work>/scores/<spec_id>/<sha256(sequence)>.json`, `spec_id` = hash of model name +
  fair-esm version + weights sha256. Every read re-checks the spec, the sequence and its hash,
  and every value; a mismatched, corrupt or partial entry is treated as absent, so the record is
  `not_assessed` unless re-scored. A record with any required position unscored (coverage limit)
  is `not_assessed`.
* `not_assessed` has no metrics and no rank, and never takes part in a PLM disagreement.
* PR #155's `tier` and PR #160's `structure_flag` are copied unchanged into `pr155_tier` and
  `pr160_fold_flag`; nothing in this run writes to PR #155 or PR #160 outputs (input sha256s are
  in `run_manifest.json`).

## Limits

* **PLMs favour natural-looking and paralog-like sequence.** Most of these chimeras join two
  paralogs (often >80 % identical), so the chimera looks like a family member and scores close to
  its parents almost by construction; a near-zero delta is weak evidence. Conversely, a chimera
  between divergent partners can score low simply for looking unusual.
* **Scores rank; they do not prove function or stability.** No measure here has been calibrated
  against lab data for chimeric domains; there is no functional cutoff, and the bottom-decile
  label is only a rank within this run.
* A 35 M-parameter model is at the small end of ESM-2. The 150 M candidate was rejected after
  the earlier worker's timing pilot; it was not used for the committed metrics. No accuracy gain
  from a larger model was measured in this run.
* Masked PLL scores each residue with all others visible, so it is insensitive to long-range
  incompatibilities that only show up in 3D (which is what PR #160's fold flag is for).
* Windowed context (4 records over 1,022 residues) is not identical to full context.
* Coverage beyond the priority records is a time-limited prefix and reverse-prefix of a seeded
  random ordering, not the full set. Runtime, record length and shared cached sequences can affect
  which records finish before the deadline. Priority cases are deliberately enriched, so the pooled
  rank reference is not an unbiased estimate of the full collection.

## Tools not run

| tool | status | reason |
|---|---|---|
| PROVEAN | not run | [Official requirements](https://www.jcvi.org/research/provean): NCBI BLAST, CD-HIT and NCBI `nr`. A full `nr` installation exceeds the 5 GB budget. A reduced database would change the supporting homologs and would need separate calibration. `blastp`, `psiblast`, `cd-hit` and PROVEAN were absent on PATH; no installation attempted. |
| INDELi | not run | [Official single-protein workflow](https://github.com/lehner-lab/deep_indel_mutagenesis/blob/main/single_protein_prediction/README.md) predicts single-amino-acid insertion/deletion tolerance using ESM-1v substitution scores, a PDB, rSASA and R packages. These chimeras are not single-residue variants of either parent, so there is no faithful input. Not installed; no claim that the complete dependency stack fits within 5 GB. |
| FoldX | not run | not installed (academic licence required); not required by the contract. |
| structure-based stability proxy on the PR #160 ColabFold models | not run | No additional energy model was run; the PLM metrics are sequence plausibility measures, not folding free energies. FoldX is absent and no Rosetta installation was attempted. |

## Reproduce

Use the existing PLM environment (do not reinstall):

```bash
.venv/bin/python scripts/domain_stability_plm.py analyze --out runs/domain_stability_plm_2026-11-13
```

Analysis validates the installed model against the external sequence cache; it does not score
new positions. Without the matching cache, records are `not_assessed`. Model weights and raw
score cache stay outside the committed artifacts. Both CPU and MPS workers use the same model
and fp32 settings, and share cache entries; small device-dependent numerical differences are
possible. `run_manifest.json` records the model hash, input hashes, scoring sessions and order seed.


## Verification

[validation_audit.json](validation_audit.json) records the output checks: all 3,039 upstream tiers preserved, all 27 #160 flags preserved, unchanged upstream run directories versus `09207fc`, complete priority coverage, finite assessed metrics/profiles, blank not-assessed metrics, matching input/log hashes and absent scoring processes. The junction profile file contains 4479 residue rows for 224 records. Regression tests cover missing output/model, non-finite scores, wrong positions/length, stale model/weight/sequence caches, truncated domains, delayed shutdown and preservation of tiers/flags.

Final-commit gates (results reported in the PR):

```bash
.venv/bin/pre-commit run --all-files
.venv/bin/ruff check .
.venv/bin/ruff format --check .
.venv/bin/mypy src scripts/domain_stability_plm.py
.venv/bin/pytest -m "not network"
```
