# Domain stability and protein-language-model scores

Run date: 2026-11-13. Stacked on PR #160 (`polly/domain-structure`).

Work in progress: both inherited scoring jobs are still running. Final coverage, named-case results, disagreements and validation will be added only after both jobs exit.

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
