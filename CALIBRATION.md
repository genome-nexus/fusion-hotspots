# Offline calibration and held-out ranking checks

Run `python -m cfh.stats.calibration --replicates 50 --n-permutations 99 --seed 42` to simulate null labels, uneven breakpoint density, intronic boundary snapping, repeated patient observations, and planted cutpoint/window signals. The JSON reports rejection rates and Wilson 95% intervals at the strict `p < 0.05` decision rule. The cutpoint and window tests permute outcome labels and use the maximum scanned statistic. These simulations do not test genomic breakpoint density, cross-gene FDR, or real clinical validity. The repeated-patient null intentionally violates event-level exchangeability; a high rejection rate there is a warning, not a corrected estimate.

An empirical scan p-value from `N` permutations cannot be smaller than `1/(N+1)`. For a predeclared family of `m` hypotheses, pass `correction_family_size=m` with a fixed `n_permutations` to the algorithm. The budget resolver rejects adaptive mode and budgets too small for the *raw* first-rank BH or Bonferroni threshold to fall strictly below `0.05`. This is a precision check. It does not establish family-wide FDR control, valid exchangeability, or independence between patients. Candidate positions within one scan are already handled by that scan's maximum-statistic permutation and should not be counted as the across-gene family.

For an external ranking comparison, prepare two JSON files and run:

```sh
python -m cfh.stats.independent_validation discovery.json validation.json --k 5
```

`discovery.json` declares patients used to fit both frozen scores:

```json
{"patient_id_namespace": "registry-v1", "patient_ids": ["D1", "D2"]}
```

`validation.json` declares patients used for external labels and the same candidate set for both scores:

```json
{
  "patient_id_namespace": "registry-v1",
  "patient_ids": ["V1", "V2"],
  "candidates": [
    {"candidate_id": "PARTNER_A", "label": 1, "composite_score": 0.8, "recurrence_score": 0.4},
    {"candidate_id": "PARTNER_B", "label": 0, "composite_score": 0.2, "recurrence_score": 0.6}
  ]
}
```

Candidate labels must be independently sourced binary truth labels, and scores must be frozen from discovery patients before examining validation labels. The tool checks complete, comparable, disjoint patient IDs and finite paired scores. It reports descriptive average precision and precision at `k`, plus differences. Equal scores are grouped for average precision, and boundary ties receive fractional credit for precision at `k`; candidate order does not decide ties. The tool cannot verify the origins of the scores or labels. Without externally labeled, patient-disjoint data, no independent validation result can be claimed.

`cfh compare-cohorts --require-patient-disjoint` adds a stricter guard to the pooled cohort comparison: every supplied event must carry a patient ID, both artifacts must declare the same `patient_id_namespace`, and neither patient nor sample IDs may overlap. Existing saved artifacts without that metadata can still produce nominal descriptive comparisons, but cannot pass the strict guard. Matching strings alone do not prove that external registries resolved patient identity correctly.
