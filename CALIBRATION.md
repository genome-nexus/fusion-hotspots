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

## Reproducible pilot

The [saved pilot](runs/calibration_20260929_seed42.json) used 50 replicates per
scenario, 24 patients, 99 permutations, seed 42, and a 100-aa window. Rates below
are scan-level rejections at `p < 0.05`, not across-gene FDR estimates.

| Simulation | Cutpoint | Window |
| --- | ---: | ---: |
| Independent null | 3/50 (6%) | 1/50 (2%) |
| Uneven-position null | 4/50 (8%) | 1/50 (2%) |
| Planted cutpoint | 48/50 (96%) | 11/50 (22%) |
| Planted narrow window | 3/50 (6%) | 7/50 (14%) |
| Mapped boundary pile-up null | 2/50 (4%) | 2/50 (4%) |
| Repeated-patient null | 30/50 (60%) | 44/50 (88%) |

The JSON includes Wilson 95% intervals. With only 50 replicates, these are
preliminary estimates for the specified simulations. For example, the repeated-
patient null intervals are 46.2–72.4% and 76.2–94.4%; the independent-null intervals
are 2.1–16.2% and 0.4–10.5%.

The repeated-patient scenario copies each patient's position and label three
times, then tests those copies with the existing event-level label permutation.
The result demonstrates a violated exchangeability assumption. **Patient-aware
inference is still required**; retaining patient IDs, sample-level expression
deduplication, and overlap guards do not repair that permutation null. The weak
planted-window sensitivity also argues for assessing power under realistic event
counts and effect sizes before interpreting an absent signal. No external ranking
validation was performed because independent labeled inputs were not supplied.

## Patient-level observation unit

The live pipeline now collapses repeated observations of one fusion in one
patient before any inferential test (`observation_unit="patient"`, the default;
see `cfh.stats.observation_units`). The same collapse is available to the
simulations:

```sh
python -m cfh.stats.calibration --replicates 200 --n-permutations 99 --seed 42 --observation-unit event
python -m cfh.stats.calibration --replicates 200 --n-permutations 99 --seed 42 --observation-unit patient
```

Both runs are checked in as
[`calibration_20260929_seed42_n200_event.json`](runs/calibration_20260929_seed42_n200_event.json)
and [`calibration_20260929_seed42_n200_patient.json`](runs/calibration_20260929_seed42_n200_patient.json).
They use identical seeds, so every scenario except the repeated-patient null is
unchanged. Cells are rejections/200 at `p < 0.05` (Wilson 95% interval):

| Simulation | Cutpoint, event unit | Cutpoint, patient unit | Window, event unit | Window, patient unit |
| --- | ---: | ---: | ---: | ---: |
| Independent null | 5 (1.1–5.7%) | 5 (1.1–5.7%) | 4 (0.8–5.0%) | 4 (0.8–5.0%) |
| Uneven-position null | 6 (1.4–6.4%) | 6 (1.4–6.4%) | 9 (2.4–8.3%) | 9 (2.4–8.3%) |
| Mapped boundary pile-up null | 8 (2.0–7.7%) | 8 (2.0–7.7%) | 7 (1.7–7.0%) | 7 (1.7–7.0%) |
| **Repeated-patient null** | **130 (58.2–71.3%)** | **11 (3.1–9.6%)** | **168 (78.3–88.4%)** | **4 (0.8–5.0%)** |
| Planted cutpoint | 192 (92.3–98.0%) | 192 (92.3–98.0%) | 53 (20.9–33.0%) | 53 (20.9–33.0%) |
| Planted narrow window | 5 (1.1–5.7%) | 5 (1.1–5.7%) | 36 (13.3–23.9%) | 36 (13.3–23.9%) |

Collapsing identical repeat observations restores the repeated-patient null to
roughly nominal rejection. It does not address dependence between *distinct*
fusions in one patient, which the simulations do not model.

## Cohort-scan FDR family without the read-support t-test

The `confidence_stats` Welch test (tumor read support, retained vs
non-retained events) is no longer part of the cross-gene BH family (#114). On
the MSK-IMPACT 50k cohort (`cfh cohort-scan msk_impact_50k_2026`, 544 genes
past the 5-patient gate, 538 analyzable), this changes the result:

| Gene | Min q, with Welch | Min q, without | Retention Fisher q, without | Welch groups (retained / not) |
| --- | ---: | ---: | ---: | ---: |
| ETV6 | 0.049 | 0.158 | 0.158 | 46 / 19 |
| NTRK3 | 0.0015 | 0.166 | 0.883 | 35 / 2 |
| ROS1 | 0.0017 | 0.166 | 1.000 | 89 / 11 |
| FLI1 | 0.023 | 0.166 | 0.701 | 99 / 7 |

With the Welch test, 4 genes were FDR-significant, 3 of them only because of
it. Without it, **no gene reaches q < 0.05**; the smallest q is ETV6's 0.158.
ETV6's q rises even though its p-value is unchanged, because the many tiny
Welch p-values no longer rank ahead of it.

These numbers come from two scans of this branch on 2026-09-29. Transient
cBioPortal 503 errors and timeouts dropped 8 genes from the first run and 17
different genes from the second, but every analyzable gene completed in at
least one run. Both runs found no significant genes (min q 0.163 and 0.158),
and every completed gene's Fisher and permutation p-values match the earlier
scan exactly. The table shows the second run.

## Power of the cutpoint and window scans

```sh
python -m cfh.stats.calibration --power-grid --replicates 100 --n-permutations 99 --seed 42
```

The run is checked in as
[`power_20260929_seed42.json`](runs/power_20260929_seed42.json). Each cell plants
one signal in independent events at uniform positions on a 1,000-aa protein:
either a cutpoint at 500 aa or a window at 400–500 aa. The effect is the absolute
difference in retained probability inside vs outside (0.8 means 90% vs 10%). The
window scan uses the production widths (25, 50, 100, 200 aa). Cells are
rejections/100 at `p < 0.05` (Wilson 95% interval):

| Events | Cutpoint, effect 0.4 | Cutpoint, 0.6 | Cutpoint, 0.8 | Window, effect 0.4 | Window, 0.6 | Window, 0.8 |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 10 | 5 (2–11%) | 15 (9–23%) | 41 (32–51%) | 0 (0–4%)\* | 1 (0–5%)\* | 0 (0–4%)\* |
| 25 | 30 (22–40%) | 75 (66–82%) | 96 (90–98%) | 7 (3–14%) | 12 (7–20%) | 30 (22–40%) |
| 50 | 57 (47–66%) | 98 (93–99%) | 100 (96–100%) | 12 (7–20%) | 30 (22–40%) | 63 (53–72%) |
| 100 | 92 (85–96%) | 100 (96–100%) | 100 (96–100%) | 20 (13–29%) | 70 (60–78%) | 91 (84–95%) |
| 250 | 100 (96–100%) | 100 (96–100%) | 100 (96–100%) | 73 (64–81%) | 100 (96–100%) | 100 (96–100%) |

\* At 10 events the window scan could not produce a result in 10, 12 and 24 of
the 100 replicates, which count as non-rejections.

The cutpoint scan has at least 75% power from 25 events at a strong effect
(0.6) and from 100 events at a moderate one (0.4). The window scan is much
weaker: below 50 events it rarely detects even a strong window, and at 100
events it needs an effect of about 0.6. In the MSK-IMPACT 50k scan the median
gene has 5 analyzable events; 490 of 538 genes have fewer than 25 and only 13
have 100 or more. For most genes, then, even the cutpoint scan is underpowered,
and an absent signal is weak evidence of absence. These are single-scan rates before cross-gene FDR,
with independent events (repeated-patient dependence is not modeled), uniform
positions (real breakpoints cluster), and no adaptive escalation beyond the
first permutation stage.
