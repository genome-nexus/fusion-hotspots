# MSK-IMPACT 50k locus-validation rerun

Live rerun completed on 2026-09-04 UTC with the same `>=5` distinct-patient gate,
adaptive 100/1,000-permutation settings, and real cBioPortal and Genome Nexus
services used by the committed 2026-09-03 scan.

## Genome-wide result

| Metric | Previously committed scan | Corrected scan |
|---|---:|---:|
| Genes with SV records | 3,919 | 3,919 |
| Genes passing the gate | 544 | 544 |
| Successfully analyzed genes | 524 | 523 |
| Genes with a numeric retention change | 0 (baseline) | 225 |
| Additional status change | none | NCOA4: `ok` to `failed` |
| Genes with a changed frame percentage while still analyzable | 0 | 0 |
| Genes with a changed Fisher p-value | 0 (baseline) | 113 |
| FDR-significant genes (`q<0.05`) | ETV6, RET | ETV6 |

In total, **226 of 544 gated gene outcomes changed materially**: 225 genes had a
different numeric domain-retention percentage (93 increased and 132 decreased),
and NCOA4 no longer produces a numeric result because Genome Nexus returned no
exon coordinates with which to validate its target locus. The complete 226-row
comparison is in `runs/cohort_scan_locus_validation_comparison_20260904.tsv`.

The large count is expected from the bug's scope: auto-configured genes were also
selecting target coordinates from site labels that are not reliable enough to
identify the target locus. Some 100-point changes are based on only one or two
events, so their percentages should not be interpreted as equally important.
Among genes with at least 20 analyzed events, illustrative large changes include
FGFR3 (39.47% to 93.42%, n=152), CCDC6 (43.40% to 1.89%, n=53), TMPRSS2
(43.37% to 6.11%, n=867), EML4 (29.15% to 53.36%, n=223), TACC3 (26.52% to
48.48%, n=132), and RET (74.23% to 92.27%, n=194).

## Recomputed headline genes

| Gene | Events | Retention, old to corrected | Fisher p, old to corrected | Genome-wide q, old to corrected | Corrected significant? |
|---|---:|---:|---:|---:|---:|
| ETV6 | 90 | 75.56% to 75.56% | 5.123e-06 to 5.123e-06 | 0.004437 to 0.004334 | yes |
| RET | 194 | 74.23% to 92.27% | 9.841e-05 to 0.0004197 | 0.04261 to 0.1198 | no |
| BRAF | 179 | 91.06% to 91.06% | 0.01337 to 0.01337 | 0.1453 to 0.2356 | no |

The old headline does **not** hold. BRAF still does not survive genome-wide BH
correction, but corrected RET also does not survive it; ETV6 is the only
FDR-significant gene in the corrected scan.

## Attribution check

To distinguish the fix from live-source drift, the immediate pre-fix commit
`5b50a55` was run against the same live 50k source and the same cached Genome
Nexus annotations. Across all 544 rows it reproduced the committed scan exactly
for status, distinct-patient count, SV count, analyzed-event count, frame
percentage, retention percentage, and Fisher p-value. The deterministic old/new
differences above are therefore attributable to locus validation, not 50k source
drift. Monte Carlo permutation p-values can vary slightly between independent
runs and are not counted as material changes by themselves.

## Individual benchmarks

| Benchmark | Old result | Corrected result | Conclusion |
|---|---|---|---|
| BRAF, 50k | 163/179 retained (91.06%), Fisher p=0.01337 | 163/179 retained (91.06%), Fisher p=0.01337 | Deterministic retention/frame result unchanged |
| RET, 50k | 144/194 retained (74.23%), Fisher p=9.841e-05 | 179/194 retained (92.27%), Fisher p=0.0004197 | Material correction |
| BRAF, 2017 | 31/35 retained (88.57%), Fisher p=0.01008 | 33/41 retained (80.49%), Fisher p=0.0006072 | Not a fix-only comparison; see below |

For BRAF 2017, running the exact old artifact commit (`a8767ee`) against today's
live services reproduces the old 35 events and all old deterministic results,
including the same event identities. Thus the earlier source-drift ambiguity is
resolved: the old artifact is reproducible and current evidence does not show
source drift. The increase from 35 to 41 fusion-annotated records comes from
subsequent ingestion changes already present before locus validation.

On identical current code/data except for locus validation, the immediate
pre-fix run maps 41/41 records and the corrected run maps 40/41. The fix rejects
one malformed PARP12-BRAF record whose two genomic coordinates are both `-1`.
In-frame and retained counts both remain 33, while the contingency table changes
from `[[31, 2], [2, 6]]` to `[[31, 2], [2, 5]]` and Fisher p changes from
0.0001575 to 0.0006072. Therefore the locus fix affects BRAF 2017's mapped count
and association statistic, but not its frame or retained count on today's
ingestion path.

## Artifacts

The specific timestamped run directories originally cited here have since
been superseded by later fix rounds and pruned per the keep-latest-run
convention (see `CONTRIBUTING.md`); the numeric claims above still match the
current committed data. The closest surviving runs for the same
gene-or-scan-kind/study_id groups are:

- Corrected full scan: `runs/cohort-scan_msk_impact_50k_2026_20260904T144201Z/cohort_scan/`
- Corrected BRAF 50k: `runs/braf_msk-impact-50k-2026_20260904T172738Z/`
- Corrected RET 50k: `runs/ret_msk-impact-50k-2026_20260904T172752Z/`
- Refreshed BRAF 2017: `runs/braf_msk-impact-2017_20260904T005539Z/`
