"""Expression-association evidence algorithm.

Does the target gene's own mRNA expression differ (a) between fusion-positive
and fusion-negative samples, and/or (b) between kinase-domain-retained and
not-retained fusion-positive samples?

Gene-agnostic and data-availability-agnostic: this algorithm never assumes
expression data exists. It is opt-in via ``params["expression_by_sample"]``
(a ``{Sample_id: expression_value}`` mapping already fetched by the caller,
e.g. ``real_benchmark.py``, from cBioPortal's molecular-data endpoint for the
target gene's Entrez id in a cohort's mRNA-expression molecular profile). A
caller running against a cohort with no such profile (e.g. a targeted DNA
panel like ``msk_impact_50k_2026``, which carries no ``MRNA_EXPRESSION``
molecular profile at all) simply omits ``expression_by_sample``, and this
algorithm produces a clean no-op result -- the same graceful-skip pattern
already used by ``confidence_stats``/``exon_retention``/``domain_disruption``
for their own optional inputs.

Two independent comparisons, each requiring >=2 non-null expression
observations per group:

``fusion_positive_vs_negative``
    The target gene's expression in samples carrying a mapped protein
    fusion for this gene (from ``events``) vs. every other profiled sample
    in the cohort. Requires ``params["cohort_sample_ids"]`` (every
    ``Sample_id`` covered by the expression molecular profile) so the
    "negative" set can be determined -- without it, only comparison (b)
    can run.

``domain_retention_split``
    Among fusion-positive samples only, expression in kinase-domain-retained
    vs. not-retained fusions, reusing exactly the same group-field/key/value
    map that
    :func:`cfh.algorithms.confidence_stats.default_confidence_stats_params`
    derives from ``gene_config.key_domains`` -- not a re-derived
    domain-collapsing rule of this module's own.

Statistical test choice
------------------------
Welch's t-test (unequal-variance, via the same
:func:`cfh.stats.ttest.welch_t_test` helper ``confidence_stats.py`` itself
calls -- not re-derived here) is the default, matching this project's
existing Welch usage elsewhere. But the small groups this algorithm
typically compares (a handful of fusion-positive cases) frequently violate
the approximate-normality assumption Welch's test still relies on, and
Shapiro-Wilk normality testing is not meaningfully informative below n=3.
So: when both groups have n>=3, a Shapiro-Wilk check (``scipy.stats.shapiro``,
alpha=0.05) on each group decides the test -- Welch's t-test when neither
group's normality is rejected, Mann-Whitney U (two-sided,
``scipy.stats.mannwhitneyu``) otherwise. Below n=3 in either group,
Mann-Whitney U is used unconditionally, since Shapiro-Wilk cannot assess
normality that small and Welch's t-test's Type-I-error control degrades
sharply at those sizes; Mann-Whitney U makes no distributional assumption
and is the standard fallback for small-sample expression comparisons in
genomics.
"""

from __future__ import annotations

from datetime import datetime, timezone
from typing import Any, Optional

from scipy.stats import mannwhitneyu, shapiro

from cfh.algorithms.base import Algorithm
from cfh.algorithms.confidence_stats import default_confidence_stats_params
from cfh.algorithms.registry import register
from cfh.genes.registry import GeneConfig
from cfh.model.algorithm_result import AlgorithmResult
from cfh.model.fusion_event import FusionEvent
from cfh.model.fusion_feature import FusionFeature
from cfh.stats.ttest import welch_t_test

ALGORITHM_NAME = "expression_association"
ALGORITHM_VERSION = "0.2.0"

_SHAPIRO_MIN_N = 3
_NORMALITY_ALPHA = 0.05
_MIN_GROUP_N = 2


def _run_two_group_test(values_a: list[float], values_b: list[float]) -> dict[str, Any]:
    """Compare two numeric samples, picking Welch's t-test or Mann-Whitney U.

    See the module docstring for the full normality-based decision rule.
    """
    if len(values_a) < _MIN_GROUP_N or len(values_b) < _MIN_GROUP_N:
        raise ValueError(f"each group must have at least {_MIN_GROUP_N} observations")

    normal_enough = False
    if len(values_a) >= _SHAPIRO_MIN_N and len(values_b) >= _SHAPIRO_MIN_N:
        try:
            normal_enough = (
                shapiro(values_a).pvalue >= _NORMALITY_ALPHA
                and shapiro(values_b).pvalue >= _NORMALITY_ALPHA
            )
        except ValueError:
            # A degenerate (e.g. constant-valued) group makes Shapiro-Wilk
            # raise rather than return a meaningful p-value; treat that as
            # "not established as normal" rather than let the check fail.
            normal_enough = False

    if normal_enough:
        welch = welch_t_test(values_a, values_b)
        return {
            "test": "welch_t_test",
            "statistic": welch["t_statistic"],
            "p_value": welch["p_value"],
            "mean_a": welch["mean_a"],
            "mean_b": welch["mean_b"],
        }

    result = mannwhitneyu(values_a, values_b, alternative="two-sided")
    return {
        "test": "mann_whitney_u",
        "statistic": float(result.statistic),
        "p_value": float(result.pvalue),
        "mean_a": sum(values_a) / len(values_a),
        "mean_b": sum(values_b) / len(values_b),
    }


def _merged_row(event: FusionEvent, feature: Optional[FusionFeature]) -> dict[str, Any]:
    """Join one event with its matching feature, feature fields taking
    precedence on a name collision -- the same join semantics
    ``confidence_stats._build_rows`` uses for a matched (event, feature)
    pair, applied here to one already-matched pair at a time."""
    if feature is None:
        return event.model_dump()
    return {**event.model_dump(), **feature.model_dump()}


@register(ALGORITHM_NAME)
class ExpressionAssociationAlgorithm(Algorithm):
    """Compare the target gene's own mRNA expression across fusion status
    and, when configured, kinase-domain-retention status.

    See the module docstring for the full comparison/test-selection design.

    Expected ``params`` keys:
        expression_by_sample (dict, required to run any comparison): maps
            ``Sample_id`` to the target gene's expression value in that
            sample. Omit (or pass empty) for a clean no-op result -- the
            expected outcome for a cohort with no mRNA-expression molecular
            profile at all.
        cohort_sample_ids (list, optional): every ``Sample_id`` assessed for
            both fusion status and expression in this cohort. Required only for the
            ``fusion_positive_vs_negative`` comparison, to determine the
            fusion-negative sample set; without it, that comparison is
            skipped (with a warning) and only ``domain_retention_split`` can
            run.
        expression_field (str, optional): human-readable label for the
            expression measurement, echoed into ``Summary`` and used by the
            report sentence (e.g. ``"mRNA expression z-score"``). Purely
            descriptive; does not affect any computation.
    """

    def run(
        self,
        events: list[FusionEvent],
        features: list[FusionFeature],
        gene_config: Optional[GeneConfig],
        params: dict,
    ) -> AlgorithmResult:
        params = params or {}
        gene_name = str(
            ((gene_config.gene_symbol or gene_config.gene_pair) if gene_config else None)
            or "This gene"
        )
        expression_by_sample: dict[str, float] = params.get("expression_by_sample") or {}
        if not expression_by_sample:
            return self._no_op_result(
                f"No mRNA expression data was available for {gene_name} in this cohort; "
                "expression-association analysis was skipped."
            )

        expression_field = params.get("expression_field", "mRNA expression")
        warnings: list[str] = []
        summary: dict[str, Any] = {"expression_field": expression_field, "analysis_unit": "sample"}

        summary_fp_fn, warning_fp_fn = self._fusion_positive_vs_negative(
            events, expression_by_sample, params.get("cohort_sample_ids"), gene_name
        )
        if summary_fp_fn is not None:
            summary["fusion_positive_vs_negative"] = summary_fp_fn
        if warning_fp_fn:
            warnings.append(warning_fp_fn)

        summary_split, warning_split = self._domain_retention_split(
            events, features, gene_config, expression_by_sample, gene_name
        )
        if summary_split is not None:
            summary["domain_retention_split"] = summary_split
        if warning_split:
            warnings.append(warning_split)

        # Expression values are keyed by sample. Repeated biopsies from a
        # known patient therefore remain separate observations; flag the
        # resulting dependence instead of presenting them as independent
        # patients.
        samples_by_patient: dict[tuple[str, str], set[str]] = {}
        for event in events:
            if event.Patient_id and event.Sample_id in expression_by_sample:
                samples_by_patient.setdefault((event.Cohort, event.Patient_id), set()).add(
                    event.Sample_id
                )
        repeated_patient_count = sum(len(samples) > 1 for samples in samples_by_patient.values())
        if repeated_patient_count:
            warnings.append(
                f"{gene_name}: expression comparisons use one observation per sample; "
                f"{repeated_patient_count} known "
                f"{'patient has' if repeated_patient_count == 1 else 'patients have'} "
                "multiple profiled samples, "
                "so observations may not be independent."
            )

        if "fusion_positive_vs_negative" not in summary and "domain_retention_split" not in summary:
            summary = {}

        return AlgorithmResult(
            Algorithm=ALGORITHM_NAME,
            Algorithm_version=ALGORITHM_VERSION,
            Parameters={"expression_field": expression_field},
            Summary=summary,
            Warnings=warnings,
            Created_at=datetime.now(timezone.utc),
        )

    @staticmethod
    def _fusion_positive_vs_negative(
        events: list[FusionEvent],
        expression_by_sample: dict[str, float],
        cohort_sample_ids: Optional[list[str]],
        gene_name: str,
    ) -> tuple[Optional[dict[str, Any]], Optional[str]]:
        if not cohort_sample_ids:
            return None, (
                f"{gene_name}: fusion-positive-vs-negative expression comparison skipped; "
                "no cohort_sample_ids was supplied to determine the fusion-negative group."
            )
        positive_sample_ids = {event.Sample_id for event in events if event.Sample_id} & set(
            cohort_sample_ids
        )
        negative_sample_ids = set(cohort_sample_ids) - positive_sample_ids
        positive_values = [
            expression_by_sample[sample_id]
            for sample_id in positive_sample_ids
            if sample_id in expression_by_sample
        ]
        negative_values = [
            expression_by_sample[sample_id]
            for sample_id in negative_sample_ids
            if sample_id in expression_by_sample
        ]
        if len(positive_values) < _MIN_GROUP_N or len(negative_values) < _MIN_GROUP_N:
            return None, (
                f"{gene_name}: fusion-positive-vs-negative expression comparison skipped; "
                f"each group needs >={_MIN_GROUP_N} samples with expression data "
                f"(fusion-positive={len(positive_values)}, fusion-negative={len(negative_values)})."
            )
        test_result = _run_two_group_test(positive_values, negative_values)
        return {
            "analysis_unit": "sample",
            "n_fusion_positive": len(positive_values),
            "n_fusion_negative": len(negative_values),
            **test_result,
        }, None

    @staticmethod
    def _domain_retention_split(
        events: list[FusionEvent],
        features: list[FusionFeature],
        gene_config: Optional[GeneConfig],
        expression_by_sample: dict[str, float],
        gene_name: str,
    ) -> tuple[Optional[dict[str, Any]], Optional[str]]:
        domain_params = default_confidence_stats_params(gene_config) if gene_config else {}
        group_field = domain_params.get("group_field")
        if not group_field:
            return None, (
                f"{gene_name}: domain-retention expression-split comparison skipped; "
                "no key_domains configured for this gene."
            )
        group_key = domain_params.get("group_key")
        group_value_map = domain_params.get("group_value_map") or {}

        features_by_event = {feature.Event_id: feature for feature in features}
        statuses_by_sample: dict[str, set[str]] = {}
        for event in events:
            if not event.Sample_id or event.Sample_id not in expression_by_sample:
                continue
            row = _merged_row(event, features_by_event.get(event.Event_id))
            raw_value = row.get(group_field)
            if group_key is not None:
                raw_value = raw_value.get(group_key) if isinstance(raw_value, dict) else None
            mapped_value = group_value_map.get(raw_value, raw_value)
            if mapped_value not in ("retained", "not_retained"):
                continue
            statuses_by_sample.setdefault(event.Sample_id, set()).add(mapped_value)

        ambiguous_samples = {
            sample_id for sample_id, statuses in statuses_by_sample.items() if len(statuses) > 1
        }
        retained_values = [
            expression_by_sample[sample_id]
            for sample_id, statuses in statuses_by_sample.items()
            if statuses == {"retained"}
        ]
        not_retained_values = [
            expression_by_sample[sample_id]
            for sample_id, statuses in statuses_by_sample.items()
            if statuses == {"not_retained"}
        ]
        ambiguous_note = (
            f" {len(ambiguous_samples)} "
            f"{'sample' if len(ambiguous_samples) == 1 else 'samples'} "
            "with both retained and not-retained "
            f"events {'was' if len(ambiguous_samples) == 1 else 'were'} excluded."
            if ambiguous_samples
            else ""
        )

        if len(retained_values) < _MIN_GROUP_N or len(not_retained_values) < _MIN_GROUP_N:
            return None, (
                f"{gene_name}: domain-retention expression-split comparison skipped; "
                f"each group needs >={_MIN_GROUP_N} fusion-positive samples with expression "
                f"data (retained={len(retained_values)}, not_retained={len(not_retained_values)})."
                f"{ambiguous_note}"
            )
        test_result = _run_two_group_test(retained_values, not_retained_values)
        return {
            "analysis_unit": "sample",
            "group_a_label": "retained",
            "group_b_label": "not_retained",
            "n_a": len(retained_values),
            "n_b": len(not_retained_values),
            "n_ambiguous_samples_excluded": len(ambiguous_samples),
            **test_result,
        }, (
            f"{gene_name}: domain-retention expression split:{ambiguous_note}"
            if ambiguous_samples
            else None
        )

    @staticmethod
    def _no_op_result(warning: str) -> AlgorithmResult:
        return AlgorithmResult(
            Algorithm=ALGORITHM_NAME,
            Algorithm_version=ALGORITHM_VERSION,
            Parameters={},
            Summary={},
            Warnings=[warning],
            Created_at=datetime.now(timezone.utc),
        )
