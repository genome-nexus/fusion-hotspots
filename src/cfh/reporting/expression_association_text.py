"""Templated report sentence for the ``expression_association`` algorithm.

Pure, deterministic string formatting from already-computed
``AlgorithmResult.Summary`` fields -- never free-text/LLM-generated, and
never states a claim without a real computed number backing it (mirroring
``cfh.reporting.domain_names.domain_interpretation_sentence``). Returns
``None`` when nothing was computed, so the caller can state absence
explicitly (e.g. by falling back to the algorithm's own warning) rather than
silently omit the section.
"""

from __future__ import annotations

from typing import Any, Mapping, Optional


def _comparison_sentence(
    lead_in: str,
    label_a: str,
    label_b: str,
    block: Mapping[str, Any],
    n_a_key: str,
    n_b_key: str,
    *,
    alpha: float = 0.05,
) -> Optional[str]:
    p_value = block.get("p_value")
    mean_a = block.get("mean_a")
    mean_b = block.get("mean_b")
    test = block.get("test")
    n_a = block.get(n_a_key)
    n_b = block.get(n_b_key)
    if p_value is None or mean_a is None or mean_b is None or test is None:
        return None
    direction = "higher" if mean_a > mean_b else ("lower" if mean_a < mean_b else "no different")
    significance = "significantly " if p_value < alpha else "not significantly "
    if direction == "no different":
        comparison = f"is {significance}different"
    else:
        comparison = f"is {significance}{direction}"
    return (
        f"{lead_in} {comparison} in {label_a} (n={n_a}) than {label_b} (n={n_b}) "
        f"({test}, p={p_value:.3g})."
    )


def expression_association_sentence(
    gene_symbol: str, expression_summary: Optional[Mapping[str, Any]]
) -> Optional[str]:
    """Render the expression-association finding(s) as templated sentences.

    ``expression_summary`` is an ``expression_association`` result's
    ``Summary`` dict. Returns ``None`` when neither comparison produced a
    p-value (e.g. a no-op result, or every comparison was skipped for
    insufficient sample size) -- the caller is expected to fall back to the
    algorithm's own warning text in that case, not silently drop the section.
    """
    if not expression_summary:
        return None
    sentences: list[str] = []

    fp_vs_fn = expression_summary.get("fusion_positive_vs_negative")
    if fp_vs_fn:
        sentence = _comparison_sentence(
            f"{gene_symbol} mRNA expression",
            "fusion-positive samples",
            "fusion-negative samples",
            fp_vs_fn,
            "n_fusion_positive",
            "n_fusion_negative",
        )
        if sentence:
            sentences.append(sentence)

    split = expression_summary.get("domain_retention_split")
    if split:
        sentence = _comparison_sentence(
            f"Among fusion-positive samples, {gene_symbol} expression",
            "kinase-domain-retained fusions",
            "not-retained fusions",
            split,
            "n_a",
            "n_b",
        )
        if sentence:
            sentences.append(sentence)

    return " ".join(sentences) if sentences else None
