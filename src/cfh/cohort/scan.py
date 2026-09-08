"""Genome-wide cohort scan orchestration.

Ties together: cohort-wide recurrence gating (:mod:`cfh.cohort.recurrence`),
auto-generated ``GeneConfig`` construction for genes with no curated config
(:mod:`cfh.cohort.auto_config`), the existing per-gene analysis pipeline
(:mod:`cfh.real_benchmark`) run through the existing registered-algorithm
orchestrator (:mod:`cfh.orchestrator.run`), and cross-gene FDR correction
(:mod:`cfh.stats.multiple_testing`, via :mod:`cfh.gene_comparison`).

A single gene's failure (a malformed record, an unmappable breakpoint, a
transient network error for just that gene) is always caught and recorded
as that gene's own outcome -- it never aborts the whole scan.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import requests

from cfh.algorithms.confidence_stats import default_confidence_stats_params
from cfh.algorithms.registry import list_algorithms
from cfh.cohort.auto_config import (
    DEFAULT_GENOME_NEXUS_BASE_URL,
    PfamDescriptionSource,
    batch_fetch_canonical_transcripts,
    build_auto_gene_config,
)
from cfh.cohort.recurrence import (
    DEFAULT_MIN_DISTINCT_PATIENTS,
    GeneRecurrence,
    RecurrenceGateResult,
    fetch_cohort_gene_recurrence,
    gate_genes_by_recurrence,
)
from cfh.gene_comparison import collect_p_values_from_algorithm_results
from cfh.genes.registry import GeneConfig, load_gene_config
from cfh.ingestion import cbioportal_api
from cfh.mapping.genome_nexus_source import GenomeNexusClient
from cfh.real_benchmark import RealBenchmarkRun, analyze_structural_variant_calls_with_config
from cfh.stats.multiple_testing import benjamini_hochberg
from cfh.studies.registry import load_study_config

DEFAULT_N_PERMUTATIONS = 1_000
DEFAULT_N_PERMUTATIONS_SMALL = 100
DEFAULT_SIGNIFICANCE_LEVEL = 0.05
ADAPTIVE_ALGORITHMS = (
    "domain_retention",
    "domain_disruption",
    "cutpoint_detection",
    "window_detection",
)


@dataclass
class GeneScanOutcome:
    """One gene's outcome from a cohort scan: never an exception, always a
    structured record so a whole-run report can be built regardless of how
    many individual genes failed or were skipped.
    """

    gene_symbol: str
    entrez_gene_id: int | None
    distinct_patient_count: int
    total_sv_count: int
    config_source: str  # "curated" | "auto" | "unresolved"
    status: str  # "ok" | "failed"
    run: RealBenchmarkRun | None = None
    error: str | None = None
    p_value_rows: list[dict] = field(default_factory=list)


@dataclass
class CohortScanResult:
    study_id: str
    min_distinct_patients: int
    recurrence_gate: RecurrenceGateResult
    curated_gene_count: int
    auto_config_gene_count: int
    unresolved_gene_count: int
    gene_outcomes: list[GeneScanOutcome]
    fdr_rows: list[dict]
    significant_genes: list[str]
    significance_level: float = DEFAULT_SIGNIFICANCE_LEVEL
    warnings: list[str] = field(default_factory=list)

    @property
    def total_genes_before_gating(self) -> int:
        return self.recurrence_gate.total_genes

    @property
    def genes_after_gating(self) -> int:
        return self.recurrence_gate.passing_count


def _resolve_configs(
    candidate_genes: list[GeneRecurrence],
    *,
    genome_nexus_base_url: str,
    genome_nexus_cache_dir: Path | None,
    pfam_description_cache_dir: Path | None,
    session: "requests.Session | None" = None,
) -> tuple[dict[str, GeneConfig], dict[str, str], list[str]]:
    """Resolve one ``GeneConfig`` per candidate gene: curated configs always
    win; everything else is auto-generated in as few Genome Nexus batch
    calls as possible. Returns ``(config_by_gene, source_by_gene,
    unresolved_gene_symbols)``.
    """
    config_by_gene: dict[str, GeneConfig] = {}
    source_by_gene: dict[str, str] = {}
    needs_auto: list[GeneRecurrence] = []

    for gene in candidate_genes:
        try:
            config_by_gene[gene.hugo_gene_symbol] = load_gene_config(gene.hugo_gene_symbol)
            source_by_gene[gene.hugo_gene_symbol] = "curated"
        except FileNotFoundError:
            needs_auto.append(gene)

    description_source = PfamDescriptionSource(
        session=session, cache_dir=pfam_description_cache_dir
    )
    batch = batch_fetch_canonical_transcripts(
        [gene.hugo_gene_symbol for gene in needs_auto],
        base_url=genome_nexus_base_url,
        cache_dir=genome_nexus_cache_dir,
        session=session,
    )

    unresolved: list[str] = []
    for gene in needs_auto:
        canonical = batch.by_gene_symbol.get(gene.hugo_gene_symbol)
        if canonical is None:
            unresolved.append(gene.hugo_gene_symbol)
            continue
        auto_config = build_auto_gene_config(
            gene.hugo_gene_symbol,
            gene.entrez_gene_id,
            canonical,
            description_source=description_source,
        )
        if auto_config is None:
            unresolved.append(gene.hugo_gene_symbol)
            continue
        config_by_gene[gene.hugo_gene_symbol] = auto_config
        source_by_gene[gene.hugo_gene_symbol] = "auto"

    return config_by_gene, source_by_gene, unresolved


def _adaptive_algorithm_params(adaptive: bool, n_permutations_small: int) -> dict[str, dict]:
    if not adaptive:
        return {}
    return {
        name: {"adaptive": True, "n_permutations_small": n_permutations_small}
        for name in ADAPTIVE_ALGORITHMS
    }


def run_cohort_scan(
    study_id: str,
    *,
    min_distinct_patients: int = DEFAULT_MIN_DISTINCT_PATIENTS,
    n_permutations: int = DEFAULT_N_PERMUTATIONS,
    adaptive: bool = True,
    n_permutations_small: int = DEFAULT_N_PERMUTATIONS_SMALL,
    significance_level: float = DEFAULT_SIGNIFICANCE_LEVEL,
    algorithm_names: list[str] | None = None,
    max_genes: int | None = None,
    cbioportal_base_url: str = cbioportal_api.DEFAULT_BASE_URL,
    cache_dir: str | Path | None = None,
    session: "requests.Session | None" = None,
) -> CohortScanResult:
    """Run the full genome-wide cohort scan for ``study_id``.

    Ingests cohort-wide SV data once (the recurrence call), gates to
    recurrently-altered genes, resolves each gated gene's ``GeneConfig``
    (curated taking precedence over an auto-generated one), runs the full
    registered algorithm suite per gene through the existing orchestrator,
    and applies Benjamini-Hochberg FDR correction across every scanned
    gene's applicable p-values -- not just BRAF/RET.

    ``max_genes`` caps the number of gated candidate genes actually
    analyzed (after sorting by recurrence), for a bounded test/demo run;
    the reported ``total_genes_before_gating``/``genes_after_gating``
    counts are unaffected by this cap and always describe the full cohort.
    """
    algorithm_names = algorithm_names or list_algorithms()
    cache_dir = Path(cache_dir) if cache_dir else None

    study_config = load_study_config(study_id)
    profile_id = (
        study_config.molecular_profile_id(study_id)
        if study_config
        else f"{study_id}_structural_variants"
    )
    genome_nexus_base_url = (
        study_config.genome_nexus_base_url if study_config else DEFAULT_GENOME_NEXUS_BASE_URL
    )

    recurrence = fetch_cohort_gene_recurrence(
        study_id, base_url=cbioportal_base_url, session=session
    )
    gate = gate_genes_by_recurrence(
        recurrence, min_distinct_patients=min_distinct_patients, study_id=study_id
    )

    candidate_genes = gate.passing_genes
    if max_genes is not None:
        candidate_genes = candidate_genes[:max_genes]

    genome_nexus_cache_dir = cache_dir / "genome_nexus_canonical_transcripts" if cache_dir else None
    pfam_description_cache_dir = cache_dir / "pfam_descriptions" if cache_dir else None
    config_by_gene, source_by_gene, unresolved = _resolve_configs(
        candidate_genes,
        genome_nexus_base_url=genome_nexus_base_url,
        genome_nexus_cache_dir=genome_nexus_cache_dir,
        pfam_description_cache_dir=pfam_description_cache_dir,
        session=session,
    )

    genome_nexus_client = GenomeNexusClient(base_url=genome_nexus_base_url, session=session)
    algorithm_params = _adaptive_algorithm_params(adaptive, n_permutations_small)

    outcomes: list[GeneScanOutcome] = []
    warnings: list[str] = []
    for gene in candidate_genes:
        symbol = gene.hugo_gene_symbol
        config = config_by_gene.get(symbol)
        if config is None:
            outcomes.append(
                GeneScanOutcome(
                    gene_symbol=symbol,
                    entrez_gene_id=gene.entrez_gene_id,
                    distinct_patient_count=gene.distinct_patient_count,
                    total_sv_count=gene.total_sv_count,
                    config_source="unresolved",
                    status="failed",
                    error="No canonical transcript/protein could be resolved for this gene.",
                )
            )
            continue

        try:
            if gene.entrez_gene_id is None:
                raise ValueError("cBioPortal recurrence record had no entrezGeneId")
            calls = cbioportal_api.fetch_structural_variants(
                [gene.entrez_gene_id],
                [profile_id],
                base_url=cbioportal_base_url,
                session=session,
            )
            run = analyze_structural_variant_calls_with_config(
                calls,
                config,
                study_id,
                molecular_profile_id=profile_id,
                genome_nexus_client=genome_nexus_client,
                n_permutations=n_permutations,
                algorithm_names=algorithm_names,
                algorithm_params={
                    **algorithm_params,
                    "confidence_stats": default_confidence_stats_params(config),
                },
            )
            p_value_rows = collect_p_values_from_algorithm_results(
                symbol, study_id, run.results, source=f"cohort_scan:{symbol}"
            )
            outcomes.append(
                GeneScanOutcome(
                    gene_symbol=symbol,
                    entrez_gene_id=gene.entrez_gene_id,
                    distinct_patient_count=gene.distinct_patient_count,
                    total_sv_count=gene.total_sv_count,
                    config_source=source_by_gene.get(symbol, "unresolved"),
                    status="ok",
                    run=run,
                    p_value_rows=p_value_rows,
                )
            )
        except Exception as exc:  # noqa: BLE001 - a single gene must never abort the scan
            outcomes.append(
                GeneScanOutcome(
                    gene_symbol=symbol,
                    entrez_gene_id=gene.entrez_gene_id,
                    distinct_patient_count=gene.distinct_patient_count,
                    total_sv_count=gene.total_sv_count,
                    config_source=source_by_gene.get(symbol, "unresolved"),
                    status="failed",
                    error=f"{type(exc).__name__}: {exc}",
                )
            )

    all_p_rows = [row for outcome in outcomes for row in outcome.p_value_rows]
    hypotheses = [
        (row["gene"], f"{row['algorithm']}:{row['test']}", row["raw_p"]) for row in all_p_rows
    ]
    adjusted = benjamini_hochberg(hypotheses)
    fdr_rows = [
        {**row, "bh_adjusted_q": q_value}
        for row, (_, _, _, q_value) in zip(all_p_rows, adjusted, strict=True)
    ]

    per_gene_min_q: dict[str, float] = {}
    for row in fdr_rows:
        gene = row["gene"]
        per_gene_min_q[gene] = min(per_gene_min_q.get(gene, 1.0), row["bh_adjusted_q"])
    significant_genes = sorted(
        gene for gene, q_value in per_gene_min_q.items() if q_value < significance_level
    )

    if unresolved:
        warnings.append(
            f"{len(unresolved)} gated gene(s) had no resolvable canonical transcript in "
            f"Genome Nexus and were skipped: {', '.join(sorted(unresolved))}"
        )

    return CohortScanResult(
        study_id=study_id,
        min_distinct_patients=min_distinct_patients,
        recurrence_gate=gate,
        curated_gene_count=sum(1 for source in source_by_gene.values() if source == "curated"),
        auto_config_gene_count=sum(1 for source in source_by_gene.values() if source == "auto"),
        unresolved_gene_count=len(unresolved),
        gene_outcomes=outcomes,
        fdr_rows=fdr_rows,
        significant_genes=significant_genes,
        significance_level=significance_level,
        warnings=warnings,
    )


def per_gene_min_q_value(result: CohortScanResult) -> dict[str, float]:
    """Convenience accessor: each gene's smallest (most significant)
    BH-adjusted q-value across whichever algorithms/tests it produced."""
    per_gene: dict[str, float] = {}
    for row in result.fdr_rows:
        gene = row["gene"]
        per_gene[gene] = min(per_gene.get(gene, 1.0), row["bh_adjusted_q"])
    return per_gene


def genes_needing_full_report(
    result: CohortScanResult, honorable_mention_genes: "frozenset[str] | set[str]" = frozenset()
) -> list[str]:
    """Every hand-curated gene, every FDR-significant gene, and every
    ``honorable_mention_genes`` (the highly ranked non-FDR-significant "second
    tier" a human reviewer should still see a full report for) -- the genes a
    cohort scan generates a full per-gene report for (everything else stays
    summary-only).

    Curated genes are derived from each outcome's own ``config_source``
    rather than a hardcoded gene list: a hardcoded list silently goes stale
    every time a new curated ``GeneConfig`` is added (this happened once
    already -- ALK and NTRK1 curated configs were added the commit right
    after this list was introduced, and it was never updated to match).
    """
    scanned = {outcome.gene_symbol for outcome in result.gene_outcomes if outcome.run is not None}
    curated = {
        outcome.gene_symbol
        for outcome in result.gene_outcomes
        if outcome.config_source == "curated"
    }
    names = curated | set(result.significant_genes) | set(honorable_mention_genes)
    return sorted(names & scanned)
