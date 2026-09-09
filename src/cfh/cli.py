"""Command-line entry point for cfh."""

from __future__ import annotations

from pathlib import Path

import click

from cfh.cohort.outputs import DEFAULT_HONORABLE_MENTION_COUNT, write_cohort_scan_outputs
from cfh.cohort.recurrence import DEFAULT_MIN_DISTINCT_PATIENTS
from cfh.cohort.scan import DEFAULT_N_PERMUTATIONS_SMALL, run_cohort_scan
from cfh.gene_comparison import compare_gene_runs, write_comparison_tsv
from cfh.genes.registry import available_genes, load_gene_config
from cfh.real_benchmark import RealBenchmarkError, run_analysis, run_real_benchmark, write_outputs


@click.group()
def main() -> None:
    """Cancer Fusion Hotspots command-line tools."""


@main.command("list-genes")
def list_genes() -> None:
    """List gene symbols with a registered config."""
    for symbol in available_genes():
        click.echo(symbol)


@main.command("show-gene")
@click.argument("gene_symbol")
def show_gene(gene_symbol: str) -> None:
    """Print the resolved GeneConfig for a gene symbol."""
    config = load_gene_config(gene_symbol)
    click.echo(config.model_dump_json(indent=2))


@main.command("compare-genes")
@click.argument(
    "run_artifacts",
    nargs=-1,
    required=True,
    type=click.Path(path_type=Path, exists=True),
)
@click.option(
    "--output",
    "output_path",
    required=True,
    type=click.Path(path_type=Path, dir_okay=False),
    help="Path for the adjusted-p-value TSV report.",
)
def compare_genes(run_artifacts: tuple[Path, ...], output_path: Path) -> None:
    """BH-adjust p-values from existing run directories or results.json files."""
    try:
        rows = compare_gene_runs(list(run_artifacts))
        if not rows:
            raise click.ClickException("No applicable p-values were found in the run artifacts.")
        write_comparison_tsv(rows, output_path)
    except ValueError as exc:
        raise click.ClickException(str(exc)) from None
    click.echo(f"Adjusted {len(rows)} p-values across {len(run_artifacts)} run artifacts")
    click.echo(f"report: {output_path}")


@main.command("real-benchmark")
@click.argument("gene_symbol")
@click.argument("study_id")
@click.option(
    "--output-dir",
    type=click.Path(path_type=Path, file_okay=False),
    default=Path("runs"),
    show_default=True,
)
@click.option("--n-permutations", type=click.IntRange(min=1), default=1_000, show_default=True)
@click.option("--output-stem", help="Override the output filename stem.")
@click.option(
    "--pdf/--no-pdf",
    default=True,
    show_default=True,
    help="Also render a self-contained report.pdf for human reviewers.",
)
def real_benchmark(
    gene_symbol: str,
    study_id: str,
    output_dir: Path,
    n_permutations: int,
    output_stem: str | None,
    pdf: bool,
) -> None:
    """Run the prototype live cBioPortal/Genome Nexus benchmark."""
    try:
        run = run_real_benchmark(gene_symbol, study_id, n_permutations=n_permutations)
        paths = write_outputs(
            run,
            output_dir,
            output_stem=output_stem,
            pdf=pdf,
            cli_args=[
                "real-benchmark",
                gene_symbol,
                study_id,
                "--output-dir",
                str(output_dir),
                "--n-permutations",
                str(n_permutations),
            ],
        )
    except RealBenchmarkError as exc:
        raise click.ClickException(str(exc)) from None
    except Exception as exc:
        raise click.ClickException(f"Benchmark failed: {type(exc).__name__}: {exc}") from None
    fisher_p_value = run.summary["fisher_p_value"]
    fisher_display = "unavailable" if fisher_p_value is None else f"{fisher_p_value:.6g}"
    if run.is_gene_pair:
        gene5, gene3 = run.summary["gene_pair"]
        click.echo(
            f"Analyzed {run.summary['eligible_event_count']} {run.gene_symbol} pair-eligible "
            f"fusions; observed {gene5}->{gene3}={run.summary['observed_count']}, "
            f"expected={run.summary['expected_count']:.2f}, Fisher p={fisher_display}"
        )
    else:
        click.echo(
            f"Analyzed {run.summary['total_fusions']} {run.gene_symbol} fusions; "
            f"mapped={run.summary['mapped_fusions']}, "
            f"in-frame={run.summary['in_frame_count']}, "
            f"domain-retained={run.summary['kinase_retained_count']}, "
            f"Fisher p={fisher_display}"
        )
    for warning in run.warnings:
        click.echo(f"Warning: {warning}", err=True)
    for kind, path in paths.items():
        click.echo(f"{kind}: {path}")


@main.command("analyze")
@click.argument("gene_symbol")
@click.argument("study_id")
@click.option(
    "--output-dir",
    type=click.Path(path_type=Path, file_okay=False),
    default=Path("runs"),
    show_default=True,
)
@click.option("--n-permutations", type=click.IntRange(min=1), default=1_000, show_default=True)
@click.option(
    "--pdf/--no-pdf",
    default=True,
    show_default=True,
    help="Also render a self-contained report.pdf for human reviewers.",
)
def analyze(
    gene_symbol: str, study_id: str, output_dir: Path, n_permutations: int, pdf: bool
) -> None:
    """Run all registered algorithms for a configured gene and live study."""
    try:
        run = run_analysis(gene_symbol, study_id, n_permutations=n_permutations)
        paths = write_outputs(
            run,
            output_dir,
            pdf=pdf,
            cli_args=[
                "analyze",
                gene_symbol,
                study_id,
                "--output-dir",
                str(output_dir),
                "--n-permutations",
                str(n_permutations),
            ],
        )
    except RealBenchmarkError as exc:
        raise click.ClickException(str(exc)) from None
    except Exception as exc:
        raise click.ClickException(f"Analysis failed: {type(exc).__name__}: {exc}") from None
    if run.is_gene_pair:
        gene5, gene3 = run.summary["gene_pair"]
        click.echo(
            f"Analyzed {run.summary['eligible_event_count']} {run.gene_symbol} pair-eligible "
            f"fusion events; observed {gene5}->{gene3}={run.summary['observed_count']}, "
            f"expected={run.summary['expected_count']:.2f}"
        )
    else:
        click.echo(
            f"Analyzed {run.summary['total_fusions']} {run.gene_symbol} fusions with "
            f"{len(run.results)} registered algorithms"
        )
    for warning in run.warnings:
        click.echo(f"Warning: {warning}", err=True)
    for kind, path in paths.items():
        click.echo(f"{kind}: {path}")


@main.command("cohort-scan")
@click.argument("study_id")
@click.option(
    "--output-dir",
    type=click.Path(path_type=Path, file_okay=False),
    default=Path("runs"),
    show_default=True,
)
@click.option(
    "--min-patients",
    "min_distinct_patients",
    type=click.IntRange(min=0),
    default=DEFAULT_MIN_DISTINCT_PATIENTS,
    show_default=True,
    help="Recurrence gate: minimum distinct patients with an SV record in a gene.",
)
@click.option("--n-permutations", type=click.IntRange(min=1), default=1_000, show_default=True)
@click.option(
    "--adaptive/--no-adaptive",
    default=True,
    show_default=True,
    help="Start permutation tests at a small budget and only escalate to "
    "--n-permutations when the small-N result is borderline.",
)
@click.option(
    "--n-permutations-small",
    type=click.IntRange(min=1),
    default=DEFAULT_N_PERMUTATIONS_SMALL,
    show_default=True,
)
@click.option(
    "--max-genes",
    type=click.IntRange(min=1),
    default=None,
    help="Cap the number of gated candidate genes actually analyzed (for a bounded run); "
    "reported before/after-gating counts are unaffected.",
)
@click.option(
    "--cache-dir",
    type=click.Path(path_type=Path, file_okay=False),
    default=None,
    help="Directory to cache batch Genome Nexus/Pfam-description lookups on disk. "
    "Defaults to <output-dir>/.cohort_scan_cache.",
)
@click.option(
    "--pdf/--no-pdf",
    default=True,
    show_default=True,
    help="Also render summary.pdf and full per-gene report.pdf files.",
)
@click.option(
    "--honorable-mention-count",
    type=click.IntRange(min=0),
    default=DEFAULT_HONORABLE_MENTION_COUNT,
    show_default=True,
    help="Size of the 'honorable mentions' highly ranked non-FDR-significant tier: the top N "
    "genes by raw Fisher p-value among genes that did NOT survive genome-wide FDR correction.",
)
def cohort_scan(
    study_id: str,
    output_dir: Path,
    min_distinct_patients: int,
    n_permutations: int,
    adaptive: bool,
    n_permutations_small: int,
    max_genes: int | None,
    cache_dir: Path | None,
    pdf: bool,
    honorable_mention_count: int,
) -> None:
    """Genome-wide fusion-hotspot scan: gate cohort-wide SV recurrence, run
    the full algorithm suite for every gated gene (auto-configuring genes
    with no curated config), and FDR-correct across every scanned gene."""
    try:
        result = run_cohort_scan(
            study_id,
            min_distinct_patients=min_distinct_patients,
            n_permutations=n_permutations,
            adaptive=adaptive,
            n_permutations_small=n_permutations_small,
            max_genes=max_genes,
            cache_dir=cache_dir or (output_dir / ".cohort_scan_cache"),
        )
        paths = write_cohort_scan_outputs(
            result, output_dir, pdf=pdf, honorable_mention_count=honorable_mention_count
        )
    except Exception as exc:
        raise click.ClickException(f"Cohort scan failed: {type(exc).__name__}: {exc}") from None

    ok_count = sum(1 for outcome in result.gene_outcomes if outcome.status == "ok")
    failed_count = len(result.gene_outcomes) - ok_count
    click.echo(
        f"{result.total_genes_before_gating} genes had SV records in {study_id}; "
        f"{result.genes_after_gating} passed the >= {min_distinct_patients}-patient gate "
        f"({result.curated_gene_count} curated, {result.auto_config_gene_count} auto-configured, "
        f"{result.unresolved_gene_count} unresolved)."
    )
    click.echo(f"Analyzed {ok_count} genes successfully, {failed_count} failed/skipped.")
    click.echo(f"FDR-significant genes (q<0.05): {len(result.significant_genes)}")
    for warning in result.warnings:
        click.echo(f"Warning: {warning}", err=True)
    for kind in (
        "run_directory",
        "summary_tsv",
        "summary_json",
        "summary_markdown",
        "summary_pdf",
        "manuscript_markdown",
        "manuscript_pdf",
    ):
        if kind in paths:
            click.echo(f"{kind}: {paths[kind]}")
    for gene_symbol, gene_paths in paths.get("gene_reports", {}).items():
        click.echo(f"gene_report[{gene_symbol}]: {gene_paths['run_directory']}")


if __name__ == "__main__":
    main()
