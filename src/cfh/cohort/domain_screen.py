"""Screen every fusion in a cohort for domains rebuilt at the junction.

``any_domain`` evaluates every Pfam family the two partners share.
``key_domain`` restricts that to each gene's key domain family, the family the
per-gene pipeline checks. The two modes differ only in which families are
examined. Output is one row per event and family, plus a summary.
"""

from __future__ import annotations

import csv
import json
from collections import Counter
from collections.abc import Iterable, Sequence
from pathlib import Path
from typing import Any, Callable

import pandas as pd

from cfh.cohort.auto_config import PfamDescriptionSource, select_key_domain
from cfh.genes.registry import load_gene_config
from cfh.mapping.domain_reconstitution import pfam_hit_for_interval
from cfh.mapping.genome_nexus_source import CanonicalTranscript
from cfh.mapping.junction_event import (
    JUNCTION_CALL_COLUMNS,
    MODE_ANY_DOMAIN,
    MODE_KEY_DOMAIN,
    NO_SHARED_FAMILY,
    READ_SUPPORT_COLUMNS,
    SEQUENCE_CONTEXT_COLUMNS,
    assess_event,
    read_support_columns,
)
from cfh.mapping.transcript_models import TranscriptModelSource, source_annotated_transcripts
from cfh.model.fusion_event import FusionEvent
from cfh.normalization.event_normalizer import normalize
from cfh.real_benchmark import _target_breakpoint

SCREEN_COLUMNS = (
    "sample_id",
    "patient_id",
    "fusion_name",
    "five_prime_gene",
    "three_prime_gene",
    "five_prime_breakpoint",
    "three_prime_breakpoint",
    "five_prime_chromosome",
    "three_prime_chromosome",
    *JUNCTION_CALL_COLUMNS,
    *READ_SUPPORT_COLUMNS,
    "source_annotation",
    *SEQUENCE_CONTEXT_COLUMNS,
)


def is_screenable_fusion(event: FusionEvent) -> bool:
    """A two-gene protein fusion with a known 5' and 3' partner."""
    five, three = event.Five_prime_gene, event.Three_prime_gene
    return bool(
        five
        and three
        and five != three
        and event.Is_protein_fusion
        and not event.Is_antisense
        and "fusion" in str(event.Event_info or "").lower()
    )


def key_domain_family(
    gene: str,
    canonical: CanonicalTranscript | None,
    description_source: PfamDescriptionSource | None = None,
) -> str | None:
    """The Pfam family of ``gene``'s key domain, curated first, then automatic."""
    if canonical is None:
        return None
    try:
        config = load_gene_config(gene)
    except FileNotFoundError:
        config = None
    if config is not None and config.key_domains:
        key = config.key_domains[0]
        if key.accession and key.accession.startswith("PF"):
            return key.accession
        hit = pfam_hit_for_interval(canonical.pfam_domains, key.start_aa, key.end_aa)
        if hit is not None:
            return hit.pfam_id
    auto = select_key_domain(canonical.pfam_domains, description_source)
    return auto.accession if auto is not None else None


def _partner_breakpoint(
    row: dict[str, Any], gene: str, span: tuple[int, int] | None
) -> tuple[int | None, str | None]:
    if span is None:
        return None, f"{gene} transcript has no exon coordinates"
    try:
        return _target_breakpoint(row, gene, span), None
    except ValueError as exc:
        return None, str(exc)


def _site_chromosome(row: dict[str, Any], gene: str, position: int | None) -> str | None:
    """Chromosome of the site reporting ``position``, preferring ``gene``'s site."""
    if position is None:
        return None
    matches = []
    for site in ("Site1", "Site2"):
        try:
            same_position = int(row.get(f"{site}_Position")) == position  # type: ignore[arg-type]
        except (TypeError, ValueError):
            continue
        chromosome = _clean(row.get(f"{site}_Chromosome"))
        if same_position and chromosome is not None:
            matches.append((row.get(f"{site}_Hugo_Symbol") == gene, str(chromosome)))
    if not matches:
        return None
    named = [chromosome for is_gene, chromosome in matches if is_gene]
    candidates = set(named) or {chromosome for _, chromosome in matches}
    return candidates.pop() if len(candidates) == 1 else None


def _clean(value: Any) -> Any:
    return None if value is None or (isinstance(value, float) and pd.isna(value)) else value


def screen_cohort(
    raw_rows: pd.DataFrame,
    source: TranscriptModelSource,
    *,
    cohort: str,
    mode: str = MODE_ANY_DOMAIN,
    key_family: Callable[[str], str | None] | None = None,
    progress: Callable[[str], None] | None = None,
    passthrough_columns: Sequence[str] = (),
    keep_unshared_rows: bool = True,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Screen ``raw_rows`` (see ``sv_parser``); returns ``(rows, summary)``.

    ``key_family`` is required in key-domain mode and maps a gene to its key
    family (see :func:`key_domain_family`). ``passthrough_columns`` are copied
    from each input row onto its output rows. With ``keep_unshared_rows``
    false, ``no_shared_family`` rows are counted in the summary but not
    returned, which bounds memory on very large tables.
    """
    if mode not in (MODE_ANY_DOMAIN, MODE_KEY_DOMAIN):
        raise ValueError(f"unknown screen mode {mode!r}")
    if mode == MODE_KEY_DOMAIN and key_family is None:
        raise ValueError("key-domain mode requires key_family")

    events = normalize(raw_rows, None, cohort)
    pairs: list[tuple[dict[str, Any], FusionEvent]] = [
        (row.where(pd.notnull(row), None).to_dict(), event)
        for (_, row), event in zip(raw_rows.iterrows(), events, strict=True)
        if is_screenable_fusion(event)
    ]
    source.prefetch_canonical(
        sorted({g for _, e in pairs for g in (e.Five_prime_gene, e.Three_prime_gene) if g})
    )

    rows: list[dict[str, Any]] = []
    class_counts: Counter[str] = Counter()
    for index, (row, event) in enumerate(pairs, start=1):
        five_gene, three_gene = event.Five_prime_gene, event.Three_prime_gene
        assert five_gene is not None and three_gene is not None
        annotation = row.get("Annotation")
        five = source.resolve(five_gene, source_annotated_transcripts(annotation, five_gene))
        three = source.resolve(three_gene, source_annotated_transcripts(annotation, three_gene))
        five_breakpoint, five_error = _partner_breakpoint(row, five_gene, five.span)
        three_breakpoint, three_error = _partner_breakpoint(row, three_gene, three.span)

        family_filter = None
        if mode == MODE_KEY_DOMAIN:
            assert key_family is not None
            allowed = {
                family
                for family in (key_family(five_gene), key_family(three_gene))
                if family is not None
            }
            family_filter = allowed.intersection

        calls = assess_event(
            five_prime_gene=five_gene,
            five_prime_breakpoint=five_breakpoint,
            three_prime_gene=three_gene,
            three_prime_breakpoint=three_breakpoint,
            annotation=annotation,
            genome_build=_clean(row.get("NCBI_Build")),
            source=source,
            family_filter=family_filter,
        )
        locus_note = "; ".join(error for error in (five_error, three_error) if error)
        for call in calls:
            class_counts[call["reconstitution_class"]] += 1
            if not keep_unshared_rows and call["reconstitution_class"] == NO_SHARED_FAMILY:
                continue
            if locus_note:
                call["no_product_reason"] = locus_note
                if call["reconstitution_class"] != NO_SHARED_FAMILY:
                    call["class_reason"] = locus_note
            rows.append(
                {
                    "sample_id": _clean(row.get("Sample_Id")),
                    "patient_id": event.Patient_id,
                    "fusion_name": f"{five_gene}::{three_gene}",
                    "five_prime_gene": five_gene,
                    "three_prime_gene": three_gene,
                    "five_prime_breakpoint": five_breakpoint,
                    "three_prime_breakpoint": three_breakpoint,
                    "five_prime_chromosome": _site_chromosome(row, five_gene, five_breakpoint),
                    "three_prime_chromosome": _site_chromosome(row, three_gene, three_breakpoint),
                    **call,
                    **read_support_columns(row),
                    "source_annotation": annotation,
                    **{column: _clean(row.get(column)) for column in passthrough_columns},
                }
            )
        if progress is not None and index % 200 == 0:
            progress(f"screened {index}/{len(pairs)} fusions")

    summary = summarize_screen(
        rows, screened_events=len(pairs), mode=mode, class_counts=class_counts
    )
    return rows, summary


def summarize_screen(
    rows: Iterable[dict[str, Any]],
    *,
    screened_events: int,
    mode: str,
    class_counts: Counter[str] | None = None,
) -> dict[str, Any]:
    """Counts by class; ``no_shared_family`` rows are counted but not listed.

    ``class_counts`` overrides the counts taken from ``rows``, for a screen
    that did not retain its ``no_shared_family`` rows.
    """
    rows = list(rows)
    classified = [row for row in rows if row["reconstitution_class"] != NO_SHARED_FAMILY]
    return {
        "mode": mode,
        "screened_fusion_events": screened_events,
        "events_with_a_shared_family": len(
            {
                (row["sample_id"], row["five_prime_gene"], row["three_prime_gene"])
                for row in classified
            }
        ),
        "family_rows": len(classified),
        "class_counts": dict(
            sorted(
                (
                    class_counts
                    if class_counts is not None
                    else Counter(row["reconstitution_class"] for row in rows)
                ).items()
            )
        ),
        "mechanism_counts": dict(
            sorted(Counter(str(row["junction_mechanism"]) for row in classified).items())
        ),
    }


def write_screen_outputs(
    rows: list[dict[str, Any]],
    summary: dict[str, Any],
    output_dir: str | Path,
    extra_columns: Sequence[str] = (),
) -> dict[str, Path]:
    """Write ``domain_screen.tsv`` and ``domain_screen_summary.json``.

    ``extra_columns`` (e.g. a screen's passthrough columns) follow
    :data:`SCREEN_COLUMNS`.
    """
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    table = destination / "domain_screen.tsv"
    with table.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[*SCREEN_COLUMNS, *extra_columns],
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)
    summary_path = destination / "domain_screen_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return {"domain_screen_tsv": table, "domain_screen_summary": summary_path}
