"""Gene-agnostic genomic-position recurrence, distinct from protein-position
cutpoint detection.

``cutpoint_detection``/``window_detection`` scan a gene's breakpoints in
*protein* coordinate space (``FusionFeature.Junction_position_aa``). Two
independent events that land on the same intron get clamped onto the same
exon boundary and therefore look identical in protein space, even when their
real DNA-level breakpoints sit tens of kilobases apart inside that intron --
or, conversely, genuinely share one exact genomic coordinate. Protein-space
recurrence alone cannot tell those two cases apart.

This algorithm answers that question directly by binning/clustering a gene's
real, committed breakpoints in *genomic* coordinate space and reporting
where they pile up, then cross-referencing that against the protein-position
picture: for every protein junction position shared by two or more events,
it reports whether those events' genomic breakpoints are the exact same DNA
coordinate (a real hotspot) or scattered across a wider genomic span
(independent intronic breakpoints that only look recurrent once mapped/
clamped onto the same protein/exon boundary).

Deliberately simple by design: a recurrence table (exact-position counts and
fixed-width genomic bins) plus one descriptive cross-reference note. This is
not a scan-statistic/permutation search over genomic windows -- that is
future scope, mirroring ``window_detection``'s relationship to
``cutpoint_detection`` but one step earlier; nothing here produces a
p-value.

No new network call and no new core-model field is read here beyond what the
caller already resolved: genomic breakpoints are supplied explicitly via
``params["genomic_breakpoints"]`` as ``{event_id: {"chromosome", "position",
"build"}}``, the same caller-supplied-sidecar-dict pattern
``window_detection`` already uses for ``mapping_sensitivity``. This is
deliberate, not a shortcut: a raw SV row's ``Site1``/``Site2`` labels are not
a reliable indicator of which site belongs to the target gene (see
``real_benchmark._target_breakpoint``'s docstring), so this algorithm always
uses the same already-locus-validated genomic position that produced the
event's ``Junction_position_aa`` in the first place -- anything else would
make the genomic-vs-protein cross-reference below incoherent.

Build handling is explicit and honest (never silently pooled): among records
with a nonempty reported build, only the most frequent build is pooled. Empty
or unknown build provenance is pooled only when every usable record is
unknown, and is never presented as a validated assembly. All exclusions and
unknown counts are surfaced in ``Summary``/``Warnings``.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from datetime import datetime, timezone
from typing import Any, Optional

from cfh.algorithms.base import Algorithm
from cfh.algorithms.composite_score import _failed, _results_by_name
from cfh.algorithms.registry import register
from cfh.genes.registry import GeneConfig
from cfh.model.algorithm_result import AlgorithmResult
from cfh.model.fusion_event import FusionEvent
from cfh.model.fusion_feature import FusionFeature

ALGORITHM_NAME = "genomic_position_recurrence"
ALGORITHM_VERSION = "0.1.0"
DEFAULT_BIN_SIZE_BP = 1_000


def _target_features(features: list[FusionFeature], gene_config: GeneConfig):
    target = (gene_config.gene_symbol or "").upper()
    return (feature for feature in features if feature.Gene.upper() == target)


def _genomic_records(
    events: list[FusionEvent],
    features: list[FusionFeature],
    gene_config: GeneConfig,
    genomic_breakpoints: dict[str, Any],
) -> tuple[list[dict], int]:
    """Join caller-supplied genomic breakpoints onto this gene's mapped events.

    Only events that are both (a) among this gene's already-mapped
    ``FusionFeature`` records and (b) present in ``genomic_breakpoints`` with
    a known chromosome and position are returned. Anything else (an event
    with no genomic-breakpoint entry, or one missing chromosome/position)
    is counted in the second return value rather than silently vanishing.
    An empty build is retained as unknown provenance; it is pooled only when
    every usable record has unknown provenance, and is never presented as a
    reported reference build.
    """
    event_ids = {event.Event_id for event in events}
    target_event_ids = {
        feature.Event_id
        for feature in _target_features(features, gene_config)
        if feature.Event_id in event_ids
    }
    records: list[dict] = []
    for event_id in sorted(target_event_ids):
        info = genomic_breakpoints.get(event_id)
        chromosome = info.get("chromosome") if isinstance(info, dict) else None
        position = info.get("position") if isinstance(info, dict) else None
        build = info.get("build") if isinstance(info, dict) else None
        if chromosome is None or position is None:
            continue
        records.append(
            {
                "event_id": event_id,
                "chromosome": str(chromosome),
                "position": int(position),
                "build": str(build).strip() if build is not None and str(build).strip() else None,
            }
        )
    skipped = len(target_event_ids) - len(records)
    return records, skipped


def _select_pooled_build(
    records: list[dict],
) -> tuple[Optional[str], list[dict], dict[str | None, int]]:
    """Pick the single reference build to pool genomic positions within.

    Coordinates from two different genome assemblies are not comparable
    numbers even on the same-numbered chromosome, so a record reporting a
    build other than this run's most frequent nonempty reported build is
    excluded from pooling here -- never silently mixed in. ``build_counts``
    names exactly what was seen (including builds excluded), so the exclusion is always visible in
    ``Summary``, not just implied by a smaller count.
    """
    build_counts = Counter(record["build"] for record in records)
    if not build_counts:
        return None, [], {}
    known_build_counts = Counter({build: count for build, count in build_counts.items() if build})
    if not known_build_counts:
        # Unknown-provenance coordinates can be described among themselves,
        # but are never labelled as a reported reference build.
        return None, [record for record in records if record["build"] is None], dict(build_counts)
    dominant_build = known_build_counts.most_common(1)[0][0]
    # Do not mix unknown provenance with reported-build coordinates, even if
    # unknown records happen to be the numerical majority.
    pooled = [record for record in records if record["build"] == dominant_build]
    return dominant_build, pooled, dict(build_counts)


def _bin_start(position: int, bin_size_bp: int) -> int:
    return (position // bin_size_bp) * bin_size_bp


def _genomic_bin_recurrence_table(records: list[dict], bin_size_bp: int) -> list[dict]:
    bins: dict[tuple[str, int], list[str]] = defaultdict(list)
    for record in records:
        key = (record["chromosome"], _bin_start(record["position"], bin_size_bp))
        bins[key].append(record["event_id"])
    table: list[dict[str, Any]] = [
        {
            "chromosome": chromosome,
            "bin_start": bin_start,
            "bin_end": bin_start + bin_size_bp,
            "n_events": len(event_ids),
            "event_ids": sorted(event_ids),
        }
        for (chromosome, bin_start), event_ids in bins.items()
    ]
    table.sort(key=lambda row: (-int(row["n_events"]), row["chromosome"], row["bin_start"]))
    return table


def _exact_position_recurrence_table(records: list[dict]) -> list[dict]:
    """Events sharing the *exact same* genomic coordinate (bin-independent).

    Only positions shared by two or more events are listed -- a singleton
    genomic position carries no recurrence signal.
    """
    counts: dict[tuple[str, int], list[str]] = defaultdict(list)
    for record in records:
        counts[(record["chromosome"], record["position"])].append(record["event_id"])
    table: list[dict[str, Any]] = [
        {
            "chromosome": chromosome,
            "position": position,
            "n_events": len(event_ids),
            "event_ids": sorted(event_ids),
        }
        for (chromosome, position), event_ids in counts.items()
        if len(event_ids) > 1
    ]
    table.sort(key=lambda row: (-int(row["n_events"]), row["chromosome"], row["position"]))
    return table


def _protein_position_by_event(
    features: list[FusionFeature], gene_config: GeneConfig
) -> dict[str, int]:
    return {
        feature.Event_id: feature.Junction_position_aa
        for feature in _target_features(features, gene_config)
        if feature.Junction_position_aa is not None
    }


def _protein_position_genomic_spread_table(
    records: list[dict], protein_position_by_event: dict[str, int]
) -> list[dict]:
    """For each protein junction position shared by >=2 events with a known
    genomic breakpoint, report whether those events' genomic breakpoints are
    one shared DNA coordinate (a real hotspot) or scattered across several
    distinct positions (independent breakpoints that only look recurrent
    because they map/clamp onto the same protein position). This is the
    direct genomic-vs-protein cross-reference the module docstring describes.
    """
    record_by_event = {record["event_id"]: record for record in records}
    events_by_aa: dict[int, list[str]] = defaultdict(list)
    for event_id, aa in protein_position_by_event.items():
        if event_id in record_by_event:
            events_by_aa[aa].append(event_id)

    rows: list[dict] = []
    for aa, event_ids in events_by_aa.items():
        if len(event_ids) < 2:
            continue
        genomic = [record_by_event[event_id] for event_id in event_ids]
        chromosomes = {record["chromosome"] for record in genomic}
        distinct_positions = {(record["chromosome"], record["position"]) for record in genomic}
        positions_only = [record["position"] for record in genomic]
        rows.append(
            {
                "protein_position_aa": aa,
                "n_events": len(event_ids),
                "n_distinct_genomic_positions": len(distinct_positions),
                "genomic_span_bp": (
                    max(positions_only) - min(positions_only) if len(chromosomes) == 1 else None
                ),
                "spans_multiple_chromosomes": len(chromosomes) > 1,
            }
        )
    rows.sort(key=lambda row: (-int(row["n_events"]), row["protein_position_aa"]))
    return rows


def _clustering_note(protein_spread_table: list[dict], cutpoint_summary: Optional[dict]) -> str:
    """Plain-language summary of whether genomic clustering differs from
    protein-position clustering, preferring the row at ``cutpoint_detection``'s
    inferred cutpoint (when available and determinable) so the note directly
    answers "is *this* protein-recurrence a real DNA hotspot or not", and
    falling back to the most-recurrent protein position otherwise.
    """
    if not protein_spread_table:
        return (
            "No protein-junction position was shared by two or more events with a known "
            "genomic breakpoint, so genomic-vs-protein clustering could not be compared."
        )
    row = protein_spread_table[0]
    if cutpoint_summary and cutpoint_summary.get("determinable"):
        cutpoint_aa = cutpoint_summary.get("inferred_cutpoint_aa")
        matched = next(
            (r for r in protein_spread_table if r["protein_position_aa"] == cutpoint_aa),
            None,
        )
        if matched is not None:
            row = matched
    if row["spans_multiple_chromosomes"]:
        return (
            f"The {row['n_events']} events sharing protein position "
            f"{row['protein_position_aa']} aa have genomic breakpoints on more than one "
            "chromosome; genomic clustering cannot be summarized as a single span."
        )
    if row["n_distinct_genomic_positions"] == 1:
        return (
            f"All {row['n_events']} events sharing protein position "
            f"{row['protein_position_aa']} aa share the exact same genomic breakpoint -- "
            "consistent with a real DNA-level recurrent breakpoint, not just a "
            "protein-coordinate coincidence."
        )
    return (
        f"The {row['n_events']} events sharing protein position {row['protein_position_aa']} aa "
        f"use {row['n_distinct_genomic_positions']} distinct genomic positions spanning "
        f"{row['genomic_span_bp']} bp -- the protein-position recurrence is broader than any "
        "single genomic hotspot, consistent with independent intronic breakpoints mapped (or "
        "clamped) onto the same protein/exon boundary rather than one shared DNA lesion."
    )


def _cutpoint_summary(params: dict) -> Optional[dict]:
    results_by_name = _results_by_name(params.get("algorithm_results"))
    cutpoint_result = results_by_name.get("cutpoint_detection")
    if cutpoint_result is None or _failed(cutpoint_result):
        return None
    return cutpoint_result.Summary


_EMPTY_TABLES: dict[str, list] = {
    "genomic_bin_recurrence": [],
    "exact_position_recurrence": [],
    "protein_position_genomic_spread": [],
}


@register(ALGORITHM_NAME)
class GenomicPositionRecurrenceAlgorithm(Algorithm):
    """Bin/cluster a gene's real, committed breakpoints by genomic (DNA)
    position and cross-reference that against ``cutpoint_detection``'s
    protein-position clustering result.

    Expected ``params`` keys:
        genomic_breakpoints (dict[str, dict], required for a determinable
            result): ``{event_id: {"chromosome": str, "position": int,
            "build": str | None}}`` for this gene's already-mapped events --
            the same genomic breakpoint (already locus-validated against the
            target gene, not a raw/possibly-mislabeled ``Site1``/``Site2``
            value) that produced each event's ``Junction_position_aa``.
            Missing entirely, or missing for every event, yields a graceful
            ``determinable: False`` result, never a crash.
        bin_size_bp (int, optional): fixed genomic bin width in base pairs
            for ``Tables["genomic_bin_recurrence"]``, default 1000.

    ``DEPENDS_ON = ("cutpoint_detection",)`` so its inferred protein cutpoint
    (when available) is used to select which protein-position row the
    ``genomic_vs_protein_clustering_note`` describes; a run that didn't also
    request ``cutpoint_detection`` still gets a note, just anchored on the
    most-recurrent protein position instead.

    Never pools genomic positions across differing reported reference-genome
    builds: only the most frequent nonempty reported build is pooled. Unknown
    provenance is pooled only when all usable records are unknown, and every
    exclusion is reported in both ``Summary`` and ``Warnings``. The aggregate
    ``n_events_excluded_other_build`` includes unknown-provenance exclusions;
    ``n_events_excluded_other_reported_build`` isolates other nonempty builds.
    """

    DEPENDS_ON = ("cutpoint_detection",)

    def run(
        self,
        events: list[FusionEvent],
        features: list[FusionFeature],
        gene_config: GeneConfig,
        params: dict,
    ) -> AlgorithmResult:
        params = params or {}
        bin_size_bp = int(params.get("bin_size_bp", DEFAULT_BIN_SIZE_BP))
        genomic_breakpoints = params.get("genomic_breakpoints") or {}
        warnings: list[str] = []
        event_ids = {event.Event_id for event in events}
        eligible_event_ids = {
            feature.Event_id
            for feature in _target_features(features, gene_config)
            if feature.Event_id in event_ids
        }

        if not genomic_breakpoints:
            warnings.append(
                "genomic-position recurrence was not computed: no genomic breakpoint "
                "positions were supplied for this run (params['genomic_breakpoints'])."
            )
            return AlgorithmResult(
                Algorithm=ALGORITHM_NAME,
                Algorithm_version=ALGORITHM_VERSION,
                Parameters={"bin_size_bp": bin_size_bp},
                Summary={
                    "determinable": False,
                    "reason": "no genomic breakpoint positions were supplied",
                    "n_events_analyzed": 0,
                    "n_events_excluded_missing_position": len(eligible_event_ids),
                    "reference_build": None,
                    "bin_size_bp": bin_size_bp,
                },
                Tables=dict(_EMPTY_TABLES),
                Warnings=warnings,
                Created_at=datetime.now(timezone.utc),
            )

        records, skipped_no_position = _genomic_records(
            events, features, gene_config, genomic_breakpoints
        )
        if not records:
            reason = (
                "no target-gene event had both a mapped breakpoint and a known genomic position"
            )
            warnings.append(reason)
            return AlgorithmResult(
                Algorithm=ALGORITHM_NAME,
                Algorithm_version=ALGORITHM_VERSION,
                Parameters={"bin_size_bp": bin_size_bp},
                Summary={
                    "determinable": False,
                    "reason": reason,
                    "n_events_analyzed": 0,
                    "n_events_excluded_missing_position": skipped_no_position,
                    "reference_build": None,
                    "bin_size_bp": bin_size_bp,
                },
                Tables=dict(_EMPTY_TABLES),
                Warnings=warnings,
                Created_at=datetime.now(timezone.utc),
            )

        build, pooled_records, build_counts = _select_pooled_build(records)
        excluded_other_build = len(records) - len(pooled_records)
        unknown_build_count = build_counts.get(None, 0)
        if unknown_build_count:
            if build is None:
                warnings.append(
                    f"{unknown_build_count} event(s) have unknown or empty reference-build "
                    "provenance; positions were grouped only with other unknown-build "
                    "records and reference_build remains null (not reported as an assembly)."
                )
            else:
                warnings.append(
                    f"{unknown_build_count} event(s) have unknown or empty reference-build "
                    f"provenance and were excluded from pooling with reported build {build!r}."
                )
        if len(build_counts) > 1:
            warnings.append(
                "Events reported more than one reference genome build "
                f"({build_counts}); only the most frequent nonempty reported build ({build!r}) "
                f"was pooled for genomic-position recurrence -- {excluded_other_build} event(s) "
                "with other or unknown provenance were excluded from pooling rather than "
                "silently mixed in."
            )

        genomic_bin_table = _genomic_bin_recurrence_table(pooled_records, bin_size_bp)
        exact_position_table = _exact_position_recurrence_table(pooled_records)
        protein_position_by_event = _protein_position_by_event(features, gene_config)
        protein_spread_table = _protein_position_genomic_spread_table(
            pooled_records, protein_position_by_event
        )
        note = _clustering_note(protein_spread_table, _cutpoint_summary(params))

        summary = {
            "determinable": True,
            "reason": None,
            "n_events_analyzed": len(pooled_records),
            "n_events_excluded_missing_position": skipped_no_position,
            "n_events_excluded_other_build": excluded_other_build,
            "n_events_excluded_other_reported_build": (
                excluded_other_build - (unknown_build_count if build is not None else 0)
            ),
            "n_events_unknown_build": unknown_build_count,
            "n_events_excluded_unknown_build": (unknown_build_count if build is not None else 0),
            "reference_build": build,
            "reference_build_counts": build_counts,
            "bin_size_bp": bin_size_bp,
            "n_distinct_genomic_positions": len(
                {(record["chromosome"], record["position"]) for record in pooled_records}
            ),
            "top_genomic_bin": genomic_bin_table[0] if genomic_bin_table else None,
            "n_positions_shared_by_multiple_events": len(exact_position_table),
            "genomic_vs_protein_clustering_note": note,
        }
        return AlgorithmResult(
            Algorithm=ALGORITHM_NAME,
            Algorithm_version=ALGORITHM_VERSION,
            Parameters={"bin_size_bp": bin_size_bp},
            Summary=summary,
            Tables={
                "genomic_bin_recurrence": genomic_bin_table,
                "exact_position_recurrence": exact_position_table,
                "protein_position_genomic_spread": protein_spread_table,
            },
            Warnings=warnings,
            Created_at=datetime.now(timezone.utc),
        )
