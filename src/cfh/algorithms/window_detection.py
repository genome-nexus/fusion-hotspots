"""Gene-agnostic, frequency/recurrence-only sliding-WINDOW detection.

Generalizes ``cutpoint_detection`` (WP11) from a single best-separating
*position* along a gene's protein axis to a variable sliding *window*
``[start, start + width]``: not every gene's functionally critical region
sits at a protein terminus the way BRAF/RET/ALK/NTRK1's kinase domains do
-- some future gene may need an internal region bounded on both sides
instead of a one-sided cutpoint. See :mod:`cfh.stats.window_scan` for the
underlying max-statistic permutation-correction scan this algorithm wraps
(the same pattern ``cutpoint_detection``/:mod:`cfh.stats.cutpoint_scan`
already use, extended from a 1-D to a 2-D ``(start, width)`` grid).

Like ``cutpoint_detection``, this uses no external oncogenicity label and
no OncoKB dependency: the outcome is domain-retention status (retained vs.
lost/disrupted), already computed per-event by the domain-retention
algorithm and stored on ``FusionFeature.Domain_retention_flags``. It runs
fully offline; an optional ``GenomeNexusClient`` may be passed in ``params``
purely to compare the inferred window against Pfam domain boundaries, the
same as ``cutpoint_detection`` already supports.

Mapping-sensitivity caveat: intronic breakpoints are commonly clamped to
the nearest exon boundary rather than mapped to an exact genomic position
(see the ``directional-intronic-breakpoint-snapping`` fix this module's
branch is built on). A window built substantially out of clamped/
approximate positions can look precise while resting on a mapping
approximation, so every result here reports what fraction of its
constituent events had a known-exact vs. clamped/approximate vs. unknown
mapping. That per-event quality signal is not part of ``FusionFeature`` --
it is supplied by the caller via ``params["mapping_sensitivity"]`` as an
``{event_id: is_intronic_breakpoint}`` mapping (``True`` = clamped/
approximate, ``False`` = exact, absent/``None`` = unknown), the same
optional, additive pattern ``cutpoint_detection`` already uses for
``domain_boundaries``/``genome_nexus_client``. Omitting it never crashes --
the caveat is simply reported as unavailable.
"""

from __future__ import annotations

from datetime import datetime, timezone
from typing import Optional

from cfh.algorithms.base import Algorithm
from cfh.algorithms.cutpoint_detection import _known_domain_boundaries, _nearest_boundary_comparison
from cfh.algorithms.registry import register
from cfh.genes.registry import GeneConfig
from cfh.model.algorithm_result import AlgorithmResult
from cfh.model.fusion_event import FusionEvent
from cfh.model.fusion_feature import FusionFeature
from cfh.stats.adaptive_permutation import is_borderline, resolve_permutation_budget
from cfh.stats.breakpoint_tests import gene_breakpoint_domain_status_event_records
from cfh.stats.window_scan import DEFAULT_MIN_EVENTS_PER_WINDOW, DEFAULT_WIDTHS, detect_window

ALGORITHM_NAME = "window_detection"
ALGORITHM_VERSION = "0.1.0"


def _mapping_sensitivity_summary(event_ids: list[str], mapping_sensitivity: Optional[dict]) -> dict:
    """Summarize what fraction of ``event_ids`` used a clamped/approximate
    breakpoint mapping vs. an exact one, so a human reviewer isn't misled by
    a result built substantially on a mapping approximation.

    ``mapping_sensitivity`` is an optional ``{event_id: is_intronic}``
    mapping (``True`` = clamped/approximate, ``False`` = exact,
    missing/``None`` = unknown). Never crashes on a missing or partial
    mapping -- unresolved events are simply counted as unknown, and
    ``available`` is ``False`` only when no mapping was supplied at all.
    """
    if not mapping_sensitivity:
        return {
            "available": False,
            "n_events": len(event_ids),
            "n_exact": 0,
            "n_clamped_or_approximate": 0,
            "n_unknown": len(event_ids),
            "fraction_clamped_or_approximate": None,
        }

    n_exact = 0
    n_clamped = 0
    n_unknown = 0
    for event_id in event_ids:
        flag = mapping_sensitivity.get(event_id)
        if flag is True:
            n_clamped += 1
        elif flag is False:
            n_exact += 1
        else:
            n_unknown += 1

    known = n_exact + n_clamped
    return {
        "available": True,
        "n_events": len(event_ids),
        "n_exact": n_exact,
        "n_clamped_or_approximate": n_clamped,
        "n_unknown": n_unknown,
        "fraction_clamped_or_approximate": (n_clamped / known) if known else None,
    }


@register(ALGORITHM_NAME)
class WindowDetectionAlgorithm(Algorithm):
    """Scan a gene's observed breakpoints for the protein-region *window*
    that best separates domain-retained from lost/disrupted events, with a
    permutation-corrected empirical p-value for having scanned many
    candidate ``(start, width)`` windows.

    Expected ``params`` keys (all optional):
        seed (int): permutation RNG seed, default 42.
        n_permutations (int): number of label permutations, default 10000.
        widths (list[int]): predeclared window widths (aa) to scan, default
            :data:`~cfh.stats.window_scan.DEFAULT_WIDTHS` (25/50/100/200).
        min_events_per_window (int): minimum event count required on both
            sides of a candidate window's split, default
            :data:`~cfh.stats.window_scan.DEFAULT_MIN_EVENTS_PER_WINDOW`.
        mapping_sensitivity (dict[str, bool | None]): optional
            ``{event_id: is_intronic_breakpoint}`` map used only to report
            the mapping-sensitivity caveat (see module docstring); never
            required.
        domain_boundaries (list[int]): known domain boundary aa positions to
            compare the inferred window's edges against.
        genome_nexus_client (GenomeNexusClient): optional fallback source
            for ``domain_boundaries``, same as ``cutpoint_detection``.
        uniprot_source: optional override for the UniProt fallback used by
            ``genome_nexus_client``-based boundary resolution (mainly for
            tests); ignored otherwise.

    Gracefully no-ops (``determinable: False`` with a ``reason``, never a
    crash) whenever the gene/data shape doesn't support a window scan --
    too few mapped events, a single distinct breakpoint position, a single
    outcome class, or no candidate window surviving the
    ``min_events_per_window`` guard -- exactly mirroring
    ``cutpoint_detection``'s no-op contract.
    """

    def run(
        self,
        events: list[FusionEvent],
        features: list[FusionFeature],
        gene_config: GeneConfig,
        params: dict,
    ) -> AlgorithmResult:
        params = params or {}
        seed = params.get("seed", 42)
        widths = tuple(params.get("widths", DEFAULT_WIDTHS))
        min_events_per_window = int(
            params.get("min_events_per_window", DEFAULT_MIN_EVENTS_PER_WINDOW)
        )
        budget = resolve_permutation_budget(params, default_full_n=10_000)

        records = gene_breakpoint_domain_status_event_records(events, features, gene_config)
        event_ids = [event_id for event_id, _, _ in records]
        positions = [position for _, position, _ in records]
        statuses = [status for _, _, status in records]

        n_permutations = budget["full_n"] if not budget["adaptive"] else budget["small_n"]
        scan_result = detect_window(
            positions,
            statuses,
            event_ids,
            seed=seed,
            n_permutations=n_permutations,
            widths=widths,
            min_events_per_window=min_events_per_window,
        )
        escalated = False
        if (
            budget["adaptive"]
            and scan_result["determinable"]
            and is_borderline(
                scan_result["corrected_p_value"],
                threshold=budget["threshold"],
                factor=budget["factor"],
            )
        ):
            escalated = True
            n_permutations = budget["full_n"]
            scan_result = detect_window(
                positions,
                statuses,
                event_ids,
                seed=seed,
                n_permutations=n_permutations,
                widths=widths,
                min_events_per_window=min_events_per_window,
            )

        warnings: list[str] = []
        boundary_comparison = None
        best_window_mapping_sensitivity = None
        mapping_sensitivity_param = params.get("mapping_sensitivity")
        overall_mapping_sensitivity = _mapping_sensitivity_summary(
            event_ids, mapping_sensitivity_param
        )
        if not overall_mapping_sensitivity["available"]:
            warnings.append(
                "mapping-sensitivity information was not supplied; the fraction of "
                "events using a clamped/approximate breakpoint position vs. an "
                "exact-genomic-mapped one could not be computed"
            )

        if scan_result["determinable"]:
            boundaries = _known_domain_boundaries(gene_config, params)
            best_window = scan_result["best_window"]
            if boundaries:
                boundary_comparison = {
                    "start_aa": _nearest_boundary_comparison(best_window["start_aa"], boundaries),
                    "end_aa": _nearest_boundary_comparison(best_window["end_aa"], boundaries),
                }
            best_window_event_ids = next(
                row["event_ids_inside"]
                for row in scan_result["scan"]
                if row["start_aa"] == best_window["start_aa"]
                and row["end_aa"] == best_window["end_aa"]
                and row["width_aa"] == best_window["width_aa"]
            )
            best_window_mapping_sensitivity = _mapping_sensitivity_summary(
                list(best_window_event_ids), mapping_sensitivity_param
            )
        else:
            warnings.append(scan_result["reason"])

        summary = {
            "determinable": scan_result["determinable"],
            "reason": scan_result["reason"],
            "n_events_analyzed": len(records),
            "widths_tested_aa": list(widths),
            "min_events_per_window": min_events_per_window,
            "n_candidate_windows": len(scan_result["scan"]),
            "n_distinct_candidate_windows_by_event_mask": len(scan_result["top_windows"]),
            "best_window": scan_result["best_window"],
            "observed_statistic_neg_log10_p": scan_result["observed_statistic"],
            "observed_p_value": scan_result["observed_p_value"],
            "observed_odds_ratio": scan_result["observed_odds_ratio"],
            "corrected_p_value": scan_result["corrected_p_value"],
            "known_domain_boundary_comparison": boundary_comparison,
            "mapping_sensitivity": overall_mapping_sensitivity,
            "best_window_mapping_sensitivity": best_window_mapping_sensitivity,
        }
        if budget["adaptive"]:
            summary["adaptive_permutations"] = {
                "enabled": True,
                "n_permutations_small": budget["small_n"],
                "n_permutations_full": budget["full_n"],
                "n_permutations_used": n_permutations,
                "escalated_to_full": escalated,
            }

        return AlgorithmResult(
            Algorithm=ALGORITHM_NAME,
            Algorithm_version=ALGORITHM_VERSION,
            Parameters={
                "seed": seed,
                "n_permutations": n_permutations,
                "widths": list(widths),
                "min_events_per_window": min_events_per_window,
                "adaptive": budget["adaptive"],
            },
            Summary=summary,
            Tables={
                "window_scan": scan_result["scan"],
                "top_windows": scan_result["top_windows"],
            },
            Warnings=warnings,
            Created_at=datetime.now(timezone.utc),
        )
