"""Assemble every metric and the heuristic tier for one chimera record."""

from __future__ import annotations

from typing import Any

from pyhmmer.plan7 import HMM

from cfh.domain_function.afdb import ParentStructure
from cfh.domain_function.metrics import (
    JUNCTION_HALF_WINDOW,
    NOT_ASSESSED,
    align_parents,
    junction_structure_report,
    key_residue_report,
    map_structure,
    place_chimera,
    schema_report,
)
from cfh.domain_function.pfam_hmm import (
    best_overlapping,
    gathering_thresholds,
    overlapping_hits,
    residues_in_inserts,
    scan,
    union_match_coverage,
)
from cfh.domain_function.record import ChimeraRecord
from cfh.domain_function.tiers import TIER_LABEL, assign_tier, disagreements
from cfh.domain_function.uniprot import UniProtEntry


def hmm_report(
    record: ChimeraRecord, hmm: HMM | None, junction_positions: list[int]
) -> dict[str, Any]:
    """HMMER metrics for the flanked chimera and parent windows.

    ``junction_positions`` are chimera positions (1-based, whole protein) whose
    insert-state placement is counted.
    """
    if hmm is None:
        return {"hmm_status": NOT_ASSESSED}
    seq_ga, dom_ga = gathering_thresholds(hmm)
    c_start, c_end = record.window()
    f_start, f_end = record.five.window()
    t_start, t_end = record.three.window()
    results = scan(
        hmm,
        [
            ("chimera", record.chimeric_protein[c_start - 1 : c_end]),
            ("five", record.five.protein[f_start - 1 : f_end]),
            ("three", record.three.protein[t_start - 1 : t_end]),
        ],
    )
    chimera_hit, chimera_hits = best_overlapping(
        results["chimera"], record.domain_start - c_start + 1, record.domain_end - c_start + 1
    )
    five_hit, _ = best_overlapping(
        results["five"],
        record.five.domain_start - f_start + 1,
        record.five.domain_end - f_start + 1,
    )
    three_hit, _ = best_overlapping(
        results["three"],
        record.three.domain_start - t_start + 1,
        record.three.domain_end - t_start + 1,
    )
    report: dict[str, Any] = {
        "hmm_status": "scored",
        "hmm_name": hmm.name.decode() if isinstance(hmm.name, bytes) else str(hmm.name),
        "hmm_accession": (
            hmm.accession.decode() if isinstance(hmm.accession, bytes) else str(hmm.accession)
        ),
        "hmm_length": hmm.M,
        "ga_sequence": seq_ga,
        "ga_domain": dom_ga,
        "chimera_domain_hits": chimera_hits,
    }
    regions = {
        "chimera": (
            results["chimera"],
            record.domain_start - c_start + 1,
            record.domain_end - c_start + 1,
        ),
        "parent5": (
            results["five"],
            record.five.domain_start - f_start + 1,
            record.five.domain_end - f_start + 1,
        ),
        "parent3": (
            results["three"],
            record.three.domain_start - t_start + 1,
            record.three.domain_end - t_start + 1,
        ),
    }
    for label, (result, start, end) in regions.items():
        hits = overlapping_hits(result, start, end)
        coverage = union_match_coverage(hits)
        report[f"{label}_hits_in_domain"] = len(hits)
        report[f"{label}_match_coverage_union"] = (
            round(coverage, 3) if coverage is not None else None
        )
    for label, hit in (("chimera", chimera_hit), ("parent5", five_hit), ("parent3", three_hit)):
        report[f"{label}_domain_bits"] = round(hit.score, 2) if hit else None
        report[f"{label}_i_evalue"] = float(f"{hit.i_evalue:.3g}") if hit else None
        report[f"{label}_ga_margin"] = round(hit.score - dom_ga, 2) if hit else None
        report[f"{label}_passes_ga"] = (hit.score >= dom_ga) if hit else False
        report[f"{label}_match_coverage"] = round(hit.match_state_coverage, 3) if hit else None
        report[f"{label}_hmm_span"] = f"{hit.hmm_from}-{hit.hmm_to}" if hit else ""
    for label, start, hit in (
        ("chimera", c_start, chimera_hit),
        ("parent5", f_start, five_hit),
        ("parent3", t_start, three_hit),
    ):
        report[f"{label}_ali_span"] = (
            f"{hit.ali_from + start - 1}-{hit.ali_to + start - 1}" if hit else ""
        )
    for label, key in (("chimera", "chimera"), ("parent5", "five"), ("parent3", "three")):
        result = results[key]
        report[f"{label}_sequence_bits"] = round(result.sequence_score, 2) if result else None
    if chimera_hit is not None:
        for parent, hit in (("5p", five_hit), ("3p", three_hit)):
            report[f"delta_bits_vs_{parent}"] = (
                round(chimera_hit.score - hit.score, 2) if hit else None
            )
        junction = record.junction_index - c_start + 1
        report["chimera_hit_spans_junction"] = (
            chimera_hit.ali_from <= junction and junction + 1 <= chimera_hit.ali_to
        )
        window_positions = [p - c_start + 1 for p in junction_positions]
        report["junction_residues"] = len(window_positions)
        report["junction_residues_in_inserts"] = residues_in_inserts(chimera_hit, window_positions)
        report["junction_residues_outside_alignment"] = sum(
            1 for p in window_positions if not chimera_hit.ali_from <= p <= chimera_hit.ali_to
        )
    else:
        for key in ("delta_bits_vs_5p", "delta_bits_vs_3p", "chimera_hit_spans_junction"):
            report[key] = None
        report["junction_residues"] = len(junction_positions)
        report["junction_residues_in_inserts"] = None
        report["junction_residues_outside_alignment"] = None
    for parent, key in (("5p", "parent5"), ("3p", "parent3")):
        chimera_seq = report["chimera_sequence_bits"]
        parent_seq = report[f"{key}_sequence_bits"]
        report[f"hmmvar_delta_vs_{parent}"] = (
            round(chimera_seq - parent_seq, 2)
            if chimera_seq is not None and parent_seq is not None
            else None
        )
    return report


def assess_record(
    record: ChimeraRecord,
    hmm: HMM | None,
    entries: dict[str, UniProtEntry | None],
    structures: dict[str, ParentStructure | None],
    structure_problems: dict[str, str] | None = None,
) -> dict[str, Any]:
    """Every metric plus the HEURISTIC tier for ``record``.

    ``entries`` / ``structures`` / ``structure_problems`` are keyed ``five`` /
    ``three``; ``structure_problems`` holds a fetch error when AFDB could not be
    reached, so a failed download is reported as such rather than as "no model".
    """
    problems = structure_problems or {}
    parents = align_parents(record)
    placed = place_chimera(record, parents)
    junction_positions = sorted(
        set(
            range(
                record.junction_index - JUNCTION_HALF_WINDOW + 1,
                record.junction_index + JUNCTION_HALF_WINDOW + 1,
            )
        )
        | set(placed.inserted)
    )
    metrics: dict[str, Any] = {
        "aligned_inserted_residues": len(placed.inserted),
        "aligned_inserted_sequence": "".join(
            record.chimeric_protein[p - 1] for p in placed.inserted
        ),
        "aligned_deleted_columns": len(placed.deleted_columns),
        "aligned_deleted_5p": len(placed.deleted_residues["five"]),
        "aligned_deleted_3p": len(placed.deleted_residues["three"]),
    }
    metrics.update(hmm_report(record, hmm, junction_positions))
    metrics.update(key_residue_report(record, parents, placed, entries))
    mapped = {}
    for which, parent in (("five", record.five), ("three", record.three)):
        prefix = "5p" if which == "five" else "3p"
        result, status, detail = map_structure(
            parent, entries.get(which), structures.get(which), problems.get(which, "")
        )
        mapped[which] = result
        metrics[f"structure_{prefix}"] = (
            f"{result.structure.model_id} v{result.structure.model_version}"
            if result
            else NOT_ASSESSED
        )
        metrics[f"structure_status_{prefix}"] = status
        metrics[f"structure_note_{prefix}"] = detail
    metrics.update(junction_structure_report(record, placed, mapped))
    metrics.update(schema_report(record, parents, placed, mapped))
    tier, rules, missing = assign_tier(metrics)
    metrics["tier_kind"] = TIER_LABEL
    metrics["tier"] = tier
    metrics["tier_rules"] = ";".join(rules)
    metrics["tier_missing_inputs"] = ";".join(missing)
    metrics["metric_disagreements"] = ";".join(disagreements(metrics, tier))
    return metrics
