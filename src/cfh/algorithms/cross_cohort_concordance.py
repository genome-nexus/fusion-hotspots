"""Offline cross-cohort analysis; deliberately outside the live event registry."""

from __future__ import annotations

import json
import math
from dataclasses import asdict
from pathlib import Path
from typing import Any

from cfh.stats.cochran_mantel_haenszel import cochran_mantel_haenszel


def compare_cohort_runs(
    run_artifacts: list[Path], *, require_patient_disjoint: bool = False
) -> dict[str, Any]:
    """Read one gene's saved summary tables, preserving counts and provenance.

    Rows are retained/not retained; columns are in-frame/other. Input runs
    must use compatible domain definitions and counting units. Sample-ID
    overlap is reported, not deduplicated: aggregate tables cannot be repaired
    without changing the requested observations. Absence of matching IDs does
    not establish independence (patient IDs may be unavailable).
    """
    if len(run_artifacts) < 2:
        raise ValueError("Provide at least two runs from distinct cohorts")
    strata: list[dict[str, Any]] = []
    genes: set[str] = set()
    studies: set[str] = set()
    sample_sets: list[set[str]] = []
    patient_sets: list[set[str]] = []
    patient_ids_complete: list[bool] = []
    patient_namespaces: list[str | None] = []
    for artifact in run_artifacts:
        path = artifact / "results.json" if artifact.is_dir() else artifact
        try:
            data = json.loads(path.read_text())
            gene, study = data["gene_symbol"], data["study_id"]
            if not isinstance(gene, str) or not gene.strip():
                raise ValueError("Missing gene symbol")
            if not isinstance(study, str) or not study.strip():
                raise ValueError("Missing study ID")
            table = data["summary"]["frame_domain_contingency_table"]
            result = cochran_mantel_haenszel([table])
            events = data.get("events", [])
            samples = {e["sample_id"] for e in events if e.get("sample_id")}
            patients = {e["patient_id"] for e in events if e.get("patient_id")}
            complete = bool(events) and all(e.get("patient_id") for e in events)
            namespace = data.get("patient_id_namespace")
        except (OSError, ValueError, KeyError, TypeError, AttributeError) as exc:
            raise ValueError(f"Invalid run artifact {path}: {exc}") from exc
        genes.add(gene.upper())
        if study in studies:
            raise ValueError(f"Duplicate cohort: {study}")
        studies.add(study)
        sample_sets.append(samples)
        patient_sets.append(patients)
        patient_ids_complete.append(complete)
        patient_namespaces.append(namespace if isinstance(namespace, str) and namespace else None)
        a, b = table[0]
        c, d = table[1]
        strata.append(
            {
                "study_id": study,
                "source": str(path),
                "frame_domain_contingency_table": table,
                "key_domains": data["summary"].get("key_domains"),
                "event_count": len(events),
                "patient_id_coverage": "complete" if complete else "incomplete_or_unavailable",
                "patient_id_namespace": patient_namespaces[-1],
                "distinct_observed_patient_count": len(patients),
                "informative": bool(result.informative_strata),
                "association_direction": (
                    "positive" if a * d > b * c else "negative" if a * d < b * c else "neutral"
                ),
            }
        )
    if len(genes) != 1:
        raise ValueError("All runs must describe the same gene")
    overlaps = [
        {
            "study_ids": [strata[i]["study_id"], strata[j]["study_id"]],
            "shared_sample_id_count": len(sample_sets[i] & sample_sets[j]),
        }
        for i in range(len(strata))
        for j in range(i + 1, len(strata))
        if sample_sets[i] & sample_sets[j]
    ]
    patient_overlaps = [
        {
            "study_ids": [strata[i]["study_id"], strata[j]["study_id"]],
            "shared_patient_id_count": len(patient_sets[i] & patient_sets[j]),
        }
        for i in range(len(strata))
        for j in range(i + 1, len(strata))
        if patient_sets[i] & patient_sets[j]
    ]
    if require_patient_disjoint:
        if not all(patient_ids_complete):
            raise ValueError("Patient-disjoint validation requires patient IDs for every event")
        if None in patient_namespaces or len(set(patient_namespaces)) != 1:
            raise ValueError("Patient-disjoint validation requires a shared patient ID namespace")
        if patient_overlaps:
            raise ValueError("Patient-disjoint validation failed: shared patient IDs")
        if overlaps:
            raise ValueError("Patient-disjoint validation failed: shared sample IDs")
    pooled = cochran_mantel_haenszel([s["frame_domain_contingency_table"] for s in strata])
    cmh = asdict(pooled)
    # Keep machine-readable output strict JSON, including complete separation.
    if pooled.common_odds_ratio is not None and math.isinf(pooled.common_odds_ratio):
        cmh["common_odds_ratio"] = "infinity"
    significant = pooled.p_value is not None and pooled.p_value < 0.05
    warnings = [
        "CMH tests pooled conditional association, not equality of cohort odds ratios "
        "or significance in every cohort.",
        "Inference assumes independent observations within and between cohorts and compatible "
        "domain definitions/counting units; absence of shared sample IDs does not verify this.",
    ]
    if pooled.informative_strata < pooled.total_strata:
        warnings.append("Zero-variance strata contribute no evidence to the pooled test.")
    if overlaps:
        warnings.append(
            "Shared sample IDs violate between-cohort independence: the nominal p-value "
            "cannot establish independent cross-cohort replication. Tables were not deduplicated."
        )
    if patient_overlaps:
        warnings.append("Shared patient IDs violate between-cohort independence.")
    if not all(patient_ids_complete):
        warnings.append(
            "Patient ID coverage is incomplete; patient-disjoint validation is unknown."
        )
    elif (
        not patient_overlaps
        and not overlaps
        and len(set(patient_namespaces)) == 1
        and None not in patient_namespaces
    ):
        warnings.append(
            "Observed patient IDs are disjoint within supplied events and a declared shared "
            "namespace; cohort provenance and counting units still require review."
        )
    else:
        warnings.append("Patient-disjoint validation is not established by supplied metadata.")
    return {
        "gene_symbol": next(iter(genes)),
        "method": "CMH chi-square(1), no continuity correction",
        "row_labels": ["domain retained", "domain not retained"],
        "column_labels": ["in-frame protein fusion", "other"],
        "strata": strata,
        "cmh": cmh,
        "alpha": 0.05,
        "nominal_pooled_significant": significant,
        "all_cohorts_positive_and_informative": all(
            s["informative"] and s["association_direction"] == "positive" for s in strata
        ),
        "sample_overlaps": overlaps,
        "patient_overlaps": patient_overlaps,
        "patient_disjoint_observed_events": (
            all(patient_ids_complete)
            and not patient_overlaps
            and not overlaps
            and None not in patient_namespaces
            and len(set(patient_namespaces)) == 1
        ),
        "warnings": warnings,
        "interpretation": (
            "Nominal pooled association is significant at alpha=0.05. "
            if significant
            else "Nominal pooled association is not significant/testable at alpha=0.05. "
        )
        + "CMH alone does not confirm concordance across all cohorts; see stratum information "
        "and independence limitations.",
    }
