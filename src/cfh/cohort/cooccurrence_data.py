"""Cohort-wide inputs for :mod:`cfh.cohort.cooccurrence_discovery`.

Fetches, once per cohort scan: each sample's gene panel for the SV, mutation,
and discrete-CNA profiles; each panel's gene list; every panel gene's
mutated, amplified, and deep-deleted samples; and each sample's OncoTree
code. Mutations use cBioPortal's ``ID`` projection, so every imported
mutation record counts (cBioPortal study imports normally omit silent and
non-coding variants, but this module does not re-filter them).
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from typing import Any

import requests

from cfh.cohort.cooccurrence_discovery import SampleRecord
from cfh.ingestion import cbioportal_api

DEFAULT_GENE_BATCH_SIZE = 50
UNKNOWN_STRATUM = "UNKNOWN"
_CNA_TYPES = {2: "amplification", -2: "deep_deletion"}


@dataclass
class DiscoveryInputs:
    samples: list[SampleRecord]
    panel_genes: dict[str, set[str]]
    alterations: dict[tuple[str, str], set[str]]
    mutation_count_by_sample: Counter
    warnings: list[str] = field(default_factory=list)


def _post(session: requests.Session, url: str, body: dict, timeout: float) -> Any:
    response = session.post(url, json=body, timeout=timeout)
    response.raise_for_status()
    return response.json()


def fetch_panel_assignments(
    profile_id: str,
    sample_list_id: str,
    *,
    base_url: str,
    session: requests.Session,
    cache: dict,
    timeout: float = 120,
) -> dict[str, str | None]:
    """Sample -> gene panel for one profile; unprofiled samples map to ``None``.

    Shares ``cache`` keys with
    :func:`cfh.ingestion.cbioportal_api.fetch_gene_panel_eligibility`.
    """
    base = base_url.rstrip("/")
    key = ("gene_panel_data", base, profile_id, sample_list_id)
    if key not in cache:
        records = _post(
            session,
            f"{base}/molecular-profiles/{profile_id}/gene-panel-data/fetch",
            {"sampleListId": sample_list_id},
            timeout,
        )
        if not isinstance(records, list):
            raise ValueError("Gene-panel data must be a list")
        cache[key] = records
    assignments: dict[str, str | None] = {}
    for record in cache[key]:
        panel = record.get("genePanelId")
        profiled = record.get("profiled") is True and panel and panel != "NA"
        assignments[str(record["sampleId"])] = panel if profiled else None
    return assignments


def fetch_panel_gene_symbols(
    panel_id: str, *, base_url: str, session: requests.Session, cache: dict, timeout: float = 60
) -> dict[int, str]:
    """Entrez ID -> HUGO symbol for every gene on a panel."""
    base = base_url.rstrip("/")
    key = ("gene_panel_symbols", base, panel_id)
    if key not in cache:
        response = session.get(f"{base}/gene-panels/{panel_id}", timeout=timeout)
        response.raise_for_status()
        genes = (response.json() or {}).get("genes") or []
        cache[key] = {
            int(gene["entrezGeneId"]): str(gene["hugoGeneSymbol"]).upper()
            for gene in genes
            if gene.get("entrezGeneId") is not None and gene.get("hugoGeneSymbol")
        }
    return cache[key]


def load_discovery_inputs(
    study_id: str,
    *,
    sv_profile_id: str,
    mutation_profile_id: str | None,
    cna_profile_id: str | None,
    sample_list_id: str,
    base_url: str = cbioportal_api.DEFAULT_BASE_URL,
    session: requests.Session | None = None,
    panel_cache: dict | None = None,
    gene_batch_size: int = DEFAULT_GENE_BATCH_SIZE,
    timeout: float = 300,
) -> DiscoveryInputs:
    """Fetch every cohort-wide input co-occurrence discovery needs."""
    session = session or requests.Session()
    cache = {} if panel_cache is None else panel_cache
    base = base_url.rstrip("/")
    warnings: list[str] = []

    sv_panels = fetch_panel_assignments(
        sv_profile_id, sample_list_id, base_url=base, session=session, cache=cache
    )
    profile_panels: dict[str, dict[str, str | None]] = {}
    for kind, profile_id in (("mutation", mutation_profile_id), ("cna", cna_profile_id)):
        if profile_id is None:
            warnings.append(f"No {kind} profile is configured; {kind} comparators were skipped.")
            profile_panels[kind] = {}
            continue
        try:
            profile_panels[kind] = fetch_panel_assignments(
                profile_id, sample_list_id, base_url=base, session=session, cache=cache
            )
        except (requests.RequestException, ValueError) as exc:
            warnings.append(
                f"{kind} gene-panel coverage for {profile_id} was unavailable "
                f"({type(exc).__name__}); {kind} comparators were skipped."
            )
            profile_panels[kind] = {}

    panel_ids = {
        panel
        for assignments in (sv_panels, *profile_panels.values())
        for panel in assignments.values()
        if panel
    }
    symbols_by_panel = {
        panel: fetch_panel_gene_symbols(panel, base_url=base, session=session, cache=cache)
        for panel in sorted(panel_ids)
    }
    panel_genes = {panel: set(symbols.values()) for panel, symbols in symbols_by_panel.items()}
    symbol_by_entrez = {
        entrez: symbol
        for symbols in symbols_by_panel.values()
        for entrez, symbol in symbols.items()
    }

    def comparator_entrez_ids(kind: str) -> list[int]:
        panels = {panel for panel in profile_panels[kind].values() if panel}
        return sorted({entrez for panel in panels for entrez in symbols_by_panel[panel]})

    alterations: dict[tuple[str, str], set[str]] = {}
    mutation_counts: Counter = Counter()
    if profile_panels["mutation"] and mutation_profile_id:
        ids = comparator_entrez_ids("mutation")
        for start in range(0, len(ids), gene_batch_size):
            records = _post(
                session,
                f"{base}/molecular-profiles/{mutation_profile_id}/mutations/fetch?projection=ID",
                {
                    "sampleListId": sample_list_id,
                    "entrezGeneIds": ids[start : start + gene_batch_size],
                },
                timeout,
            )
            for record in records:
                symbol = symbol_by_entrez.get(record.get("entrezGeneId"))
                sample_id = record.get("sampleId")
                if symbol and sample_id:
                    alterations.setdefault((symbol, "mutation"), set()).add(sample_id)
                    mutation_counts[sample_id] += 1
    if profile_panels["cna"] and cna_profile_id:
        ids = comparator_entrez_ids("cna")
        for start in range(0, len(ids), gene_batch_size):
            records = _post(
                session,
                f"{base}/molecular-profiles/{cna_profile_id}/discrete-copy-number/fetch"
                "?projection=ID&discreteCopyNumberEventType=HOMDEL_AND_AMP",
                {
                    "sampleListId": sample_list_id,
                    "entrezGeneIds": ids[start : start + gene_batch_size],
                },
                timeout,
            )
            for record in records:
                symbol = symbol_by_entrez.get(record.get("entrezGeneId"))
                sample_id = record.get("sampleId")
                alteration = record.get("alteration")
                cna_type = _CNA_TYPES.get(alteration)
                if symbol and sample_id and cna_type:
                    alterations.setdefault((symbol, cna_type), set()).add(sample_id)

    response = session.get(
        f"{base}/studies/{study_id}/clinical-data",
        params={
            "clinicalDataType": "SAMPLE",
            "attributeId": "ONCOTREE_CODE",
            "pageSize": "10000000",
        },
        timeout=timeout,
    )
    response.raise_for_status()
    oncotree = {
        str(record["sampleId"]): str(record["value"])
        for record in response.json()
        if record.get("sampleId") and record.get("value")
    }
    missing_stratum = sum(1 for sample_id in sv_panels if sample_id not in oncotree)
    if missing_stratum:
        warnings.append(
            f"{missing_stratum} samples had no ONCOTREE_CODE and form one "
            f"{UNKNOWN_STRATUM} stratum."
        )
    samples = [
        SampleRecord(
            sample_id=sample_id,
            stratum=oncotree.get(sample_id, UNKNOWN_STRATUM),
            sv_panel=sv_panel,
            mutation_panel=profile_panels["mutation"].get(sample_id),
            cna_panel=profile_panels["cna"].get(sample_id),
        )
        for sample_id, sv_panel in sorted(sv_panels.items())
    ]
    return DiscoveryInputs(
        samples=samples,
        panel_genes=panel_genes,
        alterations=alterations,
        mutation_count_by_sample=mutation_counts,
        warnings=warnings,
    )
