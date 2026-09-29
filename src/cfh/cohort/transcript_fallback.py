"""Identity-verified fallback for genes the batch canonical-transcript call misses.

Genome Nexus's batch ``/ensembl/canonical-transcript/hgnc`` response is keyed
by each payload's ``hugoSymbols``. For some genes that list is stale or names
an alias-colliding gene (observed live: TCF3 labelled ``TCF7L1``, ERF labelled
``ETF1``, SEPTIN14 labelled ``SEPT14``), so the batch path correctly refuses
to attach those payloads and the gene is reported unresolved.

This fallback asks HGNC (by Entrez ID when known, else by symbol) for the
gene's approved symbol, Ensembl gene ID, and locus group, then:

* non-protein-coding loci are classified ``non_coding`` -- domain analysis
  does not apply, which is different from a failed lookup;
* otherwise the single-gene Genome Nexus canonical transcript (the same
  endpoint the per-gene pipeline later uses) is accepted only when its
  ``geneId`` equals HGNC's Ensembl gene ID and it has a protein.

Ensembl gene IDs are stable across genome builds, so this check does not mix
GRCh37 coordinates (Genome Nexus's default instance, MSK-IMPACT) with the
GRCh38 annotation HGNC reports. No coordinates are taken from HGNC.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import requests

from cfh.mapping.genome_nexus_source import (
    CanonicalTranscript,
    GenomeNexusClient,
    GenomeNexusGeneNotFound,
    parse_canonical_transcript,
)

HGNC_BASE_URL = "https://rest.genenames.org"
PROTEIN_CODING_LOCUS_GROUP = "protein-coding gene"


@dataclass(frozen=True)
class GeneIdentity:
    symbol: str
    ensembl_gene_id: str | None
    locus_group: str | None


@dataclass(frozen=True)
class FallbackResolution:
    status: str  # "resolved" | "non_coding" | "unresolved"
    reason: str
    canonical: CanonicalTranscript | None = None
    identity: GeneIdentity | None = None


def fetch_hgnc_identity(
    gene_symbol: str,
    entrez_gene_id: int | None,
    *,
    session: requests.Session | None = None,
    base_url: str = HGNC_BASE_URL,
    cache_dir: str | Path | None = None,
    timeout: float = 30,
) -> GeneIdentity | None:
    """Return HGNC's record for a gene, or ``None`` when HGNC has none.

    Network failures raise ``requests.RequestException``; callers decide how
    to degrade. Successful lookups (including "not found") are cached on disk
    when ``cache_dir`` is supplied.
    """
    field, value = (
        ("entrez_id", str(entrez_gene_id))
        if entrez_gene_id is not None
        else ("symbol", gene_symbol.upper())
    )
    cache_file = Path(cache_dir) / f"{field}_{value}.json" if cache_dir else None
    payload: dict | None = None
    if cache_file is not None and cache_file.exists():
        try:
            payload = json.loads(cache_file.read_text())
        except (OSError, json.JSONDecodeError):
            payload = None
    if payload is None:
        session = session or requests.Session()
        response = session.get(
            f"{base_url.rstrip('/')}/fetch/{field}/{value}",
            headers={"Accept": "application/json"},
            timeout=timeout,
        )
        response.raise_for_status()
        payload = response.json()
        if cache_file is not None:
            cache_file.parent.mkdir(parents=True, exist_ok=True)
            cache_file.write_text(json.dumps(payload))
    docs = (payload.get("response") or {}).get("docs") if isinstance(payload, dict) else None
    if not docs:
        return None
    doc = docs[0]
    if not isinstance(doc, dict) or not doc.get("symbol"):
        return None
    return GeneIdentity(
        symbol=str(doc["symbol"]),
        ensembl_gene_id=doc.get("ensembl_gene_id"),
        locus_group=doc.get("locus_group"),
    )


def resolve_with_fallback(
    gene_symbol: str,
    entrez_gene_id: int | None,
    *,
    genome_nexus_client: GenomeNexusClient,
    session: requests.Session | None = None,
    hgnc_base_url: str = HGNC_BASE_URL,
    hgnc_cache_dir: str | Path | None = None,
) -> FallbackResolution:
    """Resolve one batch-missed gene, never raising for data-source problems."""
    try:
        identity = fetch_hgnc_identity(
            gene_symbol,
            entrez_gene_id,
            session=session,
            base_url=hgnc_base_url,
            cache_dir=hgnc_cache_dir,
        )
    except (requests.RequestException, ValueError) as exc:
        return FallbackResolution("unresolved", f"HGNC lookup failed: {type(exc).__name__}")
    if identity is None:
        return FallbackResolution("unresolved", "HGNC has no record for this gene")
    if identity.locus_group != PROTEIN_CODING_LOCUS_GROUP:
        return FallbackResolution(
            "non_coding",
            f"HGNC locus group is {identity.locus_group!r}; domain analysis does not apply",
            identity=identity,
        )
    if not identity.ensembl_gene_id:
        return FallbackResolution(
            "unresolved", "HGNC record has no Ensembl gene ID to verify against", identity=identity
        )
    try:
        payload = genome_nexus_client.fetch_canonical_transcript(gene_symbol)
    except GenomeNexusGeneNotFound:
        return FallbackResolution(
            "unresolved",
            "Genome Nexus has no canonical transcript for this symbol",
            identity=identity,
        )
    except (requests.RequestException, ValueError) as exc:
        return FallbackResolution(
            "unresolved", f"Genome Nexus lookup failed: {type(exc).__name__}", identity=identity
        )
    if payload.get("geneId") != identity.ensembl_gene_id:
        return FallbackResolution(
            "unresolved",
            f"Genome Nexus transcript belongs to {payload.get('geneId')!r}, "
            f"not HGNC's {identity.ensembl_gene_id!r}",
            identity=identity,
        )
    try:
        canonical = parse_canonical_transcript(payload)
    except (KeyError, TypeError) as exc:
        return FallbackResolution(
            "unresolved", f"Malformed Genome Nexus payload: {exc}", identity=identity
        )
    if not canonical.protein_id:
        return FallbackResolution(
            "unresolved",
            "Genome Nexus canonical transcript has no protein annotation",
            identity=identity,
        )
    return FallbackResolution(
        "resolved",
        "Genome Nexus single-gene transcript verified by HGNC Ensembl gene ID",
        canonical=canonical,
        identity=identity,
    )
