"""Batch auto-generated :class:`~cfh.genes.registry.GeneConfig` construction
for recurrence-gated genes with no hand-curated YAML config.

Genome Nexus's ``/ensembl/canonical-transcript/hgnc`` endpoint accepts a
*list* of gene symbols and returns canonical transcript + protein length +
Pfam domain coordinates for all of them in one batch call, so an entire
genome-wide candidate set can be resolved in as few network round-trips as
possible (chunked only when the candidate count exceeds ``batch_size``).

A curated :mod:`cfh.genes.registry` YAML config (BRAF, RET, ...) always
takes precedence -- this module is only ever consulted for a gene that has
none.
"""

from __future__ import annotations

import json
import time
from dataclasses import dataclass
from pathlib import Path

import requests

from cfh.genes.registry import GeneConfig, KeyDomain, derive_gene_config_defaults
from cfh.mapping.genome_nexus_source import (
    CanonicalTranscript,
    PfamDomain,
    parse_canonical_transcript,
)

DEFAULT_GENOME_NEXUS_BASE_URL = "https://www.genomenexus.org"
DEFAULT_BATCH_SIZE = 800
"""Chunk size for the batch canonical-transcript call. Verified live with
836 genes in a single request; chunking only guards against a candidate set
larger than that, not a documented hard limit."""

DEFAULT_PFAM_DESCRIPTION_BASE_URL = "https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam"
KINASE_DOMAIN_KEYWORDS = ("kinase", "catalytic")

_RETRYABLE_STATUS_CODES = {429, 500, 502, 503, 504}


class PfamDescriptionSource:
    """Look up a Pfam family's human-readable name/description by accession.

    Genome Nexus's canonical-transcript payload gives only the bare Pfam
    accession (e.g. ``"PF07714"``), never a description, so the
    kinase/catalytic key-domain heuristic needs this separate lookup
    against InterPro's Pfam entry API. Results are cached in-memory and, if
    ``cache_dir`` is supplied, on disk (one small JSON file per accession)
    -- the same on-disk caching pattern already used by
    :class:`~cfh.mapping.domain_source.UniProtDomainSource`.

    Any failure (network error, non-200 response, unexpected payload shape)
    degrades to ``None`` rather than raising: a domain description is an
    optional enrichment for the auto-config heuristic, which always has a
    largest-domain fallback available.
    """

    def __init__(
        self,
        session: "requests.Session | None" = None,
        cache_dir: str | Path | None = None,
        base_url: str = DEFAULT_PFAM_DESCRIPTION_BASE_URL,
        max_retries: int = 2,
        backoff_seconds: float = 0.3,
    ):
        self.session = session or requests.Session()
        self.cache_dir = Path(cache_dir) if cache_dir else None
        self.base_url = base_url.rstrip("/")
        self.max_retries = max_retries
        self.backoff_seconds = backoff_seconds
        self._cache: dict[str, str | None] = {}

    def _cache_file(self, accession: str) -> Path | None:
        return self.cache_dir / f"{accession}.json" if self.cache_dir else None

    def describe(self, accession: str) -> str | None:
        """Return a description for ``accession``, or ``None`` if unavailable."""
        if accession in self._cache:
            return self._cache[accession]

        cache_file = self._cache_file(accession)
        if cache_file and cache_file.exists():
            try:
                cached = json.loads(cache_file.read_text())
                self._cache[accession] = cached.get("description")
                return self._cache[accession]
            except (OSError, json.JSONDecodeError):
                pass

        description = self._fetch(accession)
        self._cache[accession] = description
        if cache_file:
            cache_file.parent.mkdir(parents=True, exist_ok=True)
            cache_file.write_text(json.dumps({"accession": accession, "description": description}))
        return description

    def _fetch(self, accession: str) -> str | None:
        url = f"{self.base_url}/{accession}"
        attempt = 0
        while True:
            try:
                response = self.session.get(url, timeout=15)
            except requests.RequestException:
                return None
            if response.status_code == 200:
                break
            if response.status_code not in _RETRYABLE_STATUS_CODES or attempt >= self.max_retries:
                return None
            time.sleep(self.backoff_seconds * (2**attempt))
            attempt += 1
        try:
            name = response.json()["metadata"]["name"]
        except (ValueError, KeyError, TypeError):
            return None
        if isinstance(name, dict):
            return name.get("name")
        return str(name) if name is not None else None


def select_key_domain(
    pfam_domains: list[PfamDomain],
    description_source: "PfamDescriptionSource | None" = None,
) -> KeyDomain | None:
    """Pick one Pfam domain to use as a gene's auto-generated key domain.

    Prefers a domain whose Pfam/InterPro description contains "kinase" or
    "catalytic" (case-insensitive); among any such matches, and as the sole
    fallback when none match (or no ``description_source`` is supplied),
    picks the largest annotated domain by amino-acid span. Returns ``None``
    when there are no Pfam domains at all to choose from.
    """
    if not pfam_domains:
        return None

    described = [
        (domain, (description_source.describe(domain.pfam_id) or "") if description_source else "")
        for domain in pfam_domains
    ]
    kinase_matches = [
        domain
        for domain, description in described
        if any(keyword in description.lower() for keyword in KINASE_DOMAIN_KEYWORDS)
    ]
    candidates = kinase_matches or [domain for domain, _ in described]
    chosen = max(candidates, key=lambda domain: (domain.end_aa - domain.start_aa, domain.pfam_id))
    chosen_description = next(
        (description for domain, description in described if domain is chosen), ""
    )
    return KeyDomain(
        name=chosen_description or chosen.pfam_id,
        source="genome_nexus",
        key="kinase" if chosen in kinase_matches else "auto_key_domain",
        accession=chosen.pfam_id,
    )


def build_auto_gene_config(
    gene_symbol: str,
    entrez_gene_id: int | None,
    canonical: CanonicalTranscript,
    *,
    description_source: "PfamDescriptionSource | None" = None,
) -> GeneConfig | None:
    """Build a minimal auto-generated ``GeneConfig`` from a resolved canonical
    transcript: canonical transcript id, protein id, and (when available) a
    single best-guess ``key_domains`` entry.

    Coordinate-derived defaults for ``disruption_required_domains`` and
    ``expected_retained_exon_hint`` are populated from the same complete
    canonical-transcript response. ``gene_pair`` remains unset.

    Returns ``None`` when the canonical transcript has no protein id (a
    ``GeneConfig`` requires one), which callers must treat as "this gene
    cannot be auto-configured" rather than crash on.
    """
    if not canonical.protein_id or not canonical.transcript_id:
        return None
    key_domain = select_key_domain(canonical.pfam_domains, description_source)
    config = GeneConfig(
        gene_symbol=gene_symbol.upper(),
        canonical_transcript_id=canonical.refseq_mrna_id or canonical.transcript_id,
        protein_id=canonical.protein_id,
        key_domains=[key_domain] if key_domain else [],
        entrez_gene_id=entrez_gene_id,
    )
    return derive_gene_config_defaults(
        config,
        canonical,
        domain_name_resolver=description_source.describe if description_source else None,
    )


@dataclass
class BatchCanonicalTranscriptResult:
    by_gene_symbol: dict[str, CanonicalTranscript]
    unresolved_gene_symbols: list[str]
    """Requested symbols Genome Nexus had no canonical-transcript mapping
    for at all -- distinct from a batch-call network failure, which raises
    instead of silently omitting genes."""


def batch_fetch_canonical_transcripts(
    gene_symbols: list[str],
    *,
    base_url: str = DEFAULT_GENOME_NEXUS_BASE_URL,
    session: "requests.Session | None" = None,
    cache_dir: str | Path | None = None,
    batch_size: int = DEFAULT_BATCH_SIZE,
    max_retries: int = 3,
    backoff_seconds: float = 0.5,
    timeout: float = 60,
    isoform_override_source: str | None = None,
) -> BatchCanonicalTranscriptResult:
    """Resolve canonical transcripts for many gene symbols in as few Genome
    Nexus batch calls as possible.

    Each gene's raw canonical-transcript payload is cached on disk (one
    JSON file per symbol under ``cache_dir``, when supplied), so a repeated
    cohort scan does not re-fetch genes it already resolved. Only genes not
    already cached are included in the batch POST body(ies); the body list
    is chunked at ``batch_size`` only when more candidates are requested at
    once than that.

    ``isoform_override_source`` asks Genome Nexus for another isoform set
    (e.g. ``"mskcc"``). Left ``None``, the request is unchanged; callers
    that set it must use a separate ``cache_dir``.
    """
    session = session or requests.Session()
    cache_paths: dict[str, Path] = {}
    result: dict[str, CanonicalTranscript] = {}
    to_fetch: list[str] = []

    for symbol in gene_symbols:
        cache_path = Path(cache_dir) / f"{symbol.upper()}.json" if cache_dir else None
        if cache_path is not None:
            cache_paths[symbol] = cache_path
        if cache_path is not None and cache_path.exists():
            try:
                payload = json.loads(cache_path.read_text())
                result[symbol] = parse_canonical_transcript(payload)
                continue
            except (OSError, json.JSONDecodeError, KeyError):
                pass
        to_fetch.append(symbol)

    url = f"{base_url.rstrip('/')}/ensembl/canonical-transcript/hgnc"
    for start in range(0, len(to_fetch), batch_size):
        chunk = to_fetch[start : start + batch_size]
        if not chunk:
            continue
        attempt = 0
        while True:
            response = session.post(
                url,
                json=chunk,
                params=(
                    {"isoformOverrideSource": isoform_override_source}
                    if isoform_override_source
                    else None
                ),
                timeout=timeout,
            )
            if response.status_code not in _RETRYABLE_STATUS_CODES or attempt >= max_retries:
                break
            time.sleep(backoff_seconds * (2**attempt))
            attempt += 1
        response.raise_for_status()
        payloads = response.json()

        payload_by_symbol: dict[str, dict] = {}
        for payload in payloads:
            for symbol in payload.get("hugoSymbols") or []:
                payload_by_symbol[symbol.upper()] = payload

        for symbol in chunk:
            payload = payload_by_symbol.get(symbol.upper())
            if payload is None:
                continue
            result[symbol] = parse_canonical_transcript(payload)
            cache_path = cache_paths.get(symbol)
            if cache_path:
                cache_path.parent.mkdir(parents=True, exist_ok=True)
                cache_path.write_text(json.dumps(payload))

    unresolved = [symbol for symbol in gene_symbols if symbol not in result]
    return BatchCanonicalTranscriptResult(by_gene_symbol=result, unresolved_gene_symbols=unresolved)
