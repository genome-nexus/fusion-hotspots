"""Parent key-residue annotations from UniProtKB/Swiss-Prot (REST API, cached per gene).

Each gene symbol maps to the reviewed human entry whose primary gene name is
that symbol. Key features are UniProt's own annotations of these types:

* point-like: ``Active site``, ``Binding site`` (includes metal binding),
  ``Site`` (except translocation-breakpoint sites);
* spans: ``DNA binding``, ``Motif``, ``Region`` (except ``Disordered``).

A key feature whose start or end is not an exact integer position (UniProt
modifiers ``UNKNOWN`` / ``UNSURE`` / ``OUTSIDE``, or a missing value) cannot be
mapped. It is kept as an :class:`IncompleteFeature` rather than dropped, so a
domain it may overlap is reported as ``incomplete_annotation`` (not assessed),
never as "assessed, nothing annotated".
"""

from __future__ import annotations

import json
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import requests

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
FEATURE_FIELDS = "ft_act_site,ft_binding,ft_dna_bind,ft_site,ft_motif,ft_region"
KEY_FEATURE_TYPES = ("Active site", "Binding site", "Site", "DNA binding", "Motif", "Region")
POINT_FEATURE_TYPES = ("Active site", "Binding site", "Site")


@dataclass(frozen=True)
class KeyFeature:
    type: str
    start: int
    end: int
    description: str

    @property
    def label(self) -> str:
        text = f"{self.type}:{self.start}" + (f"-{self.end}" if self.end != self.start else "")
        return f"{text}({self.description})" if self.description else text


@dataclass(frozen=True)
class IncompleteFeature:
    """A key feature with a non-exact start and/or end (``None`` = position not usable)."""

    type: str
    start: int | None
    end: int | None
    description: str
    raw_location: str

    @property
    def label(self) -> str:
        return f"{self.type}:{self.raw_location}" + (
            f"({self.description})" if self.description else ""
        )

    def may_overlap(self, low: int, high: int) -> bool:
        """Whether the feature could touch UniProt positions ``low..high``."""
        return (self.start is None or self.start <= high) and (self.end is None or self.end >= low)


@dataclass(frozen=True)
class UniProtEntry:
    gene: str
    accession: str
    sequence: str
    features: tuple[KeyFeature, ...]
    incomplete: tuple[IncompleteFeature, ...] = ()
    """Key features whose coordinates are not exact (see module docstring)."""


def is_key_feature(raw: dict[str, Any]) -> bool:
    kind = raw.get("type")
    description = str(raw.get("description") or "")
    if kind not in KEY_FEATURE_TYPES:
        return False
    if kind == "Region" and description.lower().startswith("disordered"):
        return False
    if kind == "Site" and "breakpoint" in description.lower():
        return False
    return True


def _exact(bound: dict[str, Any]) -> int | None:
    value = bound.get("value")
    modifier = bound.get("modifier", "EXACT")
    if isinstance(value, bool) or not isinstance(value, int) or modifier != "EXACT":
        return None
    return value


def _describe(bound: dict[str, Any]) -> str:
    return f"{bound.get('modifier', '?')}:{bound.get('value')}"


def parse_entry(gene: str, raw: dict[str, Any]) -> UniProtEntry:
    """A :class:`UniProtEntry` from a UniProt JSON result.

    Key features with exact integer coordinates become :class:`KeyFeature`; the
    rest become :class:`IncompleteFeature` (never silently dropped).
    """
    features = []
    incomplete = []
    for item in raw.get("features", []):
        if not is_key_feature(item):
            continue
        location = item.get("location") or {}
        start_raw, end_raw = location.get("start") or {}, location.get("end") or {}
        start, end = _exact(start_raw), _exact(end_raw)
        description = str(item.get("description") or "")
        ligand = (item.get("ligand") or {}).get("name")
        if ligand:
            description = f"{description} {ligand}".strip()
        if start is None or end is None or end < start:
            # An unknown bound stays None, so may_overlap() errs towards overlapping;
            # an inverted exact pair is unusable on both sides.
            inverted = start is not None and end is not None and end < start
            incomplete.append(
                IncompleteFeature(
                    str(item["type"]),
                    None if inverted else start,
                    None if inverted else end,
                    description,
                    f"{_describe(start_raw)}..{_describe(end_raw)}",
                )
            )
            continue
        features.append(KeyFeature(str(item["type"]), start, end, description))
    return UniProtEntry(
        gene,
        raw["primaryAccession"],
        raw["sequence"]["value"],
        tuple(sorted(set(features), key=_order)),
        tuple(incomplete),
    )


def _order(feature: KeyFeature) -> tuple[int, int, str]:
    return feature.start, feature.end, feature.type


def _primary_gene(raw: dict[str, Any]) -> str | None:
    genes = raw.get("genes") or []
    if not genes or "geneName" not in genes[0]:
        return None
    return str(genes[0]["geneName"]["value"])


class UniProtClient:
    """Batched, cached gene -> reviewed human UniProt entry lookups."""

    def __init__(self, cache_dir: Path, *, session: requests.Session | None = None):
        self.cache_dir = cache_dir
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.session = session or requests.Session()
        self.release: str | None = None

    def _path(self, gene: str) -> Path:
        return self.cache_dir / f"{gene.upper()}.json"

    def prefetch(self, genes: Iterable[str], batch_size: int = 40) -> None:
        pending = sorted({g for g in genes if g and not self._path(g).exists()})
        for start in range(0, len(pending), batch_size):
            batch = pending[start : start + batch_size]
            query = "(" + " OR ".join(f'gene_exact:"{g}"' for g in batch) + ")"
            query += " AND organism_id:9606 AND reviewed:true"
            results = self._search(query)
            by_gene: dict[str, dict[str, Any]] = {}
            for raw in results:
                primary = _primary_gene(raw)
                if primary in batch and primary not in by_gene:
                    by_gene[primary] = raw
            for gene in batch:
                payload = {"gene": gene, "release": self.release, "entry": by_gene.get(gene)}
                self._path(gene).write_text(json.dumps(payload))

    def _search(self, query: str) -> list[dict[str, Any]]:
        params = {
            "query": query,
            "fields": f"accession,gene_primary,sequence,{FEATURE_FIELDS}",
            "format": "json",
            "size": "500",
        }
        last: Exception | None = None
        for attempt in range(4):
            try:
                response = self.session.get(UNIPROT_SEARCH, params=params, timeout=120)
                response.raise_for_status()
                self.release = response.headers.get("X-UniProt-Release", self.release)
                return list(response.json()["results"])
            except (requests.RequestException, ValueError, KeyError) as error:
                last = error
                time.sleep(3 * (attempt + 1))
        raise RuntimeError(f"UniProt search failed: {last}")

    def entry(self, gene: str) -> UniProtEntry | None:
        """The cached entry; None when UniProt has no reviewed human entry for ``gene``."""
        path = self._path(gene)
        if not path.exists():
            self.prefetch([gene])
        payload = json.loads(path.read_text())
        if payload.get("release"):
            self.release = self.release or payload["release"]
        raw = payload.get("entry")
        return parse_entry(gene, raw) if raw else None
