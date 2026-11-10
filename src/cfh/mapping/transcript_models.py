"""Build :class:`~cfh.mapping.junction_translation.TranscriptModel` objects.

Exon structure and Pfam domains come from Genome Nexus, through the same
client and on-disk cache the cohort scan uses. Coding sequence comes from
Ensembl REST, because Genome Nexus does not serve it. Both are cached per
transcript, so a repeated screen does not re-fetch.

Each partner uses the transcript the source call annotated when one exists,
and is reported with which transcript was used and whether it is the
canonical one.
"""

from __future__ import annotations

import json
import re
import time
from dataclasses import dataclass
from pathlib import Path

import requests

from cfh.cohort.auto_config import batch_fetch_canonical_transcripts
from cfh.mapping.genome_nexus_source import (
    DEFAULT_ISOFORM_OVERRIDE_SOURCE,
    CanonicalTranscript,
    GenomeNexusClient,
    GenomeNexusGeneNotFound,
    cds_bounds_from_utrs,
    parse_canonical_transcript,
)
from cfh.mapping.junction_translation import (
    TranscriptDomain,
    TranscriptExon,
    TranscriptModel,
    normalize_genome_build,
)

SOURCE_ANNOTATED = "source_annotated"
CANONICAL_NO_ANNOTATION = "canonical_no_source_annotation"
CANONICAL_ANNOTATION_UNMATCHED = "canonical_source_annotation_not_found"

_ANNOTATION_TRANSCRIPT = re.compile(
    r"([A-Za-z0-9._-]+)\s*\(\s*((?:NM_\d+(?:\.\d+)?|ENST\d+(?:\.\d+)?)"
    r"(?:\s*,\s*(?:NM_\d+(?:\.\d+)?|ENST\d+(?:\.\d+)?))*)\s*\)"
)

_ENSEMBL_HOST_BY_BUILD = {
    "GRCh37": "https://grch37.rest.ensembl.org",
    "GRCh38": "https://rest.ensembl.org",
}
_RETRYABLE_STATUS_CODES = {429, 500, 502, 503, 504}


def genome_build_of_genome_nexus(base_url: str) -> str:
    """The assembly a Genome Nexus instance serves: its GRCh38 host, else GRCh37."""
    return "GRCh38" if "grch38" in base_url.lower() else "GRCh37"


def source_annotated_transcripts(annotation: str | None, gene: str) -> list[str]:
    """Transcript ids the free-text annotation names for ``gene``, in order.

    Matches ``GENE (NM_x)`` and ``GENE ( NM_x, NM_y)``, and isoform-qualified
    symbols such as ``CDKN2Ap14ARF (NM_058195)``. Versions are removed.
    """
    if not annotation:
        return []
    symbol = re.compile(rf"^{re.escape(gene)}(?:p\d+\w*)?$", re.IGNORECASE)
    found: list[str] = []
    for match in _ANNOTATION_TRANSCRIPT.finditer(annotation):
        if not symbol.match(match.group(1)):
            continue
        for transcript in match.group(2).split(","):
            transcript = transcript.strip().split(".")[0]
            if transcript not in found:
                found.append(transcript)
    return found


@dataclass(frozen=True)
class ResolvedTranscript:
    """Which transcript a partner uses, before any sequence is fetched."""

    gene: str
    transcript: CanonicalTranscript | None
    source: str | None
    canonical_transcript_id: str | None
    differs_from_canonical: bool | None
    annotated_transcripts: tuple[str, ...] = ()
    note: str | None = None

    @property
    def transcript_id(self) -> str | None:
        return self.transcript.transcript_id if self.transcript else None

    @property
    def refseq_id(self) -> str | None:
        return self.transcript.refseq_mrna_id if self.transcript else None

    @property
    def families(self) -> set[str]:
        return {d.pfam_id for d in self.transcript.pfam_domains} if self.transcript else set()

    @property
    def span(self) -> tuple[int, int] | None:
        if not self.transcript or not self.transcript.exons:
            return None
        return (
            min(e.start for e in self.transcript.exons),
            max(e.end for e in self.transcript.exons),
        )


class EnsemblCdsClient:
    """Coding sequence by Ensembl transcript id, for one genome build."""

    def __init__(
        self,
        genome_build: str,
        *,
        session: requests.Session | None = None,
        cache_dir: str | Path | None = None,
        max_retries: int = 3,
        backoff_seconds: float = 0.5,
    ):
        self.genome_build = genome_build
        self.session = session or requests.Session()
        self.cache_dir = Path(cache_dir) if cache_dir else None
        self.max_retries = max_retries
        self.backoff_seconds = backoff_seconds

    def _cache_file(self, transcript_id: str) -> Path | None:
        if self.cache_dir is None:
            return None
        return self.cache_dir / f"{transcript_id}.json"

    def fetch_cds(self, transcript_id: str) -> str | None:
        """The transcript's CDS, or ``None`` when Ensembl has none for it."""
        cache_file = self._cache_file(transcript_id)
        if cache_file is not None and cache_file.exists():
            try:
                return json.loads(cache_file.read_text()).get("seq")
            except (OSError, json.JSONDecodeError):
                pass
        host = _ENSEMBL_HOST_BY_BUILD.get(self.genome_build)
        if host is None:
            return None
        attempt = 0
        while True:
            response = self.session.get(
                f"{host}/sequence/id/{transcript_id}",
                params={"type": "cds"},
                headers={"Content-Type": "application/json"},
                timeout=30,
            )
            if response.status_code not in _RETRYABLE_STATUS_CODES or attempt >= self.max_retries:
                break
            time.sleep(self.backoff_seconds * (2**attempt))
            attempt += 1
        if response.status_code in (400, 404):
            return None
        response.raise_for_status()
        sequence = response.json().get("seq")
        if cache_file is not None and sequence:
            cache_file.parent.mkdir(parents=True, exist_ok=True)
            cache_file.write_text(json.dumps({"id": transcript_id, "seq": sequence}))
        return sequence


class TranscriptModelSource:
    """Resolves each partner's transcript and builds its model on demand."""

    def __init__(
        self,
        genome_build: str,
        *,
        genome_nexus_client: GenomeNexusClient | None = None,
        cds_client: EnsemblCdsClient | None = None,
        cache_dir: str | Path | None = None,
        session: requests.Session | None = None,
    ):
        build = normalize_genome_build(genome_build)
        if build is None:
            raise ValueError("TranscriptModelSource requires an explicit genome build")
        self.genome_build = build
        self.client = genome_nexus_client or GenomeNexusClient(session=session)
        self.cache_dir = Path(cache_dir) if cache_dir else None
        self.cds_client = cds_client or EnsemblCdsClient(
            build,
            session=session,
            cache_dir=self.cache_dir / "ensembl_cds" / build if self.cache_dir else None,
        )
        self._canonical: dict[str, CanonicalTranscript | None] = {}
        self._isoforms: dict[str, list[dict]] = {}
        self._models: dict[str, tuple[TranscriptModel | None, str | None]] = {}

    @property
    def _canonical_cache_dir(self) -> Path | None:
        return (
            self.cache_dir / "genome_nexus_canonical_transcripts_mskcc" if self.cache_dir else None
        )

    def prefetch_canonical(self, genes: list[str]) -> None:
        """Resolve many genes through the cohort scan's batch and cache path."""
        pending = sorted({gene for gene in genes if gene not in self._canonical})
        if not pending:
            return
        result = batch_fetch_canonical_transcripts(
            pending,
            base_url=self.client.base_url,
            session=self.client.session,
            cache_dir=self._canonical_cache_dir,
            isoform_override_source=DEFAULT_ISOFORM_OVERRIDE_SOURCE,
        )
        for gene, transcript in result.by_gene_symbol.items():
            self._canonical[gene] = transcript

    def canonical(self, gene: str) -> CanonicalTranscript | None:
        if gene in self._canonical:
            return self._canonical[gene]
        cache_dir = self._canonical_cache_dir
        cache_file = cache_dir / f"{gene.upper()}.json" if cache_dir else None
        if cache_file is not None and cache_file.exists():
            try:
                self._canonical[gene] = parse_canonical_transcript(
                    json.loads(cache_file.read_text())
                )
                return self._canonical[gene]
            except (OSError, json.JSONDecodeError, KeyError):
                pass
        try:
            payload = self.client.fetch_canonical_transcript(gene)
        except GenomeNexusGeneNotFound:
            self._canonical[gene] = None
            return None
        if cache_file is not None:
            cache_file.parent.mkdir(parents=True, exist_ok=True)
            cache_file.write_text(json.dumps(payload))
        self._canonical[gene] = parse_canonical_transcript(payload)
        return self._canonical[gene]

    def _gene_isoforms(self, gene: str) -> list[dict]:
        if gene in self._isoforms:
            return self._isoforms[gene]
        cache_file = (
            self.cache_dir / "genome_nexus_transcripts" / f"{gene.upper()}.json"
            if self.cache_dir
            else None
        )
        payloads: list[dict] | None = None
        if cache_file is not None and cache_file.exists():
            try:
                payloads = json.loads(cache_file.read_text())
            except (OSError, json.JSONDecodeError):
                payloads = None
        if payloads is None:
            payloads = self.client.fetch_transcripts_for_gene(gene)
            if cache_file is not None:
                cache_file.parent.mkdir(parents=True, exist_ok=True)
                cache_file.write_text(json.dumps(payloads))
        self._isoforms[gene] = payloads
        return payloads

    def resolve(self, gene: str, annotated: list[str] | None = None) -> ResolvedTranscript:
        """Use the first annotated transcript that exists, else the canonical."""
        canonical = self.canonical(gene)
        canonical_id = canonical.transcript_id if canonical else None
        annotated = annotated or []
        if annotated:
            matched = self._match_annotated(gene, annotated, canonical_id)
            if matched is not None:
                transcript, note = matched
                return ResolvedTranscript(
                    gene,
                    transcript,
                    SOURCE_ANNOTATED,
                    canonical_id,
                    transcript.transcript_id != canonical_id if canonical_id else None,
                    tuple(annotated),
                    note,
                )
        if canonical is None:
            return ResolvedTranscript(
                gene,
                None,
                None,
                None,
                None,
                tuple(annotated),
                f"Genome Nexus has no transcript for {gene}",
            )
        return ResolvedTranscript(
            gene,
            canonical,
            CANONICAL_ANNOTATION_UNMATCHED if annotated else CANONICAL_NO_ANNOTATION,
            canonical_id,
            False,
            tuple(annotated),
            (
                f"annotated transcript(s) {', '.join(annotated)} not found for {gene}; "
                "canonical used"
                if annotated
                else None
            ),
        )

    def _match_annotated(
        self, gene: str, annotated: list[str], canonical_id: str | None
    ) -> tuple[CanonicalTranscript, str | None] | None:
        """First annotated id with an isoform; among several, a repeatable pick."""
        for position, wanted in enumerate(annotated):
            hits = [
                payload
                for payload in self._gene_isoforms(gene)
                if wanted
                in (
                    (payload.get("refseqMrnaId") or "").split(".")[0],
                    (payload.get("transcriptId") or "").split(".")[0],
                )
                and payload.get("proteinLength")
                and payload.get("exons")
            ]
            if not hits:
                continue
            best = min(
                hits,
                key=lambda payload: (
                    payload.get("transcriptId") != canonical_id,
                    payload.get("ccdsId") is None,
                    -(payload.get("proteinLength") or 0),
                    payload.get("transcriptId") or "",
                ),
            )
            notes = []
            if len(hits) > 1:
                ids = ", ".join(sorted(payload["transcriptId"] for payload in hits))
                notes.append(f"{wanted} maps to {ids}; used {best['transcriptId']}")
            if len(annotated) > 1:
                notes.append(f"annotation names {', '.join(annotated)}; used {wanted}")
            if position:
                notes.append(f"{', '.join(annotated[:position])} had no isoform")
            transcript = parse_canonical_transcript(
                {
                    **best,
                    "pfamDomains": best.get("pfamDomains") or [],
                    "exons": best.get("exons") or [],
                    "utrs": best.get("utrs") or [],
                }
            )
            return transcript, "; ".join(notes) or None
        return None

    def model(self, resolved: ResolvedTranscript) -> tuple[TranscriptModel | None, str | None]:
        """``(model, problem)``; exactly one is ``None``."""
        transcript = resolved.transcript
        if transcript is None:
            return None, resolved.note or f"no transcript for {resolved.gene}"
        cached = self._models.get(transcript.transcript_id)
        if cached is not None:
            return cached
        built = self._build(resolved.gene, transcript)
        self._models[transcript.transcript_id] = built
        return built

    def _build(
        self, gene: str, transcript: CanonicalTranscript
    ) -> tuple[TranscriptModel | None, str | None]:
        if not transcript.exons:
            return None, f"{transcript.transcript_id} has no exon coordinates"
        if not transcript.protein_length:
            return None, f"{transcript.transcript_id} is not protein coding"
        sequence = self.cds_client.fetch_cds(transcript.transcript_id)
        if not sequence:
            return None, f"no CDS sequence for {transcript.transcript_id} ({self.genome_build})"
        cds_min, cds_max = cds_bounds_from_utrs(transcript.utrs)
        exons = tuple(TranscriptExon(e.rank, e.start, e.end) for e in transcript.exons)
        model = TranscriptModel(
            gene=gene,
            transcript_id=transcript.transcript_id,
            genome_build=self.genome_build,
            strand=transcript.exons[0].strand,
            exons=exons,
            cds_start=cds_min if cds_min is not None else min(e.start for e in exons),
            cds_end=cds_max if cds_max is not None else max(e.end for e in exons),
            cds_sequence=sequence,
            domains=tuple(
                TranscriptDomain(d.pfam_id, d.start_aa, d.end_aa) for d in transcript.pfam_domains
            ),
        )
        if len(model.protein) != transcript.protein_length:
            return None, (
                f"{transcript.transcript_id} translates to {len(model.protein)} aa, "
                f"Genome Nexus reports {transcript.protein_length}"
            )
        return model, None
