#!/usr/bin/env python3
"""Curation-assist tool: ask OpenEvidence for a fusion's molecular mechanism.

This is a research aid for a human curator deciding what to write into a
gene's ``mechanism_note`` (see ``cfh.genes.registry.GeneConfig``) -- it
never writes to any file itself. It prints OpenEvidence's answer plus its
citation list so the curator can read, verify against the cited primary
literature, and manually transcribe (or discard) into a YAML config.

OpenEvidence's output here is AI-synthesized and NOT independently
verified the way this project's 11 hand-curated mechanism_note entries
are (each of those was cross-checked by a human against the actual cited
papers). If a curator chooses to use this output with only light editing,
set ``mechanism_note_verified: false`` on that gene's YAML config -- the
report and static HTML viewer render an unverified note under a visibly
distinct "AI-suggested mechanism (unverified)" label, never the same
"Curated mechanism" treatment a cross-checked note gets.

Question phrasing deliberately mirrors this codebase's own prior findings
using OpenEvidence (a closed, targeted question, not an open-ended
"summarize everything" ask): a live benchmark on a sibling project
(agentic-cancer-gene-classification) found open-ended prompts attached
ungrounded specific statistics to cited PMIDs.

Requires OPENEVIDENCE_API_KEY in the environment (or a .env file loaded
before running this script). OpenEvidence access is not part of this
project's own dependencies or committed pipeline -- this script is an
optional, manually-run research aid only.

Usage:
    OPENEVIDENCE_API_KEY=... python scripts/curate_mechanism_note.py EGFR::SEPTIN14
    OPENEVIDENCE_API_KEY=... python scripts/curate_mechanism_note.py FGFR2::BICC1 \
        --tumor-type "cholangiocarcinoma"
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from dataclasses import dataclass, field

import requests

STREAMING_ANALYSIS_PATH = "/streaming/analysis"
DEFAULT_BASE_URL = "https://api.openevidence.com"
DEFAULT_MODEL = "darwin"
DEFAULT_TIMEOUT_SECONDS = 280.0
"""OpenEvidence has no fixed SLA; a live-verified smoke test against the
real API for a moderately complex question took ~220s. This script is a
manually-run, one-off tool, so a generous timeout is preferable to a
premature failure -- there is no shared budget to protect, unlike a
production pipeline serving many concurrent requests."""


class OpenEvidenceConfigurationError(RuntimeError):
    """Raised when OPENEVIDENCE_API_KEY is not set."""


@dataclass
class Citation:
    citation_key: str
    title: str = ""
    authors: str = ""
    journal: str = ""
    date: str = ""
    doi: str = ""
    url: str = ""


@dataclass
class Analysis:
    question: str
    text: str = ""
    citations: list[Citation] = field(default_factory=list)


def build_question(gene1: str, gene2: str, tumor_type: str | None) -> str:
    """A closed, pointed question asking specifically for the mechanism,
    not a free-form summary or a classification restatement."""
    cancer_context = tumor_type if tumor_type else "cancer"
    return (
        f"Based on peer-reviewed evidence, what is the molecular mechanism by which "
        f"the {gene1}::{gene2} fusion is oncogenic in {cancer_context}? State "
        "specifically whether it works by retaining a kinase/effector domain while "
        "losing an autoinhibitory or ligand-binding domain, by partner-driven "
        "dimerization/oligomerization, by promoter-swap/overexpression, by creating "
        "a chimeric (neomorphic) transcription factor, or by another named mechanism, "
        "and cite the strongest supporting evidence."
    )


def _iter_sse_payloads(raw: str) -> list[str]:
    """Split an SSE stream body into raw per-event data payloads. Per the
    SSE spec, multiple `data:` lines within one event are joined with
    "\\n" between them, not concatenated directly."""
    payloads: list[str] = []
    for block in raw.replace("\r\n", "\n").split("\n\n"):
        data_lines = [
            line[len("data:") :].strip() for line in block.splitlines() if line.startswith("data:")
        ]
        if data_lines:
            payloads.append("\n".join(data_lines))
    return payloads


def _parse_sse_events(raw: str) -> list[dict]:
    """Parse an SSE stream body into JSON event payloads. Malformed
    payloads are skipped rather than failing the whole parse."""
    events: list[dict] = []
    for payload in _iter_sse_payloads(raw):
        try:
            event = json.loads(payload)
        except json.JSONDecodeError:
            continue
        if isinstance(event, dict):
            events.append(event)
    return events


def _citation_from_event(event: dict) -> Citation | None:
    """Citation data is nested under event["reference"], with
    bibliographic fields a further level deep under
    event["reference"]["reference_detail"]."""
    reference = event.get("reference")
    if not isinstance(reference, dict):
        return None
    citation_key = reference.get("citation_key")
    if citation_key is None:
        return None
    detail = reference.get("reference_detail")
    if not isinstance(detail, dict):
        detail = {}
    return Citation(
        citation_key=str(citation_key),
        title=detail.get("title") or "",
        authors=detail.get("authors_string") or "",
        journal=detail.get("journal_name") or detail.get("journal_short_name") or "",
        date=detail.get("publication_date") or "",
        doi=detail.get("doi") or "",
        url=detail.get("url") or "",
    )


_PROGRESS_MARKER = "REACTCOMPONENT!:!InlineGenerationStep!:!"


def _strip_progress_marker(text: str) -> str:
    """OpenEvidence's first text delta is a UI progress-step JSON blob
    (e.g. "Analyzing query" / "Searching published medical literature..."),
    not part of the actual analysis prose -- strip it if present. Uses
    ``json.JSONDecoder.raw_decode`` (not a brace-counting regex) to skip
    exactly the one JSON object correctly even though it contains nested
    braces."""
    if not text.startswith(_PROGRESS_MARKER):
        return text
    remainder = text[len(_PROGRESS_MARKER) :]
    try:
        _, end_index = json.JSONDecoder().raw_decode(remainder)
    except json.JSONDecodeError:
        return text
    return remainder[end_index:].lstrip("\n")


def build_analysis(question: str, events: list[dict]) -> Analysis:
    """Accumulate prose text and dedupe citations across all events. Per
    OpenEvidence's own docs, concatenating every event's `text` field
    produces the full analysis text, citation-bearing events included
    (their `text` is typically an inline marker like "[[1]]"). `table`
    events are silently skipped -- no rendering for tabular content here."""
    text_parts: list[str] = []
    citations_by_key: dict[str, Citation] = {}
    for event in events:
        if "table" in event:
            continue
        delta = event.get("text")
        if delta:
            text_parts.append(delta)
        citation = _citation_from_event(event)
        if citation is not None:
            citations_by_key.setdefault(citation.citation_key, citation)
    text = _strip_progress_marker("".join(text_parts))
    return Analysis(question=question, text=text, citations=list(citations_by_key.values()))


def fetch_mechanism_analysis(
    gene1: str,
    gene2: str,
    *,
    tumor_type: str | None = None,
    api_key: str | None = None,
    base_url: str = DEFAULT_BASE_URL,
    model: str = DEFAULT_MODEL,
    timeout_seconds: float = DEFAULT_TIMEOUT_SECONDS,
) -> Analysis:
    api_key = api_key or os.environ.get("OPENEVIDENCE_API_KEY")
    if not api_key:
        raise OpenEvidenceConfigurationError(
            "OPENEVIDENCE_API_KEY is required (set it in the environment or a .env file)."
        )
    question = build_question(gene1, gene2, tumor_type)
    url = f"{base_url.rstrip('/')}{STREAMING_ANALYSIS_PATH}"
    headers = {"Authorization": f"Token {api_key}", "Accept": "text/event-stream"}
    payload = {"text": question, "model": model}
    response = requests.post(
        url, json=payload, headers=headers, timeout=timeout_seconds, stream=True
    )
    response.raise_for_status()
    raw = response.text
    events = _parse_sse_events(raw)
    return build_analysis(question, events)


def _print_analysis(analysis: Analysis) -> None:
    print("=" * 80)
    print("UNVERIFIED -- AI-synthesized via OpenEvidence, not independently verified.")
    print("Cross-check every claim below against the cited papers before using it in")
    print("a gene config; if used with only light editing, set mechanism_note_verified:")
    print("false so the report/viewer never present it as a fully curated note.")
    print("=" * 80)
    print()
    print("Question:", analysis.question)
    print()
    print(analysis.text)
    print()
    if analysis.citations:
        print("-" * 80)
        print(f"Citations ({len(analysis.citations)}):")
        for citation in sorted(analysis.citations, key=lambda c: int(c.citation_key)):
            print(f"[{citation.citation_key}] {citation.authors} {citation.title}")
            print(f"    {citation.journal} {citation.date}".rstrip())
            if citation.doi:
                print(f"    https://doi.org/{citation.doi}")
            elif citation.url:
                print(f"    {citation.url}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("fusion", help='Fusion in GENE1::GENE2 form, e.g. "EGFR::SEPTIN14"')
    parser.add_argument(
        "--tumor-type",
        default=None,
        help='Tumor type context, e.g. "cholangiocarcinoma" (default: generic "cancer")',
    )
    parser.add_argument("--base-url", default=DEFAULT_BASE_URL)
    parser.add_argument("--model", default=DEFAULT_MODEL)
    parser.add_argument("--timeout-seconds", type=float, default=DEFAULT_TIMEOUT_SECONDS)
    args = parser.parse_args(argv)

    if "::" not in args.fusion:
        parser.error('fusion must be in GENE1::GENE2 form, e.g. "EGFR::SEPTIN14"')
    gene1, gene2 = args.fusion.split("::", 1)

    try:
        analysis = fetch_mechanism_analysis(
            gene1,
            gene2,
            tumor_type=args.tumor_type,
            base_url=args.base_url,
            model=args.model,
            timeout_seconds=args.timeout_seconds,
        )
    except OpenEvidenceConfigurationError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    except requests.RequestException as exc:
        print(f"OpenEvidence request failed: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 1

    _print_analysis(analysis)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
