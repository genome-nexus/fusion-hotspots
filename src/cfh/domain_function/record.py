"""Per-chimera records: one per distinct chimeric domain sequence.

A record is rebuilt from the junction screen's own functions
(:func:`~cfh.mapping.junction_translation.predict_junction` and
:func:`~cfh.mapping.junction_domain_classifier.classify_family`) on the same
transcript models, so its chimeric domain is exactly the one the screen
classified. Screen rows (one per sample) that give the same chimeric domain
(same family, parent transcripts and flanked sequence) collapse into one record
that keeps every cohort, sample count, read support and artefact flag.
"""

from __future__ import annotations

import gzip
import hashlib
import json
from collections import Counter
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping

from cfh.domain_function import SCHEMA_VERSION
from cfh.mapping.junction_domain_classifier import (
    AMBIGUOUS_JUNCTION_REGISTER,
    CLEAN_RECONSTITUTION,
    GAP_OR_DELETION,
    TANDEM_DUPLICATION_OR_INSERTION,
    classify_family,
)
from cfh.mapping.junction_translation import predict_junction, translate
from cfh.mapping.transcript_models import TranscriptModelSource, source_annotated_transcripts

TRANSLATED_CLASSES = (
    CLEAN_RECONSTITUTION,
    TANDEM_DUPLICATION_OR_INSERTION,
    GAP_OR_DELETION,
    AMBIGUOUS_JUNCTION_REGISTER,
)
"""Screen classes with an in-frame chimeric domain to score."""

FLANK = 30
"""Residues of chimera/parent context kept on each side of the domain window."""

ARTEFACT_FLAG_COLUMNS = (
    "paralog_sequence_flag",
    "genomic_similarity_flag",
    "breakpoint_repeat_flag",
    "sample_burden_flag",
    "read_evidence_flag",
    "normal_recurrence_flag",
)
"""Pilot artefact flags (``pass`` / ``flagged`` / ``not_assessed``) carried per row."""


@dataclass
class ParentDomain:
    """One parent's side of the chimeric domain, in that parent's own residues."""

    gene: str
    transcript: str
    protein: str
    domain_start: int
    """1-based first residue of the parent's Pfam domain (screen annotation)."""
    domain_end: int
    junction_residue: int
    """5' parent: last residue kept. 3' parent: first complete residue kept."""

    @property
    def domain_sequence(self) -> str:
        return self.protein[self.domain_start - 1 : self.domain_end]

    def window(self, flank: int = FLANK) -> tuple[int, int]:
        """1-based inclusive window of the domain plus ``flank`` residues."""
        return max(1, self.domain_start - flank), min(len(self.protein), self.domain_end + flank)


@dataclass
class ChimeraRecord:
    record_id: str
    pfam_family: str
    junction_class: str
    class_reason: str
    five: ParentDomain
    three: ParentDomain
    hybrid_residue: str
    """Residue from a junction-spanning codon ('' when the junction falls between codons)."""
    chimeric_protein: str
    """Full translated fusion protein (5' start through the first stop)."""
    domain_start: int
    """1-based first residue of the chimeric domain in ``chimeric_protein``."""
    domain_end: int
    junction_index: int
    """1-based position of the last 5'-derived residue (hybrid included) in the chimera."""
    pilot_indel_residues: int | None
    """Screen's junction indel: + insertion, - deletion, None when clean/undetermined."""
    pilot_indel_sequence: str
    pilot_alignment_parent: str
    pilot_alignment_identity: float | None
    genome_builds: list[str] = field(default_factory=list)
    cohorts: list[str] = field(default_factory=list)
    fusions: list[str] = field(default_factory=list)
    n_rows: int = 0
    n_samples: int = 0
    rows_per_cohort: dict[str, int] = field(default_factory=dict)
    max_split_reads: int | None = None
    total_split_reads: int = 0
    max_paired_end_reads: int | None = None
    artefact_flag_counts: dict[str, dict[str, int]] = field(default_factory=dict)
    artefact_clean_rows: int = 0
    """Rows on which every pilot artefact flag is ``pass``."""
    screen_notes: list[str] = field(default_factory=list)

    @property
    def chimeric_domain(self) -> str:
        return self.chimeric_protein[self.domain_start - 1 : self.domain_end]

    def window(self, flank: int = FLANK) -> tuple[int, int]:
        return max(1, self.domain_start - flank), min(
            len(self.chimeric_protein), self.domain_end + flank
        )

    @property
    def chimeric_domain_flanked(self) -> str:
        start, end = self.window()
        return self.chimeric_protein[start - 1 : end]

    @property
    def fusion(self) -> str:
        return f"{self.five.gene}::{self.three.gene}"

    def origin(self, position: int) -> tuple[str, int | None]:
        """Where chimera residue ``position`` (1-based) came from.

        ``("five", i)`` / ``("three", i)`` give the parent residue; ``("hybrid", None)``
        is the junction-spanning codon's residue.
        """
        five_kept = self.five.junction_residue
        if position <= five_kept:
            return "five", position
        if self.hybrid_residue and position == five_kept + 1:
            return "hybrid", None
        offset = position - five_kept - (1 if self.hybrid_residue else 0)
        return "three", self.three.junction_residue + offset - 1


def record_key(family: str, five_tx: str, three_tx: str, flanked: str) -> str:
    digest = hashlib.sha256(f"{family}|{five_tx}|{three_tx}|{flanked}".encode()).hexdigest()
    return f"chim_{digest[:16]}"


def _int(value: Any) -> int | None:
    try:
        number = int(float(value))
    except (TypeError, ValueError):
        return None
    return number if number >= 0 else None


def build_core(
    row: Mapping[str, Any], source: TranscriptModelSource
) -> tuple[ChimeraRecord | None, str | None]:
    """Rebuild the chimeric domain for one screen row; ``(record, None)`` or ``(None, why)``.

    The rebuilt class must equal the row's ``reconstitution_class``; any
    disagreement is returned as a problem, never silently accepted.
    """
    five_gene, three_gene = row["five_prime_gene"], row["three_prime_gene"]
    annotation = row.get("source_annotation") or None
    five_bp, three_bp = (
        _int(row.get("five_prime_breakpoint")),
        _int(row.get("three_prime_breakpoint")),
    )
    if five_bp is None or three_bp is None:
        return None, "missing breakpoint"
    five = source.resolve(five_gene, source_annotated_transcripts(annotation, five_gene))
    three = source.resolve(three_gene, source_annotated_transcripts(annotation, three_gene))
    five_model, five_problem = source.model(five)
    three_model, three_problem = source.model(three)
    if five_model is None or three_model is None:
        return None, five_problem or three_problem or "no transcript model"
    if "*" in five_model.protein or "*" in three_model.protein:
        # Selenoproteins (UGA = selenocysteine) translate with an internal stop; the
        # translated fusion stops there, so no full chimeric domain can be scored.
        return None, "parent CDS has an internal stop codon (selenocysteine UGA?)"
    prediction = predict_junction(
        five_gene, five_bp, five_model, three_gene, three_bp, three_model, row["genome_build"]
    )
    family = row["pfam_family"]
    call = classify_family(prediction, five_model, three_model, family)
    if call.classification != row["reconstitution_class"]:
        return None, (
            f"rebuilt class {call.classification} != screen class {row['reconstitution_class']}"
        )
    if (
        call.chimeric_domain is None
        or call.five_prime_domain is None
        or call.three_prime_domain is None
        or prediction.chimeric_protein is None
        or prediction.five_prime_last_residue is None
        or prediction.three_prime_first_residue is None
    ):
        return None, "no chimeric domain"
    hybrid = translate(prediction.hybrid_codon) if prediction.hybrid_codon else ""
    five_span, three_span = call.five_prime_domain, call.three_prime_domain
    domain_start = five_span.start_aa
    five_len = prediction.five_prime_last_residue + len(hybrid)
    domain_end = five_len + (three_span.end_aa - prediction.three_prime_first_residue + 1)
    protein = prediction.chimeric_protein
    if protein[domain_start - 1 : domain_end] != call.chimeric_domain:
        return None, "chimeric domain is not a slice of the translated fusion protein"
    record = ChimeraRecord(
        record_id="",
        pfam_family=family,
        junction_class=call.classification,
        class_reason=call.reason,
        five=ParentDomain(
            five_gene,
            five_model.transcript_id,
            five_model.protein,
            five_span.start_aa,
            five_span.end_aa,
            prediction.five_prime_last_residue,
        ),
        three=ParentDomain(
            three_gene,
            three_model.transcript_id,
            three_model.protein,
            three_span.start_aa,
            three_span.end_aa,
            prediction.three_prime_first_residue,
        ),
        hybrid_residue=hybrid,
        chimeric_protein=protein,
        domain_start=domain_start,
        domain_end=domain_end,
        junction_index=five_len,
        pilot_indel_residues=call.junction_indel.length if call.junction_indel else None,
        pilot_indel_sequence=call.junction_indel.residues if call.junction_indel else "",
        pilot_alignment_parent=call.alignment_parent or "",
        pilot_alignment_identity=call.alignment_identity,
    )
    record.record_id = record_key(
        family, record.five.transcript, record.three.transcript, record.chimeric_domain_flanked
    )
    return record, None


def add_row(record: ChimeraRecord, row: Mapping[str, Any]) -> None:
    """Fold one screen row's cohort, read support and artefact flags into ``record``."""
    cohort = str(row.get("cohort") or "")
    record.n_rows += 1
    record.rows_per_cohort[cohort] = record.rows_per_cohort.get(cohort, 0) + 1
    for values, value in (
        (record.cohorts, cohort),
        (record.fusions, row.get("fusion_name") or record.fusion),
        (record.genome_builds, row.get("genome_build") or ""),
    ):
        if value and value not in values:
            values.append(value)
            values.sort()
    split, paired = _int(row.get("tumor_split_reads")), _int(row.get("tumor_paired_end_reads"))
    if split is not None:
        record.total_split_reads += split
        record.max_split_reads = max(record.max_split_reads or 0, split)
    if paired is not None:
        record.max_paired_end_reads = max(record.max_paired_end_reads or 0, paired)
    clean = True
    for column in ARTEFACT_FLAG_COLUMNS:
        value = str(row.get(column) or "not_assessed")
        counts = record.artefact_flag_counts.setdefault(column, {})
        counts[value] = counts.get(value, 0) + 1
        clean = clean and value == "pass"
    record.artefact_clean_rows += int(clean)


def build_records(
    rows: Iterable[Mapping[str, Any]],
    sources: Mapping[str, TranscriptModelSource],
    *,
    progress: Any = None,
) -> tuple[list[ChimeraRecord], Counter[str]]:
    """Records for every translated-class row; ``(records, problem_counts)``.

    ``sources`` maps a genome build to its transcript source. Rows sharing a
    junction (same genes, breakpoints, annotation, build and family) are rebuilt
    once.
    """
    records: dict[str, ChimeraRecord] = {}
    built: dict[tuple[Any, ...], tuple[ChimeraRecord | None, str | None]] = {}
    samples: dict[str, set[tuple[str, str]]] = {}
    problems: Counter[str] = Counter()
    for index, row in enumerate(rows, start=1):
        if row.get("reconstitution_class") not in TRANSLATED_CLASSES:
            continue
        key = tuple(
            row.get(column)
            for column in (
                "five_prime_gene",
                "five_prime_breakpoint",
                "three_prime_gene",
                "three_prime_breakpoint",
                "source_annotation",
                "genome_build",
                "pfam_family",
            )
        )
        if key not in built:
            source = sources.get(str(row.get("genome_build")))
            built[key] = build_core(row, source) if source else (None, "no source for genome build")
        core, problem = built[key]
        if core is None:
            problems[problem or "unknown"] += 1
            continue
        record = records.setdefault(core.record_id, _copy_core(core))
        add_row(record, row)
        samples.setdefault(core.record_id, set()).add(
            (str(row.get("cohort")), str(row.get("sample_id")))
        )
        if progress is not None and index % 2000 == 0:
            progress(f"rows {index}: {len(built)} junctions, {len(records)} records")
    for record_id, record in records.items():
        record.n_samples = len(samples[record_id])
    return sorted(records.values(), key=lambda r: r.record_id), problems


def _copy_core(core: ChimeraRecord) -> ChimeraRecord:
    data = asdict(core)
    return record_from_dict(data)


def record_to_dict(record: ChimeraRecord) -> dict[str, Any]:
    data = asdict(record)
    data["chimeric_domain"] = record.chimeric_domain
    data["chimeric_domain_flanked"] = record.chimeric_domain_flanked
    data["flank"] = FLANK
    data["window"] = list(record.window())
    data["five"]["domain_sequence"] = record.five.domain_sequence
    data["three"]["domain_sequence"] = record.three.domain_sequence
    return data


def record_from_dict(data: Mapping[str, Any]) -> ChimeraRecord:
    derived = {"chimeric_domain", "chimeric_domain_flanked", "flank", "window"}
    fields = {k: v for k, v in data.items() if k not in derived}

    def parent(raw: Mapping[str, Any]) -> ParentDomain:
        return ParentDomain(**{k: v for k, v in raw.items() if k != "domain_sequence"})

    fields["five"] = parent(fields["five"])
    fields["three"] = parent(fields["three"])
    return ChimeraRecord(**fields)


def read_records_json(path: Path) -> list[ChimeraRecord]:
    raw = Path(path).read_bytes()
    payload = json.loads(gzip.decompress(raw) if raw[:2] == b"\x1f\x8b" else raw)
    if payload.get("schema_version") != SCHEMA_VERSION:
        raise ValueError(
            f"{path}: schema {payload.get('schema_version')!r}, expected {SCHEMA_VERSION!r}"
        )
    return [record_from_dict(raw) for raw in payload["records"]]
