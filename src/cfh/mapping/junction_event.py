"""Assess one structural-variant event's junction, one row per Pfam family.

This is the adapter between a cohort's rows and the pure core in
:mod:`cfh.mapping.junction_translation`. Row fields that vary by cohort stay
here, so the core takes only genes, breakpoints, transcripts and a build.
"""

from __future__ import annotations

from typing import Any, Callable, Mapping

from cfh.mapping.junction_domain_classifier import (
    NO_TRANSCRIPT_MODEL,
    ReconstitutionCall,
    classify_family,
)
from cfh.mapping.junction_translation import (
    CONVENTION_THREE_PRIME_BASE_DROPPED,
    MECHANISM_NONE,
    JunctionPrediction,
    TranscriptModel,
    predict_junction,
)
from cfh.mapping.transcript_models import (
    ResolvedTranscript,
    TranscriptModelSource,
    source_annotated_transcripts,
)

NO_SHARED_FAMILY = "no_shared_family"
"""Screen-level outcome: the two transcripts share no (selected) Pfam family,
so there is nothing to classify."""

MODE_ANY_DOMAIN = "any_domain"
MODE_KEY_DOMAIN = "key_domain"

READ_SUPPORT_COLUMNS = (
    "tumor_split_reads",
    "tumor_paired_end_reads",
    "tumor_variant_count",
    "normal_variant_count",
    "breakpoint_precision",
)

JUNCTION_CALL_COLUMNS = (
    "pfam_family",
    "reconstitution_class",
    "class_reason",
    "genome_build",
    "five_prime_transcript",
    "three_prime_transcript",
    "five_prime_transcript_source",
    "three_prime_transcript_source",
    "five_prime_transcript_differs_from_canonical",
    "three_prime_transcript_differs_from_canonical",
    "five_prime_canonical_transcript",
    "three_prime_canonical_transcript",
    "junction_mechanism",
    "translated_frame_status",
    "alternative_convention_frame_status",
    "coordinate_convention",
    "five_prime_break_location",
    "three_prime_break_location",
    "junction_notes",
    "no_product_reason",
    "five_prime_cds_retained",
    "three_prime_cds_start",
    "junction_peptide",
    "junction_cds",
    "five_prime_domain_span",
    "three_prime_domain_span",
    "alignment_parent",
    "alignment_identity",
    "junction_indel_residues",
    "junction_indel_sequence",
    "transcript_notes",
)


def junction_call_columns(
    call: Mapping[str, Any] | None, row: Mapping[str, Any] | None = None
) -> dict[str, Any]:
    """Flat event-row columns, present (blank) on every row so TSV output keeps
    one column set. ``reconstitution_class`` keeps its name; the rest are
    prefixed ``reconstitution_``.
    """
    columns: dict[str, Any] = {}
    for name in JUNCTION_CALL_COLUMNS:
        key = name if name == "reconstitution_class" else f"reconstitution_{name}"
        columns[key] = call.get(name) if call else None
    support = read_support_columns(row) if row is not None else {}
    for name in READ_SUPPORT_COLUMNS:
        columns[f"reconstitution_{name}"] = support.get(name)
    return columns


def read_support_columns(row: Mapping[str, Any]) -> dict[str, Any]:
    """Read support and precision. Negative or absent counts become blank."""

    def count(column: str) -> int | None:
        value = row.get(column)
        try:
            number = int(value)  # type: ignore[arg-type]
        except (TypeError, ValueError):
            return None
        return number if number >= 0 else None

    precision = str(row.get("Breakpoint_Type") or "").strip().upper()
    return {
        "tumor_split_reads": count("Tumor_Split_Read_Count"),
        "tumor_paired_end_reads": count("Tumor_Paired_End_Read_Count"),
        "tumor_variant_count": count("Tumor_Variant_Count"),
        "normal_variant_count": count("Normal_Variant_Count"),
        # MSK's archive spells the imprecise class both ways.
        "breakpoint_precision": precision.replace("IMPPRECISE", "IMPRECISE") or None,
    }


def _span_text(span) -> str | None:
    if span is None:
        return None
    return f"{span.start_aa}-{span.end_aa} kept {span.retained_start_aa}-{span.retained_end_aa}"


def _flat_row(
    family: str | None,
    classification: str,
    reason: str,
    genome_build: str | None,
    five: ResolvedTranscript,
    three: ResolvedTranscript,
    prediction: JunctionPrediction | None = None,
    alternative_frame: str | None = None,
    call: ReconstitutionCall | None = None,
) -> dict[str, Any]:
    notes = list(prediction.notes) if prediction else []
    five_location = prediction.five_prime_location if prediction else None
    three_location = prediction.three_prime_location if prediction else None
    indel = call.junction_indel if call else None
    return {
        "pfam_family": family,
        "reconstitution_class": classification,
        "class_reason": reason,
        "genome_build": genome_build,
        "five_prime_transcript": five.refseq_id or five.transcript_id,
        "three_prime_transcript": three.refseq_id or three.transcript_id,
        "five_prime_transcript_source": five.source,
        "three_prime_transcript_source": three.source,
        "five_prime_transcript_differs_from_canonical": five.differs_from_canonical,
        "three_prime_transcript_differs_from_canonical": three.differs_from_canonical,
        "five_prime_canonical_transcript": five.canonical_transcript_id,
        "three_prime_canonical_transcript": three.canonical_transcript_id,
        "junction_mechanism": prediction.mechanism if prediction else MECHANISM_NONE,
        "translated_frame_status": prediction.frame_status if prediction else None,
        "alternative_convention_frame_status": alternative_frame,
        "coordinate_convention": prediction.coordinate_convention if prediction else None,
        "five_prime_break_location": five_location.describes if five_location else None,
        "three_prime_break_location": three_location.describes if three_location else None,
        "junction_notes": "; ".join(notes) or None,
        "no_product_reason": prediction.no_product_reason if prediction else reason,
        "five_prime_cds_retained": prediction.five_prime_cds_retained if prediction else None,
        "three_prime_cds_start": prediction.three_prime_cds_start if prediction else None,
        "junction_peptide": prediction.junction_peptide if prediction else None,
        "junction_cds": prediction.junction_cds if prediction else None,
        "five_prime_domain_span": _span_text(call.five_prime_domain) if call else None,
        "three_prime_domain_span": _span_text(call.three_prime_domain) if call else None,
        "alignment_parent": call.alignment_parent if call else None,
        "alignment_identity": call.alignment_identity if call else None,
        "junction_indel_residues": indel.length if indel else None,
        "junction_indel_sequence": indel.residues if indel else None,
        "transcript_notes": "; ".join(n for n in (five.note, three.note) if n) or None,
    }


def assess_event(
    *,
    five_prime_gene: str,
    five_prime_breakpoint: int | None,
    three_prime_gene: str,
    three_prime_breakpoint: int | None,
    annotation: str | None,
    genome_build: str | None,
    source: TranscriptModelSource,
    family_filter: Callable[[set[str]], set[str]] | None = None,
) -> list[dict[str, Any]]:
    """One flat row per selected Pfam family shared by the two transcripts.

    Pfam families come from the transcripts actually used, so a family present
    only on another isoform is not shared here. With no shared family this
    returns one :data:`NO_SHARED_FAMILY` row. CDS sequence is only fetched when
    a family is shared.

    Network failures propagate, so a flaky lookup is retried rather than
    recorded as a biological result.
    """
    five = source.resolve(
        five_prime_gene, source_annotated_transcripts(annotation, five_prime_gene)
    )
    three = source.resolve(
        three_prime_gene, source_annotated_transcripts(annotation, three_prime_gene)
    )
    shared = five.families & three.families
    if family_filter is not None:
        shared = family_filter(shared)
    if not shared:
        reason = (
            f"{five_prime_gene} ({five.refseq_id or five.transcript_id}) and {three_prime_gene} "
            f"({three.refseq_id or three.transcript_id}) share no selected Pfam family"
        )
        if five.transcript is None or three.transcript is None:
            reason = five.note or three.note or reason
        return [_flat_row(None, NO_SHARED_FAMILY, reason, genome_build, five, three)]

    five_model: TranscriptModel | None = None
    three_model: TranscriptModel | None = None
    five_problem: str | None = None
    three_problem: str | None = None
    if five_prime_breakpoint is None or three_prime_breakpoint is None:
        missing = five_prime_gene if five_prime_breakpoint is None else three_prime_gene
        reason = f"no breakpoint within the {missing} transcript's locus"
        return [
            _flat_row(family, NO_TRANSCRIPT_MODEL, reason, genome_build, five, three)
            for family in sorted(shared)
        ]
    five_model, five_problem = source.model(five)
    three_model, three_problem = source.model(three)
    problem = five_problem or three_problem
    if problem:
        return [
            _flat_row(family, NO_TRANSCRIPT_MODEL, problem, genome_build, five, three)
            for family in sorted(shared)
        ]
    prediction = predict_junction(
        five_prime_gene,
        five_prime_breakpoint,
        five_model,
        three_prime_gene,
        three_prime_breakpoint,
        three_model,
        genome_build,
    )
    alternative = predict_junction(
        five_prime_gene,
        five_prime_breakpoint,
        five_model,
        three_prime_gene,
        three_prime_breakpoint,
        three_model,
        genome_build,
        coordinate_convention=CONVENTION_THREE_PRIME_BASE_DROPPED,
    )
    alternative_frame = alternative.frame_status if alternative.has_product else None
    rows = []
    for family in sorted(shared):
        call = classify_family(prediction, five_model, three_model, family)
        rows.append(
            _flat_row(
                family,
                call.classification,
                call.reason,
                genome_build,
                five,
                three,
                prediction,
                alternative_frame,
                call,
            )
        )
    return rows
