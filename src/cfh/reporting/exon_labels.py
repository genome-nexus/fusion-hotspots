"""Shared, gene-agnostic phrasing for stating a breakpoint's protein exon.

Every report surface that states a specific breakpoint's amino-acid
position also states the protein exon it falls in, so a reader never has
to cross-reference the SVG schematics to know which exon a stated "aa 439"
corresponds to. Two distinct situations are handled, deliberately kept
separate rather than unified into one "guess an exon" helper:

1. A REAL, already-mapped event has its own ``Breakpoint_exon`` (and
   ``is_intronic_breakpoint``) computed once, per event, by
   ``cfh.mapping.transcript_source``/``cfh.real_benchmark`` and carried on
   its ``FusionFeature``/results row. :func:`format_breakpoint_exon`
   renders that already-computed value verbatim -- it is never re-derived
   here, since the per-event mapping already resolved any annotation-vs-
   Genome-Nexus-fallback ambiguity that a from-scratch re-derivation could
   get wrong.

2. A statistically inferred, not-tied-to-one-event position (e.g. the
   cutpoint-detection algorithm's inferred cutpoint, which is the best
   dividing candidate among many observed positions, not the mapped
   position of any single reported breakpoint) has no such per-event
   field to reuse. :func:`exon_label_for_protein_position` instead reads
   the gene's own exon-to-residue boundary structure
   (``gene_track["exon_boundaries_aa"]``, the same field
   ``cfh.reporting.fusion_schematic.exon_boundary_ticks_svg`` already uses
   for the schematics' exon-tick axis labels) to say which exon a protein
   position falls in -- and, when the position sits exactly on the shared
   coordinate between two consecutive exons (the residue-rounding artifact
   produced by a clamped intronic-breakpoint estimate; see
   ``map_genomic_breakpoint_to_protein_position``), names both flanking
   exons rather than fabricating a single one.
"""

from __future__ import annotations

from typing import Any, Optional


def format_breakpoint_exon(
    breakpoint_exon: Optional[int], is_intronic: Optional[bool] = None
) -> str:
    """Render one real event's already-computed ``Breakpoint_exon``.

    ``breakpoint_exon`` being ``None`` means the exon genuinely could not be
    resolved for this breakpoint (e.g. only the Ensembl protein-feature
    fallback was available, which has no exon-to-protein-coordinate map) --
    rendered as an explicit "exon unavailable", never omitted or invented.
    ``is_intronic`` is the paired ``is_intronic_breakpoint`` flag: when
    true, ``breakpoint_exon`` is the *nearest* exon to a breakpoint that
    genomically falls in an intron, not an exon the breakpoint literally
    sits inside (see ``map_genomic_breakpoint_to_protein_position``), so
    this is phrased as an approximation rather than a plain "exon N".
    """
    if breakpoint_exon is None:
        return "exon unavailable"
    if is_intronic:
        return f"intronic breakpoint, nearest exon {breakpoint_exon}"
    return f"exon {breakpoint_exon}"


def exon_label_for_protein_position(
    exon_boundaries: Optional[list[dict[str, Any]]], position: Optional[int]
) -> str:
    """Render which exon(s) a protein-residue ``position`` falls in/at.

    Reads ``exon_boundaries`` (``gene_track["exon_boundaries_aa"]``, a list
    of ``{"exon_rank", "start_aa", "end_aa"}`` dicts), not any single
    event's ``Breakpoint_exon`` -- for a position with no one real mapped
    event to point to (see the module docstring).

    Returns "exon N" when ``position`` falls strictly inside one exon's
    span; "exon N/M boundary" when it sits exactly on the shared coordinate
    between two consecutive exons (both spans contain it because their
    boundaries coincide there); and "exon position unavailable" when
    ``position``/``exon_boundaries`` is missing or no configured exon span
    contains ``position`` at all.
    """
    if position is None or not exon_boundaries:
        return "exon position unavailable"
    containing_ranks = sorted(
        {
            boundary["exon_rank"]
            for boundary in exon_boundaries
            if boundary.get("start_aa") is not None
            and boundary.get("end_aa") is not None
            and boundary["start_aa"] <= position <= boundary["end_aa"]
        }
    )
    if not containing_ranks:
        return "exon position unavailable"
    if len(containing_ranks) == 1:
        return f"exon {containing_ranks[0]}"
    return "exon " + "/".join(str(rank) for rank in containing_ranks) + " boundary"
