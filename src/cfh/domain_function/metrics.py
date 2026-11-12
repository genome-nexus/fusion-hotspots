"""Sequence/structure metrics of one chimeric domain against its two parents.

Coordinates: parents and the chimera use their own 1-based residue numbers.
The two parent windows (domain +/- :data:`~cfh.domain_function.record.FLANK`)
are aligned globally (Biopython ``PairwiseAligner``, BLOSUM62, gap open -11,
extend -1, free end gaps); every parent residue in its window gets an
alignment column.

Placing the chimera on those columns (``ChimeraColumns``): 5'-derived residues
take their 5' parent column; the hybrid-codon residue takes the 3' parent's
column for residue ``F-1``; 3'-derived residues take their 3' parent column.
Any residue whose column is not after the previous chimera residue's column
cannot be placed -- it is a *junction-inserted* residue (a tandem duplicate of
aligned columns). Columns skipped between the last 5'-side and first 3'-side
residue are *junction-deleted*. Standard SCHEMA assumes crossovers between
aligned blocks; inserted and deleted residues are therefore counted separately
rather than forced into the block model.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache
from typing import Any

from Bio.Align import PairwiseAligner, substitution_matrices

from cfh.domain_function.afdb import CONFIDENT_PLDDT, ParentStructure, ss_class
from cfh.domain_function.record import ChimeraRecord, ParentDomain
from cfh.domain_function.uniprot import POINT_FEATURE_TYPES, UniProtEntry

NOT_ASSESSED = "not_assessed"
KEY_RESIDUE_PROXIMITY = 5
"""Junction/indel within this many parent residues of a key site is reported."""
MAX_SITE_SPAN = 10
"""Motif scale. Spans of at most this many residues count as key *sites*: for junction
proximity (a DNA-binding region covering the whole domain would make every junction
'near') and for tiering a lost/altered span (rule U9). Longer regions are reported
in ``key_span_lost_or_altered`` only."""
MIN_UNIPROT_IDENTITY = 0.9
"""Minimum transcript-protein vs UniProt identity over the parent window to map features."""
JUNCTION_HALF_WINDOW = 3
"""Junction residues = this many chimera residues each side of the junction, plus inserts."""


def _aligner() -> PairwiseAligner:
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aligner.open_end_gap_score = 0
    aligner.extend_end_gap_score = 0
    return aligner


_ALIGNER = _aligner()


def align_map(first: str, second: str) -> list[tuple[int | None, int | None]]:
    """Alignment columns of two sequences as (0-based index or None) pairs."""
    alignment = _ALIGNER.align(first, second)[0]
    indices = alignment.indices
    return [
        (int(a) if a >= 0 else None, int(b) if b >= 0 else None)
        for a, b in zip(indices[0], indices[1], strict=True)
    ]


@dataclass
class ParentAlignment:
    """Column layout of the two parent windows."""

    columns: list[tuple[int | None, int | None]]
    """Per column: (5' parent residue, 3' parent residue), 1-based, None for a gap."""
    five_col: dict[int, int]
    three_col: dict[int, int]

    def residue(self, parent: ParentDomain, which: str, column: int) -> str:
        position = self.columns[column][0 if which == "five" else 1]
        return parent.protein[position - 1] if position is not None else "-"


def align_parents(record: ChimeraRecord) -> ParentAlignment:
    f_start, f_end = record.five.window()
    t_start, t_end = record.three.window()
    pairs = align_map(
        record.five.protein[f_start - 1 : f_end], record.three.protein[t_start - 1 : t_end]
    )
    columns = [
        (a + f_start if a is not None else None, b + t_start if b is not None else None)
        for a, b in pairs
    ]
    five_col = {a: c for c, (a, _) in enumerate(columns) if a is not None}
    three_col = {b: c for c, (_, b) in enumerate(columns) if b is not None}
    return ParentAlignment(columns, five_col, three_col)


@dataclass
class ChimeraColumns:
    chimera_col: dict[int, int]
    """Chimera position -> column, for placed residues inside the parent windows."""
    col_residue: dict[int, str]
    inserted: list[int]
    """Chimera positions of junction-inserted residues."""
    inserted_sources: list[tuple[str, int]]
    """(parent, parent position) each inserted residue was copied from."""
    deleted_columns: list[int]
    deleted_residues: dict[str, list[int]] = field(default_factory=dict)
    """Parent positions in deleted columns, per parent ('five'/'three')."""


def place_chimera(record: ChimeraRecord, parents: ParentAlignment) -> ChimeraColumns:
    start, end = record.window()
    chimera_col: dict[int, int] = {}
    inserted: list[int] = []
    sources: list[tuple[str, int]] = []
    last_col = -1
    last_five_side_col = -1
    first_three_col: int | None = None
    for position in range(start, end + 1):
        which, parent_position = record.origin(position)
        if which == "five":
            column = parents.five_col.get(parent_position)  # type: ignore[arg-type]
        elif which == "three":
            column = parents.three_col.get(parent_position)  # type: ignore[arg-type]
        else:
            column = parents.three_col.get(record.three.junction_residue - 1)
        if column is None:
            continue  # outside the aligned parent windows
        if column <= last_col:
            inserted.append(position)
            source_parent = "three" if which in ("three", "hybrid") else "five"
            source_position = (
                parent_position
                if parent_position is not None
                else record.three.junction_residue - 1
            )
            sources.append((source_parent, source_position))
            continue
        chimera_col[position] = column
        last_col = column
        if which in ("five", "hybrid"):
            last_five_side_col = column
        elif first_three_col is None:
            first_three_col = column
    deleted = (
        list(range(last_five_side_col + 1, first_three_col))
        if first_three_col is not None and last_five_side_col >= 0
        else []
    )
    deleted_residues: dict[str, list[int]] = {"five": [], "three": []}
    for column in deleted:
        five_position, three_position = parents.columns[column]
        if five_position is not None:
            deleted_residues["five"].append(five_position)
        if three_position is not None:
            deleted_residues["three"].append(three_position)
    col_residue = {
        column: record.chimeric_protein[position - 1] for position, column in chimera_col.items()
    }
    return ChimeraColumns(chimera_col, col_residue, inserted, sources, deleted, deleted_residues)


@lru_cache(maxsize=4096)
def _sequence_map(reference: str, protein: str) -> dict[int, int]:
    pairs = align_map(reference, protein)
    return {a + 1: b + 1 for a, b in pairs if a is not None and b is not None}


def uniprot_position_map(parent: ParentDomain, entry: UniProtEntry) -> tuple[dict[int, int], float]:
    """(UniProt position -> parent position, identity over the parent window)."""
    mapping = _sequence_map(entry.sequence, parent.protein)
    w_start, w_end = parent.window()
    window = [(u, p) for u, p in mapping.items() if w_start <= p <= w_end]
    same = sum(1 for u, p in window if entry.sequence[u - 1] == parent.protein[p - 1])
    identity = same / (w_end - w_start + 1)
    return mapping, identity


STATUS_ORDER = ("retained", "conserved", "swapped", "altered", "unresolved", "lost")
"""Per annotated parent residue, least to most severe:

* retained: the parent's own residue is in the chimera (decided from the
  junction coordinates, not the alignment);
* conserved: it is absent, but the residue now at its aligned column is the
  same amino acid (supplied by the other parent);
* swapped: absent; the column holds the other parent's different residue;
* altered: absent; the column holds a residue matching neither parent (the
  hybrid-codon residue), or a present residue differs from the parent;
* unresolved: the UniProt position does not align to the parent transcript, so
  survival cannot be decided (mapping uncertainty -- never counted as retained);
* lost: absent, and no chimera residue occupies its aligned column (a
  junction-deleted column, or a column where the other parent has a gap).

A feature takes its worst residue status."""

KEY_ANNOTATION_OK = ("annotated", "no_annotations_in_window")
"""Per-parent annotation statuses that count as assessed."""


def _column_residue(
    record: ChimeraRecord,
    parents: ParentAlignment,
    placed: ChimeraColumns,
    other: str,
    column: int,
) -> str | None:
    """The chimera residue at alignment ``column``, from whichever source holds it."""
    residue = placed.col_residue.get(column)
    if residue is not None:
        return residue
    other_position = parents.columns[column][0 if other == "five" else 1]
    if other_position is None:
        return None
    position = record.chimera_position(other, other_position)
    return record.chimeric_protein[position - 1] if position is not None else None


def residue_status(
    record: ChimeraRecord,
    parents: ParentAlignment,
    placed: ChimeraColumns,
    which: str,
    position: int,
) -> str:
    """Status (see :data:`STATUS_ORDER`) of parent residue ``position`` of ``which``."""
    parent = record.five if which == "five" else record.three
    other = "three" if which == "five" else "five"
    own = parent.protein[position - 1]
    present = record.chimera_position(which, position)
    if present is not None:
        return "retained" if record.chimeric_protein[present - 1] == own else "altered"
    column = (parents.five_col if which == "five" else parents.three_col).get(position)
    if column is None:
        return "unresolved"
    if column in set(placed.deleted_columns):
        return "lost"
    residue = _column_residue(record, parents, placed, other, column)
    if residue is None:
        return "lost"
    if residue == own:
        return "conserved"
    partner = parents.columns[column][0 if other == "five" else 1]
    partner_parent = record.five if other == "five" else record.three
    if partner is not None and partner_parent.protein[partner - 1] == residue:
        return "swapped"
    return "altered"


def _empty_parent_report(report: dict[str, Any], prefix: str, status: str) -> None:
    report[f"key_annotation_status_{prefix}"] = status
    report[f"key_features_{prefix}"] = NOT_ASSESSED
    report[f"junction_to_key_{prefix}"] = NOT_ASSESSED


def key_residue_report(
    record: ChimeraRecord,
    parents: ParentAlignment,
    placed: ChimeraColumns,
    entries: dict[str, UniProtEntry | None],
) -> dict[str, Any]:
    """Key-feature survival for both parents.

    Each parent gets ``key_annotation_status_*``: ``annotated`` /
    ``no_annotations_in_window`` (assessed) or ``no_uniprot_entry`` /
    ``low_uniprot_identity`` / ``incomplete_annotation`` (not assessed).
    ``key_features_assessed`` is True only when both parents are assessed.
    Features are evaluated over the parent window (domain +/- FLANK, in UniProt
    coordinates the span between the first and last window residues that align).
    """
    report: dict[str, Any] = {}
    statuses: dict[str, list[str]] = {s: [] for s in STATUS_ORDER}
    point_lost: list[str] = []
    point_altered: list[str] = []
    span_lost_altered: list[str] = []
    motif_lost_altered: list[str] = []
    incomplete_labels: list[str] = []
    assessed_all = True
    for which, parent in (("five", record.five), ("three", record.three)):
        prefix = "5p" if which == "five" else "3p"
        entry = entries.get(which)
        if entry is None:
            _empty_parent_report(report, prefix, "no_uniprot_entry")
            assessed_all = False
            continue
        mapping, identity = uniprot_position_map(parent, entry)
        report[f"uniprot_{prefix}"] = entry.accession
        report[f"uniprot_identity_{prefix}"] = round(identity, 3)
        if identity < MIN_UNIPROT_IDENTITY:
            _empty_parent_report(report, prefix, "low_uniprot_identity")
            assessed_all = False
            continue
        w_start, w_end = parent.window()
        in_window = [u for u, p in mapping.items() if w_start <= p <= w_end]
        u_low, u_high = min(in_window), max(in_window)
        relevant_incomplete = [f for f in entry.incomplete if f.may_overlap(u_low, u_high)]
        if relevant_incomplete:
            incomplete_labels += [f"{parent.gene}:{f.label}" for f in relevant_incomplete]
            _empty_parent_report(report, prefix, "incomplete_annotation")
            assessed_all = False
            continue
        n_features = 0
        key_positions: list[int] = []
        for feature in entry.features:
            uniprot_positions = range(max(feature.start, u_low), min(feature.end, u_high) + 1)
            if not uniprot_positions:
                continue
            n_features += 1
            worst = "retained"
            for u in uniprot_positions:
                if u not in mapping:
                    status = "unresolved"
                else:
                    status = residue_status(record, parents, placed, which, mapping[u])
                    if (
                        feature.type in POINT_FEATURE_TYPES
                        or feature.end - feature.start < MAX_SITE_SPAN
                    ):
                        key_positions.append(mapping[u])
                if STATUS_ORDER.index(status) > STATUS_ORDER.index(worst):
                    worst = status
            label = f"{parent.gene}:{feature.label}"
            statuses[worst].append(label)
            point = feature.type in POINT_FEATURE_TYPES
            if point and worst == "lost":
                point_lost.append(label)
            if point and worst == "altered":
                point_altered.append(label)
            if not point and worst in ("lost", "altered"):
                span_lost_altered.append(label)
                if feature.end - feature.start < MAX_SITE_SPAN:
                    motif_lost_altered.append(label)
        report[f"key_annotation_status_{prefix}"] = (
            "annotated" if n_features else "no_annotations_in_window"
        )
        report[f"key_features_{prefix}"] = n_features
        junction = parent.junction_residue
        report[f"junction_to_key_{prefix}"] = (
            min(abs(p - junction) for p in key_positions) if key_positions else None
        )
    report["key_features_assessed"] = assessed_all
    for status in STATUS_ORDER:
        report[f"key_{status}"] = len(statuses[status]) if assessed_all else NOT_ASSESSED
        report[f"key_{status}_list"] = ";".join(statuses[status])
    report["key_point_lost"] = ";".join(point_lost)
    report["key_point_altered"] = ";".join(point_altered)
    report["key_span_lost_or_altered"] = ";".join(span_lost_altered)
    report["key_motif_lost_or_altered"] = ";".join(motif_lost_altered)
    report["key_incomplete_annotations"] = ";".join(incomplete_labels)
    distances: list[int] = [
        value
        for value in (report.get("junction_to_key_5p"), report.get("junction_to_key_3p"))
        if isinstance(value, int)
    ]
    near = bool(distances) and min(distances) <= KEY_RESIDUE_PROXIMITY
    report["junction_near_key"] = near if assessed_all else NOT_ASSESSED
    return report


@dataclass
class MappedStructure:
    structure: ParentStructure
    to_parent: dict[int, int]
    """UniProt position -> parent transcript position."""
    from_parent: dict[int, int]


STRUCTURE_OK = "ok"


def map_structure(
    parent: ParentDomain,
    entry: UniProtEntry | None,
    structure: ParentStructure | None,
    fetch_problem: str = "",
) -> tuple[MappedStructure | None, str, str]:
    """``(mapped, status, detail)``; ``mapped`` is None unless ``status`` is ``ok``.

    Statuses: ``ok``, ``no_dssp`` (mapped, but the model has no secondary
    structure annotation), ``no_uniprot_entry``, ``afdb_fetch_failed``,
    ``no_afdb_model``, ``afdb_sequence_mismatch`` (model built on another
    sequence; rejected), ``low_uniprot_identity``.
    """
    if entry is None:
        return None, "no_uniprot_entry", "no reviewed UniProt entry"
    if structure is None:
        if fetch_problem:
            return None, "afdb_fetch_failed", fetch_problem
        return None, "no_afdb_model", "AlphaFold DB has no usable model"
    if structure.sequence != entry.sequence:
        return None, "afdb_sequence_mismatch", "AlphaFold DB model sequence differs from UniProt"
    mapping, identity = uniprot_position_map(parent, entry)
    if identity < MIN_UNIPROT_IDENTITY:
        return None, "low_uniprot_identity", f"transcript vs UniProt identity {identity:.2f}"
    mapped = MappedStructure(structure, mapping, {p: u for u, p in mapping.items()})
    if structure.ss is None:
        return mapped, "no_dssp", "model has no _struct_conf secondary structure"
    return mapped, STRUCTURE_OK, ""


def _ss_codes(mapped: MappedStructure, positions: list[int]) -> list[str | None]:
    return [
        mapped.structure.ss_at(mapped.from_parent[p]) if p in mapped.from_parent else None
        for p in positions
    ]


def junction_structure_report(
    record: ChimeraRecord,
    placed: ChimeraColumns,
    structures: dict[str, MappedStructure | None],
) -> dict[str, Any]:
    """Secondary structure (DSSP via AFDB) at the junction and at junction indel residues."""
    report: dict[str, Any] = {}
    involved: list[str | None] = []
    plddts: list[float] = []
    indel_present = bool(placed.inserted or placed.deleted_columns)
    for which, parent in (("five", record.five), ("three", record.three)):
        prefix = "5p" if which == "five" else "3p"
        mapped = structures.get(which)
        j = parent.junction_residue
        boundary = [j, j + 1] if which == "five" else [j - 1, j]
        indel_positions = list(placed.deleted_residues.get(which, []))
        indel_positions += [p for source, p in placed.inserted_sources if source == which]
        if mapped is None or mapped.structure.ss is None:
            report[f"junction_ss_{prefix}"] = NOT_ASSESSED
            if indel_positions:
                involved.append(None)
            continue
        codes = _ss_codes(mapped, boundary)
        classes = sorted({ss_class(code) for code in codes})
        report[f"junction_ss_{prefix}"] = "|".join(classes)
        context = list(range(j - 3, j + 4))
        report[f"junction_dssp_{prefix}"] = "".join(
            code or "?" for code in _ss_codes(mapped, context)
        )
        for code in _ss_codes(mapped, indel_positions):
            involved.append(code)
        for position in indel_positions + boundary:
            uniprot = mapped.from_parent.get(position)
            value = mapped.structure.plddt_at(uniprot) if uniprot else None
            if value is not None:
                plddts.append(value)
    if not indel_present:
        report["indel_ss"] = "no_indel"
        report["indel_in_ss"] = False
    elif any(code is None for code in involved) or not involved:
        report["indel_ss"] = NOT_ASSESSED
        report["indel_in_ss"] = NOT_ASSESSED
    else:
        classes = sorted({ss_class(code) for code in involved})
        report["indel_ss"] = "|".join(classes)
        report["indel_in_ss"] = any(c in ("helix", "strand") for c in classes)
    report["junction_min_plddt"] = round(min(plddts), 1) if plddts else None
    return report


def schema_report(
    record: ChimeraRecord,
    parents: ParentAlignment,
    placed: ChimeraColumns,
    structures: dict[str, MappedStructure | None],
) -> dict[str, Any]:
    """SCHEMA-like broken contacts (Voigt et al. 2002) per parent structure.

    For each parent model, contacts between two residues of that parent's Pfam
    domain, both with pLDDT >= 70, are placed on the alignment. A contact is
    *broken* (counts toward E) when the chimera's residue pair at those columns
    matches neither parent's pair. Contacts with a residue in a junction-deleted
    column are counted as ``lost_to_deletion`` instead; contacts touching the
    source copy of a junction-inserted residue as ``touching_inserted``.
    """
    report: dict[str, Any] = {}
    deleted = set(placed.deleted_columns)
    for which, parent in (("five", record.five), ("three", record.three)):
        prefix = "5p" if which == "five" else "3p"
        mapped = structures.get(which)
        keys = (
            "E",
            "E_core",
            "contacts",
            "lost_to_deletion",
            "touching_inserted",
            "plddt_coverage",
        )
        if mapped is None:
            for key in keys:
                report[f"schema_{key}_{prefix}"] = NOT_ASSESSED
            continue
        cols = parents.five_col if which == "five" else parents.three_col
        domain = range(parent.domain_start, parent.domain_end + 1)
        confident = {
            p
            for p in domain
            if p in mapped.from_parent
            and (mapped.structure.plddt_at(mapped.from_parent[p]) or 0.0) >= CONFIDENT_PLDDT
        }
        inserted_sources = {p for source, p in placed.inserted_sources if source == which}
        energy = core = total = lost = touching = 0
        for u_a, u_b in mapped.structure.contacts:
            a, b = mapped.to_parent.get(u_a), mapped.to_parent.get(u_b)
            if a not in confident or b not in confident:
                continue
            assert a is not None and b is not None
            total += 1
            if a in inserted_sources or b in inserted_sources:
                touching += 1
            col_a, col_b = cols.get(a), cols.get(b)
            if col_a is None or col_b is None:
                continue
            if col_a in deleted or col_b in deleted:
                lost += 1
                continue
            chimera = (placed.col_residue.get(col_a), placed.col_residue.get(col_b))
            if None in chimera:
                continue
            five_pair = (
                parents.residue(record.five, "five", col_a),
                parents.residue(record.five, "five", col_b),
            )
            three_pair = (
                parents.residue(record.three, "three", col_a),
                parents.residue(record.three, "three", col_b),
            )
            if chimera != five_pair and chimera != three_pair:
                energy += 1
                if mapped.structure.buried(u_a) and mapped.structure.buried(u_b):
                    core += 1
        report[f"schema_E_{prefix}"] = energy
        report[f"schema_E_core_{prefix}"] = core
        report[f"schema_contacts_{prefix}"] = total
        report[f"schema_lost_to_deletion_{prefix}"] = lost
        report[f"schema_touching_inserted_{prefix}"] = touching
        report[f"schema_plddt_coverage_{prefix}"] = round(len(confident) / len(domain), 3)
    return report
