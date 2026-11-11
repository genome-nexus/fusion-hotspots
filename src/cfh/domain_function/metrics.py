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
"""Proximity uses point features and spans of at most this many residues (e.g. motifs);
a DNA-binding region covering the whole domain would make every junction 'near'."""
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


STATUS_ORDER = ("retained", "swapped", "altered", "lost")
"""retained: chimera keeps the parent residue; swapped: it has the other parent's
aligned residue; altered: differs from both parents; lost: column deleted."""


def key_residue_report(
    record: ChimeraRecord,
    parents: ParentAlignment,
    placed: ChimeraColumns,
    entries: dict[str, UniProtEntry | None],
) -> dict[str, Any]:
    """Key-feature retention for both parents; fields are NOT_ASSESSED when data is missing."""
    report: dict[str, Any] = {}
    statuses: dict[str, list[str]] = {s: [] for s in STATUS_ORDER}
    point_lost: list[str] = []
    point_altered: list[str] = []
    assessed_all = True
    deleted = set(placed.deleted_columns)
    for which, parent in (("five", record.five), ("three", record.three)):
        other = "three" if which == "five" else "five"
        prefix = "5p" if which == "five" else "3p"
        entry = entries.get(which)
        if entry is None:
            report[f"uniprot_{prefix}"] = NOT_ASSESSED
            report[f"key_features_{prefix}"] = NOT_ASSESSED
            report[f"junction_to_key_{prefix}"] = NOT_ASSESSED
            assessed_all = False
            continue
        mapping, identity = uniprot_position_map(parent, entry)
        report[f"uniprot_{prefix}"] = entry.accession
        report[f"uniprot_identity_{prefix}"] = round(identity, 3)
        if identity < MIN_UNIPROT_IDENTITY:
            report[f"key_features_{prefix}"] = NOT_ASSESSED
            report[f"junction_to_key_{prefix}"] = NOT_ASSESSED
            assessed_all = False
            continue
        cols = parents.five_col if which == "five" else parents.three_col
        w_start, w_end = parent.window()
        n_features = 0
        key_positions: list[int] = []
        for feature in entry.features:
            positions = [
                mapping[u]
                for u in range(feature.start, feature.end + 1)
                if u in mapping and w_start <= mapping[u] <= w_end
            ]
            if not positions:
                continue
            n_features += 1
            if feature.type in POINT_FEATURE_TYPES or feature.end - feature.start < MAX_SITE_SPAN:
                key_positions.extend(positions)
            worst = "retained"
            for position in positions:
                column = cols.get(position)
                if column is None:
                    continue
                if column in deleted:
                    status = "lost"
                else:
                    chimera = placed.col_residue.get(column)
                    if chimera is None:
                        continue  # beyond the chimera's window
                    own = parent.protein[position - 1]
                    partner = parents.residue(
                        record.three if which == "five" else record.five, other, column
                    )
                    status = (
                        "retained"
                        if chimera == own
                        else "swapped"
                        if chimera == partner
                        else "altered"
                    )
                if STATUS_ORDER.index(status) > STATUS_ORDER.index(worst):
                    worst = status
            label = f"{parent.gene}:{feature.label}"
            statuses[worst].append(label)
            if feature.type in POINT_FEATURE_TYPES and worst == "lost":
                point_lost.append(label)
            if feature.type in POINT_FEATURE_TYPES and worst == "altered":
                point_altered.append(label)
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


def map_structure(
    parent: ParentDomain, entry: UniProtEntry | None, structure: ParentStructure | None
) -> tuple[MappedStructure | None, str]:
    """The structure on parent coordinates, or (None, reason)."""
    if entry is None:
        return None, "no UniProt entry"
    if structure is None:
        return None, "no AlphaFold DB model"
    if structure.sequence and structure.sequence != entry.sequence:
        return None, "AlphaFold DB model sequence differs from UniProt"
    mapping, identity = uniprot_position_map(parent, entry)
    if identity < MIN_UNIPROT_IDENTITY:
        return None, f"transcript vs UniProt identity {identity:.2f} in window"
    return MappedStructure(structure, mapping, {p: u for u, p in mapping.items()}), ""


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
