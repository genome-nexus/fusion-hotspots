"""Per-record confidence and comparison metrics of a predicted chimeric domain.

Coordinates are the chimera's and each parent's own 1-based residue numbers
(as in :mod:`cfh.domain_function`); a prediction covers its window starting at
``*_start``, so residue ``p`` is model index ``p - start``.

Definitions (see ``docs/domain_structure.md``):

* **domain** -- the Pfam domain span (chimera: ``domain_start..domain_end``;
  parent: its own Pfam span), not the +/- flank.
* **junction bond** -- between chimera residues ``J`` and ``J+1``
  (``J = junction_index``, the last 5'-derived residue, hybrid included); in
  the 5' parent between ``F`` and ``F+1``, in the 3' parent between ``T-1``
  and ``T``.
* **junction window** -- the 10 residues on each side of the bond
  (``J-9..J+10``), widened to keep 10 residues beyond any junction-inserted
  residue; clipped to the predicted window.
* **halves** -- chimera domain residues from the 5' parent vs from the 3'
  parent (the hybrid-codon residue is in neither). A parent's halves are its
  domain residues before / after its junction bond, so each parent gives the
  like-for-like baseline for "do these two segments pack together".
* **residue mapping** -- chimera residue -> parent residue through the
  parent-window alignment of :func:`cfh.domain_function.metrics.align_parents`
  and :func:`~cfh.domain_function.metrics.place_chimera`. Junction-inserted
  chimera residues and columns where the parent has a gap are unmapped.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Sequence

import numpy as np

from cfh.domain_function.metrics import (
    ChimeraColumns,
    ParentAlignment,
    align_map,
    align_parents,
    place_chimera,
)
from cfh.domain_function.record import ChimeraRecord, ParentDomain
from cfh.domain_structure.geometry import superposed_rmsd, tm_align, tm_score_fixed
from cfh.domain_structure.jobs import RecordJobs
from cfh.domain_structure.predictor import Prediction

JUNCTION_FLANK = 10
"""Residues on each side of the junction bond in the junction window."""
LOCAL_RADIUS = 8.0
"""CA-CA distance (A) from a junction residue that defines the local region."""
MIN_MAPPED_PAIRS = 10
"""Fewer mapped domain residues than this is a mapping failure."""
MIN_AFDB_IDENTITY = 0.9
"""Minimum transcript-window vs AFDB-sequence identity to map onto an AFDB model."""


def _r(value: float | None, digits: int = 2) -> float | None:
    return None if value is None else round(float(value), digits)


def _mean(values: Sequence[float]) -> float | None:
    return float(np.mean(values)) if len(values) else None


@dataclass(frozen=True)
class Window:
    """A prediction plus where its first residue sits in its own protein."""

    prediction: Prediction
    start: int

    @property
    def end(self) -> int:
        return self.start + len(self.prediction.sequence) - 1

    def inside(self, positions: Sequence[int]) -> list[int]:
        return [p for p in positions if self.start <= p <= self.end]

    def plddt(self, positions: Sequence[int]) -> list[float]:
        return [self.prediction.plddt[p - self.start] for p in self.inside(positions)]

    def mean_pae(self, first: Sequence[int], second: Sequence[int] | None = None) -> float | None:
        """Mean PAE over pairs (i != j); between two sets it is symmetrised."""
        a = [p - self.start for p in self.inside(first)]
        pae = np.asarray(self.prediction.pae, dtype=float)
        if second is None:
            if len(a) < 2:
                return None
            block = pae[np.ix_(a, a)]
            return float((block.sum() - np.trace(block)) / (len(a) * (len(a) - 1)))
        b = [p - self.start for p in self.inside(second)]
        if not a or not b:
            return None
        return float((pae[np.ix_(a, b)].mean() + pae[np.ix_(b, a)].mean()) / 2)


def junction_window(record: ChimeraRecord, inserted: Sequence[int]) -> list[int]:
    j = record.junction_index
    low, high = j - JUNCTION_FLANK + 1, j + JUNCTION_FLANK
    if inserted:
        low = min(low, min(inserted) - JUNCTION_FLANK)
        high = max(high, max(inserted) + JUNCTION_FLANK)
    return list(range(low, high + 1))


def parent_bond(parent: ParentDomain, which: str) -> int:
    """Last parent residue before the junction bond."""
    return parent.junction_residue if which == "five" else parent.junction_residue - 1


def parent_junction_window(parent: ParentDomain, which: str) -> list[int]:
    bond = parent_bond(parent, which)
    return list(range(bond - JUNCTION_FLANK + 1, bond + JUNCTION_FLANK + 1))


def chimera_halves(record: ChimeraRecord) -> tuple[list[int], list[int]]:
    five: list[int] = []
    three: list[int] = []
    for position in range(record.domain_start, record.domain_end + 1):
        which, _ = record.origin(position)
        if which == "five":
            five.append(position)
        elif which == "three":
            three.append(position)
    return five, three


def parent_halves(parent: ParentDomain, which: str) -> tuple[list[int], list[int]]:
    bond = parent_bond(parent, which)
    domain = range(parent.domain_start, parent.domain_end + 1)
    return [p for p in domain if p <= bond], [p for p in domain if p > bond]


def confidence(window: Window, domain: list[int], junction: list[int], halves) -> dict[str, Any]:
    domain_plddt = window.plddt(domain)
    junction_plddt = window.plddt(junction)
    return {
        "domain_plddt": _r(_mean(domain_plddt)),
        "junction_plddt_mean": _r(_mean(junction_plddt)),
        "junction_plddt_min": _r(min(junction_plddt) if junction_plddt else None),
        "domain_pae": _r(window.mean_pae(domain)),
        "interhalf_pae": _r(window.mean_pae(halves[0], halves[1])),
        "half_sizes": f"{len(window.inside(halves[0]))}/{len(window.inside(halves[1]))}",
        "ptm": _r(window.prediction.ptm, 3),
    }


def mapped_pairs(
    record: ChimeraRecord, parents: ParentAlignment, placed: ChimeraColumns, which: str
) -> list[tuple[int, int]]:
    """(chimera position, parent position) pairs with BOTH residues in their Pfam domains.

    The chimera side is restricted to the chimeric domain and the parent side
    to that parent's domain, so the pairs never outnumber the parent domain
    length that ``tm`` is normalised by (flank residues never count).
    """
    side = 0 if which == "five" else 1
    parent = record.five if which == "five" else record.three
    pairs = []
    for position in range(record.domain_start, record.domain_end + 1):
        column = placed.chimera_col.get(position)
        if column is None:
            continue
        parent_position = parents.columns[column][side]
        if (
            parent_position is not None
            and parent.domain_start <= parent_position <= parent.domain_end
        ):
            pairs.append((position, parent_position))
    return pairs


def _coords(window: Window) -> np.ndarray:
    return np.asarray(window.prediction.ca_coords(), dtype=float)


def compare(
    chimera: Window,
    target_coords: np.ndarray,
    target_start: int,
    pairs: list[tuple[int, int]],
    junction_residues: list[int],
    reference_length: int,
    chimera_length: int,
) -> dict[str, Any]:
    """Mapped TM-score/RMSD of the chimera domain and local RMSD around the junction.

    ``tm`` is normalised by ``reference_length`` (the parent domain), the
    usual convention for scoring a model against a reference: chimera residues
    with no parent counterpart (a junction insertion) do not lower it, parent
    domain residues missing from the chimera (a deletion) do. ``tm_chim_norm``
    is normalised by the chimera domain length instead.
    """
    pairs = [
        (c, p)
        for c, p in pairs
        if chimera.start <= c <= chimera.end and 0 <= p - target_start < len(target_coords)
    ]
    if len(pairs) < MIN_MAPPED_PAIRS:
        return {"status": f"mapping_failed: {len(pairs)} mapped residues"}
    if len(pairs) > min(reference_length, chimera_length):
        return {
            "status": f"mapping_failed: {len(pairs)} mapped residues exceed the domain length "
            f"({reference_length} reference, {chimera_length} chimera)"
        }
    if len({p for _, p in pairs}) != len(pairs) or len({c for c, _ in pairs}) != len(pairs):
        return {"status": "mapping_failed: residue mapping is not one-to-one"}
    ca = _coords(chimera)
    mobile = np.asarray([ca[c - chimera.start] for c, _ in pairs])
    target = np.asarray([target_coords[p - target_start] for _, p in pairs])
    tm = tm_score_fixed(mobile, target, reference_length)
    tm_chim = tm_score_fixed(mobile, target, chimera_length)
    if not (0.0 <= tm.tm_score <= 1.0 and 0.0 <= tm_chim.tm_score <= 1.0):
        return {"status": f"mapping_failed: TM-score out of range ({tm.tm_score})"}
    junction_ca = [ca[c - chimera.start] for c in chimera.inside(junction_residues)]
    local_rmsd = None
    n_local = 0
    if junction_ca:
        ref = np.asarray(junction_ca)
        near = [
            k
            for k, (c, _) in enumerate(pairs)
            if np.linalg.norm(ref - ca[c - chimera.start], axis=1).min() <= LOCAL_RADIUS
        ]
        n_local = len(near)
        if n_local >= 3:
            local_rmsd = superposed_rmsd(mobile[near], target[near])
    return {
        "status": "ok",
        "tm": tm.tm_score,
        "tm_chim_norm": tm_chim.tm_score,
        "rmsd": tm.rmsd,
        "n_mapped": tm.n_aligned,
        "local_rmsd": _r(local_rmsd, 3),
        "n_local": n_local,
    }


def _tmalign_domain(chimera: Window, record: ChimeraRecord, other: Window, parent: ParentDomain):
    c_dom = chimera.inside(list(range(record.domain_start, record.domain_end + 1)))
    p_dom = other.inside(list(range(parent.domain_start, parent.domain_end + 1)))
    if len(c_dom) < 3 or len(p_dom) < 3:
        return None
    ca, pa = _coords(chimera), _coords(other)
    result = tm_align(
        ca[[c - chimera.start for c in c_dom]],
        pa[[p - other.start for p in p_dom]],
        "".join(record.chimeric_protein[c - 1] for c in c_dom),
        "".join(parent.protein[p - 1] for p in p_dom),
    )
    return result.tm_norm_second


@dataclass(frozen=True)
class AfdbModel:
    accession: str
    version: str
    sequence: str
    ca: tuple[tuple[float, float, float], ...]
    plddt: tuple[float, ...]


def afdb_parent_map(parent: ParentDomain, model: AfdbModel) -> tuple[dict[int, int], float]:
    """Parent residue -> AFDB model index over the parent window, with window identity."""
    start, end = parent.window()
    window = parent.protein[start - 1 : end]
    pairs = align_map(window, model.sequence)
    mapping = {a + start: b for a, b in pairs if a is not None and b is not None}
    same = sum(1 for p, b in mapping.items() if parent.protein[p - 1] == model.sequence[b])
    return mapping, same / max(1, len(window))


def record_metrics(
    jobs: RecordJobs,
    predictions: dict[str, Prediction],
    failures: dict[str, str],
    afdb: dict[str, AfdbModel | str],
) -> dict[str, Any]:
    """Every reported measure for one record (missing pieces say why in ``*_status``)."""
    record = jobs.record
    row: dict[str, Any] = {}
    parents = align_parents(record)
    placed = place_chimera(record, parents)
    junction = junction_window(record, placed.inserted)
    j_residues = sorted({record.junction_index, record.junction_index + 1, *placed.inserted})
    domain = list(range(record.domain_start, record.domain_end + 1))
    row["junction_window"] = f"{junction[0]}-{junction[-1]}"
    row["n_inserted"] = len(placed.inserted)
    row["n_deleted_columns"] = len(placed.deleted_columns)

    windows: dict[str, Window | None] = {}
    for role, job_id, start in (
        ("chimera", jobs.chimera, jobs.chimera_start),
        ("parent5", jobs.parent5, jobs.parent5_start),
        ("parent3", jobs.parent3, jobs.parent3_start),
    ):
        prediction = predictions.get(job_id)
        row[f"{role}_job"] = job_id
        if prediction is None:
            row[f"{role}_prediction_status"] = failures.get(job_id, "prediction_missing")
            windows[role] = None
            continue
        row[f"{role}_prediction_status"] = "ok"
        row[f"{role}_wall_seconds"] = prediction.wall_seconds
        windows[role] = Window(prediction, start)

    chim = windows["chimera"]
    if chim is not None:
        for key, value in confidence(chim, domain, junction, chimera_halves(record)).items():
            row[f"chimera_{key}"] = value
    for role, parent, which in (
        ("parent5", record.five, "five"),
        ("parent3", record.three, "three"),
    ):
        window = windows[role]
        if window is None:
            continue
        p_domain = list(range(parent.domain_start, parent.domain_end + 1))
        values = confidence(
            window, p_domain, parent_junction_window(parent, which), parent_halves(parent, which)
        )
        for key, value in values.items():
            row[f"{role}_{key}"] = value

    for metric in ("domain_plddt", "junction_plddt_mean"):
        parent_values = [row.get(f"parent5_{metric}"), row.get(f"parent3_{metric}")]
        known = [v for v in parent_values if v is not None]
        chimera_value = row.get(f"chimera_{metric}")
        if chimera_value is not None and len(known) == 2:
            row[f"delta_{metric}_vs_better_parent"] = _r(chimera_value - max(known))

    # Chimera vs each parent prediction (like-for-like) and vs the AFDB parent models.
    chimera_length = record.domain_end - record.domain_start + 1
    for role, parent, which, short in (
        ("parent5", record.five, "five", "5p"),
        ("parent3", record.three, "three", "3p"),
    ):
        pairs = mapped_pairs(record, parents, placed, which)
        parent_length = parent.domain_end - parent.domain_start + 1
        window = windows[role]
        if chim is None or window is None:
            row[f"vs_{short}_status"] = "not_assessed: prediction missing"
        else:
            try:
                result = compare(
                    chim,
                    _coords(window),
                    window.start,
                    pairs,
                    j_residues,
                    parent_length,
                    chimera_length,
                )
                row[f"vs_{short}_status"] = result.pop("status")
                for key, value in result.items():
                    row[f"vs_{short}_{key}"] = value
                row[f"vs_{short}_tmalign_seqfree"] = _tmalign_domain(chim, record, window, parent)
            except (ValueError, np.linalg.LinAlgError) as error:
                row[f"vs_{short}_status"] = f"mapping_failed: {error}"
        model = afdb.get(f"{which}")
        if chim is None:
            row[f"vs_afdb_{short}_status"] = "not_assessed: prediction missing"
            continue
        if not isinstance(model, AfdbModel):
            row[f"vs_afdb_{short}_status"] = f"not_assessed: {model or 'no AFDB model'}"
            continue
        row[f"afdb_{short}_model"] = f"{model.accession} v{model.version}"
        mapping, identity = afdb_parent_map(parent, model)
        row[f"afdb_{short}_identity"] = _r(identity, 3)
        if identity < MIN_AFDB_IDENTITY:
            row[f"vs_afdb_{short}_status"] = "not_assessed: low identity to AFDB sequence"
            continue
        afdb_pairs = [(c, mapping[p]) for c, p in pairs if p in mapping]
        coords = np.asarray(model.ca, dtype=float)
        try:
            result = compare(chim, coords, 0, afdb_pairs, j_residues, parent_length, chimera_length)
        except (ValueError, np.linalg.LinAlgError) as error:
            result = {"status": f"mapping_failed: {error}"}
        row[f"vs_afdb_{short}_status"] = result.pop("status")
        for key, value in result.items():
            row[f"vs_afdb_{short}_{key}"] = value
        p_dom = [p for p in range(parent.domain_start, parent.domain_end + 1) if p in mapping]
        if window is not None and len(p_dom) >= MIN_MAPPED_PAIRS:
            # Sanity check of the settings: our parent-window model vs AFDB, parent domain.
            own = _coords(window)
            mobile = np.asarray([own[p - window.start] for p in p_dom])
            target = coords[[mapping[p] for p in p_dom]]
            check = tm_score_fixed(mobile, target, len(p_dom))
            row[f"parent{short[0]}_vs_afdb_tm"] = check.tm_score
            row[f"afdb_{short}_domain_plddt"] = _r(
                float(np.mean([model.plddt[mapping[p]] for p in p_dom]))
            )
    return row
