"""Parent structure features from AlphaFold DB models (one small download per accession).

For each UniProt accession the AlphaFold DB API gives the current model's
mmCIF. From it we keep only a compact derived JSON (the mmCIF is deleted):

* per-residue pLDDT (the CA B-factor column of AFDB models);
* secondary structure exactly as the file's ``_struct_conf`` category gives it
  (AFDB writes DSSP assignments there). A model without that category gives
  ``ss = None`` -- never an assumed loop;
* residue contacts: any heavy-atom pair within :data:`CONTACT_DISTANCE` A
  between residues at least :data:`MIN_SEQUENCE_SEPARATION` apart;
* burial: relative solvent accessibility from FreeSASA (Lee-Richards, default
  radii, relative to FreeSASA's tripeptide reference areas); a residue is buried
  core when RSA < :data:`BURIED_MAX_RSA` in the full-length parent model.
"""

from __future__ import annotations

import json
import time
from dataclasses import dataclass
from pathlib import Path

import freesasa
import gemmi
import requests

AFDB_API = "https://alphafold.ebi.ac.uk/api/prediction"
USER_AGENT = "cfh-domain-function/1 (research pipeline)"
CONTACT_DISTANCE = 4.5
MIN_SEQUENCE_SEPARATION = 3
BURIED_MAX_RSA = 0.20
"""A residue counts as buried core when its relative solvent accessibility is below this."""
CONFIDENT_PLDDT = 70.0

SS_CODES = {
    "HELX_RH_AL_P": "H",
    "HELX_RH_3T_P": "G",
    "HELX_RH_PI_P": "I",
    "STRN": "E",
    "HELX_LH_PP_P": "P",
    "TURN_TY1_P": "T",
    "BEND": "S",
}
HELIX_CODES = frozenset("HGI")
STRAND_CODES = frozenset("E")


@dataclass(frozen=True)
class ParentStructure:
    accession: str
    model_id: str
    model_version: str
    sequence: str
    """Full UniProt sequence the model was built on."""
    first: int
    """UniProt position of the model's first residue (AFDB fragments may start later)."""
    plddt: tuple[float, ...]
    ss: str | None
    """One DSSP-derived code per modelled residue ('-' = coil); None when not annotated."""
    contacts: tuple[tuple[int, int], ...]
    """UniProt positions (i < j)."""
    rsa: tuple[float | None, ...]
    """Relative solvent accessibility per modelled residue (FreeSASA, Lee-Richards)."""

    def index(self, position: int) -> int | None:
        offset = position - self.first
        return offset if 0 <= offset < len(self.plddt) else None

    def plddt_at(self, position: int) -> float | None:
        offset = self.index(position)
        return self.plddt[offset] if offset is not None else None

    def ss_at(self, position: int) -> str | None:
        offset = self.index(position)
        if offset is None or self.ss is None:
            return None
        return self.ss[offset]

    def buried(self, position: int) -> bool | None:
        offset = self.index(position)
        value = self.rsa[offset] if offset is not None else None
        return None if value is None else value < BURIED_MAX_RSA


def ss_class(code: str | None) -> str:
    """``helix`` / ``strand`` / ``loop`` from a DSSP code; ``not_assessed`` for None."""
    if code is None:
        return "not_assessed"
    if code in HELIX_CODES:
        return "helix"
    if code in STRAND_CODES:
        return "strand"
    return "loop"


def _is_amino_acid(name: str) -> bool:
    info = gemmi.find_tabulated_residue(name)
    return info is not None and info.is_amino_acid()


def _relative_sasa(
    structure: gemmi.Structure, path: Path, numbers: list[int]
) -> list[float | None]:
    pdb = path.with_suffix(".sasa.pdb")
    structure.write_pdb(str(pdb))
    try:
        sasa_structure = freesasa.Structure(str(pdb))
        areas = freesasa.calc(sasa_structure).residueAreas()
    finally:
        pdb.unlink(missing_ok=True)
    chain: dict = next(iter(areas.values()), {})
    values: list[float | None] = []
    for number in numbers:
        area = chain.get(str(number))
        relative = getattr(area, "relativeTotal", None) if area is not None else None
        values.append(None if relative is None else round(float(relative), 3))
    return values


def derive_from_mmcif(path: Path, accession: str, meta: dict) -> dict:
    """Compact derived features of one AFDB mmCIF (see module docstring)."""
    structure = gemmi.read_structure(str(path))
    structure.setup_entities()
    model = structure[0]
    chain = model[0]
    residues = [residue for residue in chain if _is_amino_acid(residue.name)]
    numbers = [int(residue.label_seq or 0) for residue in residues]
    if not numbers or numbers != list(range(numbers[0], numbers[0] + len(numbers))):
        raise ValueError(f"{accession}: missing or non-contiguous model numbering")
    label_first = numbers[0]
    first = int(meta.get("uniprotStart") or 1) + label_first - 1
    plddt = []
    for residue in residues:
        ca = residue.find_atom("CA", "*")
        if ca is None:
            raise ValueError(f"{accession}: residue {residue.label_seq} lacks CA")
        plddt.append(round(float(ca.b_iso), 2))
    rsa = _relative_sasa(structure, path, [int(residue.seqid.num or 0) for residue in residues])

    search = gemmi.NeighborSearch(model, structure.cell, 6).populate(include_h=False)
    contact_search = gemmi.ContactSearch(CONTACT_DISTANCE)
    contact_search.ignore = gemmi.ContactSearch.Ignore.SameResidue
    pairs: set[tuple[int, int]] = set()
    for result in contact_search.find_contacts(search):
        raw_a = result.partner1.residue.label_seq
        raw_b = result.partner2.residue.label_seq
        if raw_a is None or raw_b is None:
            continue
        a, b = int(raw_a), int(raw_b)
        if abs(a - b) < MIN_SEQUENCE_SEPARATION:
            continue
        i, j = sorted((a, b))
        pairs.add((i - label_first + first, j - label_first + first))

    document = gemmi.cif.read(str(path))
    block = document.sole_block()
    ss: str | None = None
    table = block.find("_struct_conf.", ["conf_type_id", "beg_label_seq_id", "end_label_seq_id"])
    if len(table) > 0:
        codes = ["-"] * len(residues)
        for row in table:
            code = SS_CODES.get(row[0], "-")
            for label in range(int(row[1]), int(row[2]) + 1):
                offset = label - label_first
                if 0 <= offset < len(codes):
                    codes[offset] = code
        ss = "".join(codes)
    return {
        "accession": accession,
        "model_id": meta.get("modelEntityId") or meta.get("entryId") or "",
        "model_version": str(meta.get("latestVersion") or ""),
        "sequence": meta.get("sequence") or meta.get("uniprotSequence") or "",
        "first": first,
        "plddt": plddt,
        "ss": ss,
        "contacts": sorted(pairs),
        "rsa": rsa,
    }


def structure_from_dict(data: dict) -> ParentStructure:
    return ParentStructure(
        accession=data["accession"],
        model_id=data["model_id"],
        model_version=data["model_version"],
        sequence=data["sequence"],
        first=int(data["first"]),
        plddt=tuple(float(x) for x in data["plddt"]),
        ss=data["ss"],
        contacts=tuple((int(a), int(b)) for a, b in data["contacts"]),
        rsa=tuple(None if x is None else float(x) for x in data["rsa"]),
    )


class AfdbClient:
    def __init__(self, cache_dir: Path, *, session: requests.Session | None = None):
        self.cache_dir = cache_dir
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.session = session or requests.Session()
        self.session.headers.setdefault("User-Agent", USER_AGENT)
        self.session.headers["User-Agent"] = USER_AGENT

    def _get(self, url: str) -> requests.Response:
        """GET with retries; 429/5xx back off (honouring Retry-After), 404 is returned."""
        last: Exception | None = None
        for attempt in range(6):
            try:
                response = self.session.get(url, timeout=120)
                if response.status_code == 404:
                    return response
                if response.status_code == 429 or response.status_code >= 500:
                    retry_after = response.headers.get("Retry-After", "")
                    wait = float(retry_after) if retry_after.isdigit() else 10.0 * 2**attempt
                    last = requests.HTTPError(f"{response.status_code} for {url}")
                    time.sleep(min(wait, 300.0))
                    continue
                response.raise_for_status()
                return response
            except requests.RequestException as error:
                last = error
                time.sleep(3 * (attempt + 1))
        raise RuntimeError(f"AFDB request failed for {url}: {last}")

    def structure(self, accession: str) -> ParentStructure | None:
        """Derived features for ``accession``'s AFDB model; None when AFDB has no model."""
        path = self.cache_dir / f"{accession}.json"
        if path.exists():
            data = json.loads(path.read_text())
            return structure_from_dict(data) if data.get("plddt") else None
        response = self._get(f"{AFDB_API}/{accession}")
        if response.status_code == 404:
            path.write_text(json.dumps({"accession": accession, "missing": "no AFDB model"}))
            return None
        entries = [e for e in response.json() if e.get("uniprotAccession", accession) == accession]
        entries = entries or response.json()
        meta = next(
            (e for e in entries if str(e.get("modelEntityId", "")).endswith("-F1")), entries[0]
        )
        cif = self._get(meta["cifUrl"])
        if cif.status_code == 404:
            path.write_text(json.dumps({"accession": accession, "missing": "mmCIF 404"}))
            return None
        tmp = self.cache_dir / f"{accession}.cif"
        tmp.write_bytes(cif.content)
        try:
            data = derive_from_mmcif(tmp, accession, meta)
        finally:
            tmp.unlink(missing_ok=True)
        path.write_text(json.dumps(data, separators=(",", ":")))
        return structure_from_dict(data)
