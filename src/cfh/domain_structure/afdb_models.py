"""CA coordinates of the parents' AlphaFold DB models (the models PR 1 used).

PR 1 (:class:`cfh.domain_function.afdb.AfdbClient`) keeps only derived
features, not coordinates. Here the same entry (current version, fragment F1)
is fetched again and only sequence, CA coordinates and pLDDT are cached as a
compact JSON; the mmCIF is deleted. Any failure is recorded as a reason
string -- never a model.
"""

from __future__ import annotations

import json
from pathlib import Path

import gemmi

from cfh.domain_function.afdb import AFDB_API, AfdbClient
from cfh.domain_structure.metrics import AfdbModel
from cfh.domain_structure.predictor import THREE_TO_ONE


def derive_ca(path: Path, accession: str, version: str, sequence: str) -> dict:
    structure = gemmi.read_structure(str(path))
    chain = structure[0][0]
    ca, plddt, letters = [], [], []
    for residue in chain:
        atom = residue.find_atom("CA", "*")
        if atom is None:
            raise ValueError(f"{accession}: residue {residue.seqid.num} lacks CA")
        letters.append(THREE_TO_ONE.get(residue.name, "X"))
        ca.append([round(atom.pos.x, 3), round(atom.pos.y, 3), round(atom.pos.z, 3)])
        plddt.append(round(float(atom.b_iso), 2))
    model_sequence = "".join(letters)
    if sequence and model_sequence != sequence:
        raise ValueError(f"{accession}: AFDB model sequence differs from the entry sequence")
    return {
        "accession": accession,
        "version": version,
        "sequence": model_sequence,
        "ca": ca,
        "plddt": plddt,
    }


def model_from_dict(data: dict) -> AfdbModel:
    return AfdbModel(
        accession=data["accession"],
        version=str(data["version"]),
        sequence=data["sequence"],
        ca=tuple((float(a), float(b), float(c)) for a, b, c in data["ca"]),
        plddt=tuple(float(x) for x in data["plddt"]),
    )


class AfdbCoordinates(AfdbClient):
    """Same API access and retry policy as PR 1's client; caches CA-only JSON."""

    def model(self, accession: str) -> AfdbModel | str:
        """The model, or a reason string (``no AFDB model``, ``afdb_fetch_failed: ...``)."""
        if not accession:
            return "no UniProt accession"
        path = self.cache_dir / f"{accession}.ca.json"
        if path.exists():
            data = json.loads(path.read_text())
            return model_from_dict(data) if data.get("ca") else str(data.get("missing"))
        try:
            response = self._get(f"{AFDB_API}/{accession}")
            if response.status_code == 404:
                path.write_text(json.dumps({"accession": accession, "missing": "no AFDB model"}))
                return "no AFDB model"
            entries = [
                e for e in response.json() if e.get("uniprotAccession", accession) == accession
            ] or response.json()
            meta = next(
                (e for e in entries if str(e.get("modelEntityId", "")).endswith("-F1")),
                entries[0],
            )
            cif = self._get(meta["cifUrl"])
            if cif.status_code == 404:
                return "afdb_fetch_failed: mmCIF 404"
            tmp = self.cache_dir / f"{accession}.cif"
            tmp.write_bytes(cif.content)
            try:
                data = derive_ca(
                    tmp,
                    accession,
                    str(meta.get("latestVersion") or ""),
                    meta.get("sequence") or meta.get("uniprotSequence") or "",
                )
            finally:
                tmp.unlink(missing_ok=True)
        except (RuntimeError, ValueError, KeyError, IndexError) as error:
            return f"afdb_fetch_failed: {error}"
        path.write_text(json.dumps(data, separators=(",", ":")))
        return model_from_dict(data)
