"""Structure-prediction backend: MSAs from a ColabFold MMseqs2 server, models from ColabFold.

Everything that decides *how* a structure is predicted lives in
:class:`PredictorConfig` (a JSON file can override any field), so moving to a
GPU host means pointing ``colabfold_command`` / ``msa_python`` / ``data_dir``
at that host's install -- no code changes. The settings that change a
prediction (model type, model count, recycles, MSA depth, templates, seed, ...)
are hashed into :meth:`PredictorConfig.settings_id`; cached predictions live
under that id, so changing a setting never reuses an old result.

Fail closed. A prediction is cached only after :func:`validate_outputs` has
checked it (exit status, exactly one rank-1 model and score file, pLDDT/PAE
shapes equal to the sequence length, finite pTM, model sequence equal to the
query). Anything else is a :class:`PredictionFailure` with a reason; the
failure is written next to the job for the record but is never read back as a
result, so a rerun tries again.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import time
from dataclasses import asdict, dataclass, field, fields
from pathlib import Path
from typing import Any, Callable

import gemmi

CACHE_ROOT = Path(os.environ.get("CFH_COLABFOLD_HOME", Path.home() / ".cache" / "cfh_colabfold"))
"""External (outside the git tree) home of the ColabFold env and weights."""

MSA_WORKER = Path(__file__).with_name("msa_worker.py")
USER_AGENT = "cfh-domain-structure/1"
"""Identifies the pipeline to the public MSA server (no personal contact details)."""


class DiskSpaceError(RuntimeError):
    """Free disk fell below ``min_free_disk_gb``: the run stops cleanly."""


class PredictionFailure(RuntimeError):
    """A prediction or its MSA could not be produced or validated."""


@dataclass
class PredictorConfig:
    backend: str = "colabfold_batch"
    """Only ``colabfold_batch`` (ColabFold's CLI) is implemented."""
    colabfold_command: list[str] = field(
        default_factory=lambda: [str(CACHE_ROOT / "env" / "bin" / "colabfold_batch")]
    )
    msa_python: str = str(CACHE_ROOT / "env" / "bin" / "python")
    """Python of the ColabFold env, used to run :mod:`.msa_worker`."""
    data_dir: str = str(CACHE_ROOT / "weights")
    """ColabFold ``--data``: directory holding ``params/params_model_*_ptm.npz``."""
    model_type: str = "alphafold2_ptm"
    num_models: int = 1
    model_order: str = "1"
    num_recycle: int = 1
    recycle_early_stop_tolerance: float | None = None
    templates: bool = False
    random_seed: int = 0
    num_seeds: int = 1
    max_msa: str = "256:512"
    """``max_seq:max_extra_seq``; ColabFold's default for alphafold2_ptm is 512:5120."""
    msa_mode: str = "mmseqs2_uniref_env"
    host_url: str = "https://api.colabfold.com"
    user_agent: str = USER_AGENT
    msa_batch_size: int = 8
    """Sequences per MSA-server ticket."""
    msa_pause_seconds: float = 20.0
    """Pause between MSA tickets (on top of the server's own RATELIMIT back-off)."""
    job_timeout_seconds: int = 6 * 3600
    min_free_disk_gb: float = 5.0
    extra_args: list[str] = field(default_factory=list)

    PREDICTION_SETTINGS = (
        "backend",
        "model_type",
        "num_models",
        "model_order",
        "num_recycle",
        "recycle_early_stop_tolerance",
        "templates",
        "random_seed",
        "num_seeds",
        "max_msa",
        "msa_mode",
        "host_url",
        "extra_args",
    )
    """Everything that changes a prediction, including where its MSA came from."""

    @classmethod
    def load(cls, path: Path | None) -> PredictorConfig:
        if path is None:
            return cls()
        data = json.loads(path.read_text())
        known = {f.name for f in fields(cls)}
        unknown = set(data) - known
        if unknown:
            raise ValueError(f"unknown predictor settings: {sorted(unknown)}")
        return cls(**data)

    def prediction_settings(self) -> dict[str, Any]:
        data = asdict(self)
        return {key: data[key] for key in self.PREDICTION_SETTINGS}

    def settings_id(self) -> str:
        text = json.dumps(self.prediction_settings(), sort_keys=True)
        return "set_" + hashlib.sha256(text.encode()).hexdigest()[:12]

    def msa_id(self) -> str:
        return hashlib.sha256(f"{self.msa_mode}|{self.host_url}".encode()).hexdigest()[:8]

    def msa_source(self) -> dict[str, str]:
        """MSA provenance written into every checkpoint and required when it is loaded."""
        return {
            "msa_mode": self.msa_mode,
            "host_url": self.host_url,
            "msa_cache": f"msa_{self.msa_id()}",
        }

    def batch_args(self) -> list[str]:
        if self.backend != "colabfold_batch":
            raise ValueError(f"unsupported backend {self.backend!r}")
        args = [
            "--model-type",
            self.model_type,
            "--num-models",
            str(self.num_models),
            "--model-order",
            self.model_order,
            "--num-recycle",
            str(self.num_recycle),
            "--random-seed",
            str(self.random_seed),
            "--num-seeds",
            str(self.num_seeds),
            "--max-msa",
            self.max_msa,
            "--msa-mode",
            self.msa_mode,
            "--data",
            self.data_dir,
        ]
        if self.recycle_early_stop_tolerance is not None:
            args += ["--recycle-early-stop-tolerance", str(self.recycle_early_stop_tolerance)]
        if self.templates:
            args.append("--templates")
        return args + list(self.extra_args)


def sequence_hash(sequence: str) -> str:
    return hashlib.sha256(sequence.encode()).hexdigest()[:16]


def check_disk(path: Path, min_free_gb: float) -> float:
    free = shutil.disk_usage(path).free / 1e9
    if free < min_free_gb:
        raise DiskSpaceError(f"only {free:.1f} GB free at {path} (minimum {min_free_gb} GB)")
    return free


# --------------------------------------------------------------------------- MSAs
class MsaCache:
    """a3m files per sequence; a cached a3m is never re-queried."""

    def __init__(self, root: Path, config: PredictorConfig):
        self.config = config
        self.root = root / f"msa_{config.msa_id()}"
        self.root.mkdir(parents=True, exist_ok=True)

    def path(self, sequence: str) -> Path:
        return self.root / f"{sequence_hash(sequence)}.a3m"

    def get(self, sequence: str) -> str | None:
        path = self.path(sequence)
        if not path.exists():
            return None
        text = path.read_text()
        return text if a3m_query(text) == sequence else None

    def fetch(
        self,
        sequences: list[str],
        log: Callable[[str], None] = print,
        runner: Callable[[list[str]], int] | None = None,
    ) -> dict[str, str]:
        """Query the server for every uncached sequence, in batches. Returns errors by sequence."""
        missing = sorted({s for s in sequences if self.get(s) is None}, key=len)
        errors: dict[str, str] = {}
        batch_size = max(1, self.config.msa_batch_size)
        for index in range(0, len(missing), batch_size):
            check_disk(self.root, self.config.min_free_disk_gb)
            batch = missing[index : index + batch_size]
            work = self.root / "_ticket"
            shutil.rmtree(work, ignore_errors=True)
            work.mkdir()
            request = {
                "host_url": self.config.host_url,
                "user_agent": self.config.user_agent,
                "use_env": self.config.msa_mode == "mmseqs2_uniref_env",
                "work_dir": str(work.resolve()),
                "queries": [{"key": sequence_hash(s), "sequence": s} for s in batch],
            }
            (work / "request.json").write_text(json.dumps(request))
            command = [
                self.config.msa_python,
                "-I",
                str(MSA_WORKER),
                str(work / "request.json"),
                str(work / "result.json"),
            ]
            log(f"MSA ticket {index // batch_size + 1}: {len(batch)} sequences")
            with open(self.root / "msa_worker.log", "a") as handle:
                code = (runner or _run_logged(handle))(command)
            result_path = work / "result.json"
            result = json.loads(result_path.read_text()) if result_path.exists() else {}
            if code != 0 or not result.get("ok"):
                reason = f"msa_server_error: {result.get('error') or f'worker exit {code}'}"
                for sequence in batch:
                    errors[sequence] = reason
                log(reason)
            else:
                for sequence in batch:
                    text = result["a3m"].get(sequence_hash(sequence), "")
                    if a3m_query(text) != sequence:
                        errors[sequence] = "msa_server_error: MSA query differs from sequence"
                        continue
                    tmp = self.path(sequence).with_suffix(".tmp")
                    tmp.write_text(text)
                    tmp.replace(self.path(sequence))
            shutil.rmtree(work, ignore_errors=True)
            if index + batch_size < len(missing):
                time.sleep(self.config.msa_pause_seconds)
        return errors


def _run_logged(handle: Any) -> Callable[[list[str]], int]:
    def run(command: list[str]) -> int:
        return subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT).returncode

    return run


def a3m_query(text: str) -> str | None:
    """First sequence of an a3m (the query), or None."""
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    lines = [line for line in lines if not line.startswith("#")]
    if len(lines) < 2 or not lines[0].startswith(">"):
        return None
    return lines[1]


# --------------------------------------------------------------------------- predictions
@dataclass
class Prediction:
    """One validated top-ranked model and its confidence scores."""

    job_id: str
    sequence: str
    plddt: list[float]
    pae: list[list[float]]
    ptm: float
    model_name: str
    pdb_path: Path
    wall_seconds: float | None = None
    settings_id: str = ""
    msa_source: dict[str, str] = field(default_factory=dict)
    """Recorded at prediction time: MSA mode, server and cache the a3m came from."""
    a3m_sha256: str = ""
    msa_provenance: str = ""
    """``recorded`` (written at prediction time) or ``inferred`` (written afterwards by
    :mod:`cfh.domain_structure.provenance` from documented run-level evidence)."""
    """sha256 of the exact a3m text handed to ColabFold, recorded at prediction time."""

    def ca_coords(self) -> list[tuple[float, float, float]]:
        return read_ca(self.pdb_path, self.sequence)


THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}  # fmt: skip


def read_ca(path: Path, sequence: str) -> list[tuple[float, float, float]]:
    """CA coordinates of a single-chain model; raises unless it is exactly ``sequence``."""
    if path.suffix == ".gz":
        with gzip.open(path, "rt") as handle:
            structure = gemmi.read_pdb_string(handle.read())
    else:
        structure = gemmi.read_structure(str(path))
    if len(structure) == 0 or len(structure[0]) != 1:
        raise PredictionFailure(f"parse_error: {path.name} is not a single-chain model")
    chain = structure[0][0]
    letters: list[str] = []
    coords: list[tuple[float, float, float]] = []
    for residue in chain:
        atom = residue.find_atom("CA", "*")
        if atom is None:
            raise PredictionFailure(f"parse_error: residue {residue.seqid.num} lacks CA")
        letters.append(THREE_TO_ONE.get(residue.name, "X"))
        coords.append((atom.pos.x, atom.pos.y, atom.pos.z))
    if "".join(letters) != sequence:
        raise PredictionFailure("parse_error: model sequence differs from the query")
    return coords


def validate_outputs(out_dir: Path, job_id: str, sequence: str) -> tuple[Path, Path, dict]:
    """(model pdb, score json, parsed scores) of the rank-1 model; raises PredictionFailure."""
    scores = sorted(out_dir.glob(f"{job_id}_scores_rank_001_*.json"))
    models = sorted(out_dir.glob(f"{job_id}_unrelaxed_rank_001_*.pdb"))
    if len(scores) != 1 or len(models) != 1:
        raise PredictionFailure(
            f"missing_model_file: {len(models)} rank-1 models, {len(scores)} score files"
        )
    try:
        data = json.loads(scores[0].read_text())
    except (OSError, ValueError) as error:
        raise PredictionFailure(f"parse_error: score json: {error}") from error
    check_scores(data, len(sequence))
    read_ca(models[0], sequence)
    return models[0], scores[0], data


def _finite_number(value: Any) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def check_scores(data: Any, length: int) -> None:
    """Raise PredictionFailure unless ``data`` holds finite pLDDT[L], PAE[L][L] and pTM.

    Applied identically to fresh ColabFold output and to every loaded checkpoint.
    """
    if not isinstance(data, dict):
        raise PredictionFailure("parse_error: scores are not a JSON object")
    plddt, pae, ptm = data.get("plddt"), data.get("pae"), data.get("ptm")
    if not isinstance(plddt, list) or len(plddt) != length:
        raise PredictionFailure("parse_error: pLDDT length differs from the sequence")
    if (
        not isinstance(pae, list)
        or len(pae) != length
        or any(not isinstance(row, list) or len(row) != length for row in pae)
    ):
        raise PredictionFailure("parse_error: PAE is not L x L")
    values = [ptm, *plddt, *(v for row in pae for v in row)]
    if not all(_finite_number(v) for v in values):
        raise PredictionFailure("parse_error: non-finite or non-numeric pLDDT/PAE/pTM")


class Predictor:
    """Runs/caches one ColabFold prediction per sequence under ``root/<settings_id>/``."""

    def __init__(self, root: Path, config: PredictorConfig, msas: MsaCache):
        self.config = config
        self.msas = msas
        self.root = root / config.settings_id()
        self.root.mkdir(parents=True, exist_ok=True)

    def job_dir(self, job_id: str) -> Path:
        return self.root / job_id

    def cached(self, job_id: str, sequence: str) -> Prediction | None:
        """A valid cached prediction, or None (a recorded failure is not a result)."""
        directory = self.job_dir(job_id)
        result = directory / "result.json"
        if not result.exists():
            return None
        try:
            return load_prediction(directory, sequence, self.config.msa_source())
        except PredictionFailure as error:
            if not str(error).startswith("provenance"):
                result.unlink(missing_ok=True)
            return None  # a provenance mismatch is kept on disk but never used
        except (OSError, ValueError, KeyError):
            result.unlink(missing_ok=True)
            return None

    def predict(self, job_id: str, sequence: str, log: Callable[[str], None] = print) -> Prediction:
        cached = self.cached(job_id, sequence)
        if cached is not None:
            return cached
        directory = self.job_dir(job_id)
        check_disk(self.root, self.config.min_free_disk_gb)
        a3m = self.msas.get(sequence)
        if a3m is None:
            raise self._fail(directory, "msa_missing: no cached MSA for this sequence")
        work = directory / "colabfold"
        shutil.rmtree(work, ignore_errors=True)
        (work / "in").mkdir(parents=True)
        (work / "in" / f"{job_id}.a3m").write_text(a3m)
        command = [
            *self.config.colabfold_command,
            str(work / "in"),
            str(work / "out"),
            *self.config.batch_args(),
        ]
        started = time.time()
        log(f"predict {job_id} ({len(sequence)} aa)")
        try:
            with open(directory / "colabfold.log", "w") as handle:
                handle.write(" ".join(command) + "\n")
                handle.flush()
                code = subprocess.run(
                    command,
                    stdout=handle,
                    stderr=subprocess.STDOUT,
                    timeout=self.config.job_timeout_seconds,
                ).returncode
        except subprocess.TimeoutExpired as error:
            raise self._fail(directory, f"prediction_timeout: > {error.timeout} s") from error
        except OSError as error:
            raise self._fail(directory, f"prediction_failed: {error}") from error
        wall = round(time.time() - started, 1)
        if code != 0:
            raise self._fail(directory, f"prediction_failed: colabfold exit status {code}")
        try:
            model, _, scores = validate_outputs(work / "out", job_id, sequence)
        except PredictionFailure as error:
            raise self._fail(directory, str(error)) from error
        with open(model, "rb") as src, gzip.open(directory / "model.pdb.gz", "wb") as dst:
            shutil.copyfileobj(src, dst)
        with gzip.open(directory / "scores.json.gz", "wt") as handle:
            json.dump(
                {"plddt": scores["plddt"], "pae": scores["pae"], "ptm": scores["ptm"]}, handle
            )
        meta = {
            "job_id": job_id,
            "sequence": sequence,
            "model_name": model.name.split("_rank_001_")[-1].removesuffix(".pdb"),
            "wall_seconds": wall,
            "settings_id": self.config.settings_id(),
            "msa_source": self.config.msa_source(),
            "msa_provenance": "recorded",
            "a3m_sha256": hashlib.sha256(a3m.encode()).hexdigest(),
        }
        shutil.rmtree(work, ignore_errors=True)
        (directory / "failure.json").unlink(missing_ok=True)
        tmp = directory / "result.json.tmp"
        tmp.write_text(json.dumps(meta, indent=1))
        tmp.replace(directory / "result.json")
        prediction = load_prediction(directory, sequence, self.config.msa_source())
        log(f"done {job_id} in {wall:.0f} s")
        return prediction

    def _fail(self, directory: Path, reason: str) -> PredictionFailure:
        directory.mkdir(parents=True, exist_ok=True)
        (directory / "result.json").unlink(missing_ok=True)
        (directory / "failure.json").write_text(json.dumps({"reason": reason, "at": time.time()}))
        return PredictionFailure(reason)


def load_prediction(
    directory: Path, sequence: str, expected_msa_source: dict[str, str] | None = None
) -> Prediction:
    """Load and re-validate a stored prediction directory (result.json + model + scores).

    With ``expected_msa_source`` the checkpoint must carry an MSA provenance
    record -- the same ``msa_source``, an a3m sha256 and ``msa_provenance``
    ``recorded`` (written at prediction time) or ``inferred`` (written later by
    :func:`cfh.domain_structure.provenance.record_inferred_provenance`, with its
    evidence) -- else a ``provenance: ...`` PredictionFailure is raised. An
    inferred record is accepted only for the exact source it names; a
    missing or different record is never a result.
    """
    meta = json.loads((directory / "result.json").read_text())
    if meta.get("sequence") != sequence:
        raise PredictionFailure("parse_error: cached prediction is for another sequence")
    if expected_msa_source is not None:
        recorded = meta.get("msa_source")
        if recorded != expected_msa_source:
            raise PredictionFailure(
                f"provenance: checkpoint MSA source {recorded!r} != {expected_msa_source!r}"
            )
        digest = meta.get("a3m_sha256")
        if not (isinstance(digest, str) and re.fullmatch(r"[0-9a-f]{64}", digest)):
            raise PredictionFailure("provenance: checkpoint has no a3m sha256")
        kind = meta.get("msa_provenance")
        if kind not in ("recorded", "inferred"):
            raise PredictionFailure(f"provenance: unknown msa_provenance {kind!r}")
        if kind == "inferred" and not meta.get("msa_provenance_evidence"):
            raise PredictionFailure("provenance: inferred record without evidence")
    scores_path = directory / "scores.json.gz"
    with gzip.open(scores_path, "rt") as handle:
        scores = json.load(handle)
    check_scores(scores, len(sequence))
    model = directory / "model.pdb.gz"
    read_ca(model, sequence)
    return Prediction(
        job_id=meta["job_id"],
        sequence=sequence,
        plddt=[float(x) for x in scores["plddt"]],
        pae=[[float(x) for x in row] for row in scores["pae"]],
        ptm=float(scores["ptm"]),
        model_name=meta.get("model_name", ""),
        pdb_path=model,
        wall_seconds=meta.get("wall_seconds"),
        settings_id=meta.get("settings_id", ""),
        msa_source=dict(meta.get("msa_source") or {}),
        a3m_sha256=meta.get("a3m_sha256", ""),
        msa_provenance=meta.get("msa_provenance", ""),
    )


def colabfold_version(config: PredictorConfig) -> str:
    try:
        out = subprocess.run(
            [config.msa_python, "-I", "-c", "import importlib.metadata as m;"
             "print(m.version('colabfold'), m.version('alphafold-colabfold'), m.version('jax'))"],
            capture_output=True, text=True, timeout=120,
        )  # fmt: skip
        parts = out.stdout.split()
        if len(parts) == 3:
            return f"colabfold {parts[0]}; alphafold-colabfold {parts[1]}; jax {parts[2]}"
    except (OSError, subprocess.SubprocessError):
        pass
    return "unknown"


def weights_manifest(config: PredictorConfig) -> dict[str, str]:
    """sha256 of the weight files the configured models use (params_model_<n>_ptm.npz)."""
    params = Path(config.data_dir) / "params"
    out: dict[str, str] = {}
    suffix = "_ptm" if config.model_type == "alphafold2_ptm" else ""
    for number in config.model_order.split(",")[: config.num_models]:
        path = params / f"params_model_{number.strip()}{suffix}.npz"
        if path.exists():
            digest = hashlib.sha256()
            with open(path, "rb") as handle:
                for chunk in iter(lambda: handle.read(1 << 20), b""):
                    digest.update(chunk)
            out[path.name] = digest.hexdigest()
        else:
            out[path.name] = "missing"
    return out
