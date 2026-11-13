"""Per-sequence PLM score cache, keyed by model spec + sequence hash; fail closed on read.

Layout: ``<root>/<spec_id>/<sha[:2]>/<sha>.json`` where ``spec_id`` hashes the
model name, fair-esm version and the sha256 of the weights file, and ``sha`` is
the sha256 of the sequence. Each entry repeats the full spec and the sequence;
:meth:`ScoreCache.get` re-checks them and every value, so an entry written by a
different model, for a different sequence, or holding a NaN/inf/positive
log-probability is treated as absent (stale), never as a score.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class ModelSpec:
    model: str
    esm_version: str
    weights_sha256: str

    def spec_id(self) -> str:
        text = json.dumps(asdict(self), sort_keys=True)
        return f"plm_{hashlib.sha256(text.encode()).hexdigest()[:16]}"

    @classmethod
    def from_header(cls, header: dict[str, Any]) -> "ModelSpec":
        return cls(
            str(header["model"]), str(header["esm_version"]), str(header["weights_sha256"])
        )


def sequence_sha(sequence: str) -> str:
    return hashlib.sha256(sequence.encode()).hexdigest()


def valid_logp(value: Any) -> bool:
    """A log-probability: a finite real number <= 0 (bools are rejected)."""
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return False
    return math.isfinite(value) and value <= 0.0


class ScoreCache:
    def __init__(self, root: Path, spec: ModelSpec):
        self.root = Path(root)
        self.spec = spec
        self.dir = self.root / spec.spec_id()

    def path(self, sequence: str) -> Path:
        sha = sequence_sha(sequence)
        return self.dir / sha[:2] / f"{sha}.json"

    def get(self, sequence: str) -> dict[int, float]:
        """Validated cached scores for ``sequence`` (position -> logp); {} when absent/stale."""
        path = self.path(sequence)
        if not path.exists():
            return {}
        try:
            entry = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError):
            return {}
        if (
            entry.get("spec") != asdict(self.spec)
            or entry.get("sequence") != sequence
            or entry.get("sequence_sha256") != sequence_sha(sequence)
        ):
            return {}
        out: dict[int, float] = {}
        for key, value in (entry.get("logp") or {}).items():
            try:
                position = int(key)
            except (TypeError, ValueError):
                return {}
            if not 1 <= position <= len(sequence) or not valid_logp(value):
                return {}
            out[position] = float(value)
        return out

    def put(self, sequence: str, logp: dict[int, float]) -> None:
        """Merge validated scores into the entry (atomic write). Invalid values raise."""
        for position, value in logp.items():
            if not 1 <= position <= len(sequence) or not valid_logp(value):
                raise ValueError(f"invalid score at position {position}: {value!r}")
        merged = {**self.get(sequence), **logp}
        path = self.path(sequence)
        path.parent.mkdir(parents=True, exist_ok=True)
        entry = {
            "spec": asdict(self.spec),
            "sequence": sequence,
            "sequence_sha256": sequence_sha(sequence),
            "logp": {str(p): merged[p] for p in sorted(merged)},
        }
        tmp = path.with_suffix(".tmp")
        tmp.write_text(json.dumps(entry))
        os.replace(tmp, path)
