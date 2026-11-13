"""Standalone ESM-2 masked-marginal worker (run with the external PLM env's Python).

This file imports nothing from ``cfh``: it is executed as a script by the PLM
env (``$CFH_PLM_HOME/env/bin/python``, which has torch + fair-esm) and talks
JSON lines over stdin/stdout, so the ``cfh`` package and CI never need torch.

Protocol. On start the worker loads the model and writes one header line::

    {"header": {"model": ..., "esm_version": ..., "torch_version": ...,
                "weights_sha256": ..., "device": ...}}

Then for every request line ``{"id": str, "sequence": str, "positions": [int]}``
(1-based positions) it writes ``{"id": ..., "sequence": ..., "logp": {"<pos>":
float}}`` where each value is ``log p(true residue | sequence with only that
position masked)`` (natural log). Errors are written as ``{"id": ...,
"error": str}``; the worker never substitutes a value.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path
from typing import Any


def weights_sha256(model_name: str) -> str:
    import torch

    path = Path(torch.hub.get_dir()) / "checkpoints" / f"{model_name}.pt"
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    import esm
    import torch

    model_name = sys.argv[1]
    device = sys.argv[2] if len(sys.argv) > 2 else "mps"
    batch_size = int(sys.argv[3]) if len(sys.argv) > 3 else 64
    if device == "mps" and not torch.backends.mps.is_available():
        device = "cpu"
    model, alphabet = esm.pretrained.load_model_and_alphabet(model_name)
    model = model.eval().to(device)
    converter = alphabet.get_batch_converter()
    header = {
        "model": model_name,
        "esm_version": getattr(esm, "__version__", "unknown"),
        "torch_version": torch.__version__,
        "weights_sha256": weights_sha256(model_name),
        "device": device,
        "max_residues": 1022,
    }
    print(json.dumps({"header": header}), flush=True)

    for line in sys.stdin:
        if not line.strip():
            continue
        request: dict[str, Any] = json.loads(line)
        rid = request.get("id")
        try:
            sequence = request["sequence"]
            positions = sorted({int(p) for p in request["positions"]})
            if len(sequence) > header["max_residues"]:
                raise ValueError(f"sequence longer than {header['max_residues']} residues")
            if any(p < 1 or p > len(sequence) for p in positions):
                raise ValueError("position outside the sequence")
            _, _, tokens = converter([("q", sequence)])
            tokens = tokens.to(device)
            logp: dict[str, float] = {}
            with torch.no_grad():
                for start in range(0, len(positions), batch_size):
                    chunk = positions[start : start + batch_size]
                    batch = tokens.repeat(len(chunk), 1)
                    rows = torch.arange(len(chunk), device=device)
                    cols = torch.tensor(chunk, device=device)  # BOS at 0: residue p at p
                    batch[rows, cols] = alphabet.mask_idx
                    logits = model(batch)["logits"].float()
                    log_probs = torch.log_softmax(logits, dim=-1)
                    values = log_probs[rows, cols, tokens[0, cols]].cpu().tolist()
                    logp.update({str(p): float(v) for p, v in zip(chunk, values)})
            print(json.dumps({"id": rid, "sequence": sequence, "logp": logp}), flush=True)
        except Exception as error:  # reported per request, never replaced by a value
            print(json.dumps({"id": rid, "error": f"{type(error).__name__}: {error}"}), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
