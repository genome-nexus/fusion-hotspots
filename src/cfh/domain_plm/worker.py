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

MAX_RESIDUES = 1022
"""ESM-2 positional limit (1024 tokens including BOS/EOS)."""


def context_windows(
    length: int, positions: list[int], limit: int = MAX_RESIDUES
) -> dict[int, list[int]]:
    """Group 1-based positions by the 1-based start of the context they are scored in.

    A sequence of at most ``limit`` residues is scored whole (one group, start 1).
    A longer one is scored per position in the ``limit``-residue window centred
    on that position, shifted inwards at the ends; positions sharing a window
    share a forward batch.
    """
    if length <= limit:
        return {1: sorted(positions)}
    groups: dict[int, list[int]] = {}
    for position in sorted(positions):
        start = min(max(1, position - limit // 2), length - limit + 1)
        groups.setdefault(start, []).append(position)
    return groups


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
        "max_residues": MAX_RESIDUES,
        "long_sequences": "per-position centred window of max_residues residues",
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
            if any(p < 1 or p > len(sequence) for p in positions):
                raise ValueError("position outside the sequence")
            logp: dict[str, float] = {}
            for offset, group in context_windows(len(sequence), positions).items():
                context = sequence[offset - 1 : offset - 1 + MAX_RESIDUES]
                _, _, tokens = converter([("q", context)])
                tokens = tokens.to(device)
                with torch.no_grad():
                    for start in range(0, len(group), batch_size):
                        chunk = group[start : start + batch_size]
                        batch = tokens.repeat(len(chunk), 1)
                        rows = torch.arange(len(chunk), device=device)
                        # BOS at token 0: residue p of the context is token p
                        cols = torch.tensor([p - offset + 1 for p in chunk], device=device)
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
