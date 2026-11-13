"""Fetch MSAs for a batch of sequences from a ColabFold MMseqs2 server.

Runs inside the *ColabFold* environment (it imports ``colabfold``), not the cfh
environment: the driver executes this file by path with the ColabFold Python.
It uses only the standard library and ``colabfold.colabfold.run_mmseqs2``.

Input JSON (argv[1])::

    {"host_url": ..., "user_agent": ..., "use_env": true,
     "work_dir": ..., "queries": [{"key": ..., "sequence": ...}, ...]}

All queries go to the server as ONE ticket (``run_mmseqs2`` with a list), which
is how ColabFold itself batches; the server's RATELIMIT / queue answers are
handled by ``run_mmseqs2``'s own back-off. Output JSON (argv[2])::

    {"ok": true, "a3m": {key: "<a3m text>", ...}}  or  {"ok": false, "error": ...}

Nothing is written for a failed batch, so the caller never caches a failure.
"""

from __future__ import annotations

import json
import shutil
import sys
import traceback
from pathlib import Path


def main(request_path: str, result_path: str) -> int:
    request = json.loads(Path(request_path).read_text())
    work_dir = Path(request["work_dir"])
    work_dir.mkdir(parents=True, exist_ok=True)
    queries = request["queries"]
    try:
        from colabfold.colabfold import run_mmseqs2

        lines = run_mmseqs2(
            [q["sequence"] for q in queries],
            str(work_dir / "mmseqs2"),
            use_env=bool(request.get("use_env", True)),
            use_filter=True,
            use_templates=False,
            host_url=request["host_url"],
            user_agent=request["user_agent"],
        )
        if len(lines) != len(queries):
            raise RuntimeError(f"server returned {len(lines)} MSAs for {len(queries)} queries")
        a3m = {}
        for query, text in zip(queries, lines):
            first = [line for line in text.splitlines() if line and not line.startswith("#")]
            if len(first) < 2 or first[1].strip() != query["sequence"]:
                raise RuntimeError(f"MSA for {query['key']} does not start with its query")
            a3m[query["key"]] = text
        result = {"ok": True, "a3m": a3m}
    except BaseException as error:  # noqa: BLE001 - reported to the caller, never cached
        result = {"ok": False, "error": f"{type(error).__name__}: {error}"}
        traceback.print_exc()
    finally:
        shutil.rmtree(work_dir / "mmseqs2_env", ignore_errors=True)
        shutil.rmtree(work_dir / "mmseqs2_all", ignore_errors=True)
    Path(result_path).write_text(json.dumps(result))
    return 0 if result["ok"] else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2]))
