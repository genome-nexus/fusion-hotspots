"""Committed outputs: metrics TSV, junction profiles, counts, report tables and manifest."""

from __future__ import annotations

from pathlib import Path
from typing import Any


def write_outputs(
    out: Path,
    rows: list[dict[str, Any]],
    profiles: list[tuple[str, str, int, str, float]],
    manifest: dict[str, Any],
    named: tuple[str, ...],
) -> None:
    raise NotImplementedError
