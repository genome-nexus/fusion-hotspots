"""Shared helpers for rendering configured biological domain names."""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping
from typing import Any


def domain_label_for_accession(
    domains: Iterable[Mapping[str, Any]], accession: str | None
) -> str | None:
    """Resolve a source accession to the human-readable name in ``domains``.

    The accession itself is retained as a graceful fallback for legacy payloads
    that do not carry a configured name.
    """
    if not accession:
        return None
    for domain in domains:
        if domain.get("accession") == accession:
            return domain.get("name") or accession
    return accession


def configured_domain_names(domains: Iterable[Mapping[str, Any]]) -> list[str]:
    """Return configured domain names, de-duplicated without changing order."""
    names: list[str] = []
    for domain in domains:
        name = domain.get("name")
        accession = domain.get("accession")
        if isinstance(name, str) and name and name != accession and name not in names:
            names.append(name)
    return names


def format_domain_names(names: Iterable[str]) -> str | None:
    """Join one or more domain names as a natural-language noun phrase."""
    unique = list(dict.fromkeys(name for name in names if name))
    if not unique:
        return None
    if len(unique) == 1:
        return unique[0]
    if len(unique) == 2:
        return f"{unique[0]} and {unique[1]}"
    return f"{', '.join(unique[:-1])}, and {unique[-1]}"


def domain_interpretation_sentence(
    domains: Iterable[Mapping[str, Any]],
    *,
    fisher_p_value: Any,
    fisher_odds_ratio: Any,
    effect: str,
    alpha: float = 0.05,
) -> str | None:
    """Render a configured-domain conclusion only for a supported Fisher effect."""
    if effect not in {"retention", "disruption"}:
        raise ValueError(f"unsupported domain interpretation effect: {effect}")
    if isinstance(fisher_p_value, bool) or isinstance(fisher_odds_ratio, bool):
        return None
    if not isinstance(fisher_p_value, (int, float)) or not isinstance(
        fisher_odds_ratio, (int, float)
    ):
        return None
    if (
        not math.isfinite(fisher_p_value)
        or not math.isfinite(fisher_odds_ratio)
        or fisher_p_value >= alpha
        or fisher_odds_ratio <= 1
    ):
        return None
    names = configured_domain_names(domains)
    domain = format_domain_names(names)
    if not domain:
        return None
    verb = "appears" if len(names) == 1 else "appear"
    if effect == "retention":
        return f"The {domain} {verb} to be required for retention."
    return f"The {domain} {verb} to require loss or disruption rather than retention."
