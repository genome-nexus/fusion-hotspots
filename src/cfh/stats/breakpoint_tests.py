"""Tests for frame-specific retention/disruption of configured protein domains.

The inputs deliberately remain the pipeline's typed ``FusionEvent`` and
``FusionFeature`` objects.  This keeps the statistical layer independent of a
particular cBioPortal row format while retaining the event/feature join that
was established during normalization and mapping.

The same machinery answers two symmetric questions, distinguished only by
which domain list and which target statuses ("hit statuses") are supplied:
whether a configured domain is reliably *retained* among in-frame fusions
(the ``domain_retention`` algorithm), or reliably *lost/disrupted* among
them (the ``domain_disruption`` algorithm).
"""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np
from scipy.stats import fisher_exact

from cfh.genes.registry import GeneConfig, KeyDomain
from cfh.mapping.genome_nexus_source import (
    GenomeNexusClient,
    cds_bounds_from_utrs,
    map_genomic_breakpoint_to_protein_position,
    parse_canonical_transcript,
)
from cfh.model.fusion_event import FusionEvent
from cfh.model.fusion_feature import FusionFeature

_RETAINED = "retained"
_NON_RETAINED = {"lost", "disrupted"}
_KNOWN_STATUSES = {_RETAINED, *_NON_RETAINED}

RETAINED_STATUSES = frozenset({_RETAINED})
"""Hit-status set for testing domain *retention* (the default)."""

DISRUPTED_STATUSES = frozenset(_NON_RETAINED)
"""Hit-status set for testing domain *loss/disruption* (the inverse test)."""


def fishers_frame_domain_test(
    contingency_table: Iterable[Iterable[int]],
) -> tuple[float, float]:
    """Run the pre-specified one-sided Fisher exact test.

    Rows are ``(configured hit-status, everything else)`` -- e.g. ``(domain
    retained, domain lost/disrupted)`` for the retention test, or ``(domain
    lost/disrupted, domain retained)`` for the disruption test -- and
    columns are ``(in-frame protein fusion, other)``.  A zero off-diagonal
    cell produces an infinite odds ratio in Fisher's exact test; this is
    expected for the Zehir paper's reported 33/33 pattern, not a numerical
    bug.  SciPy still calculates the exact finite p-value for that table.
    """
    table = np.asarray(list(contingency_table), dtype=int)
    if table.shape != (2, 2):
        raise ValueError("contingency_table must be a 2x2 table")
    if np.any(table < 0):
        raise ValueError("contingency_table counts must be non-negative")

    odds_ratio, p_value = fisher_exact(table, alternative="greater")
    return float(odds_ratio), float(p_value)


def build_frame_domain_contingency_table(
    events: list[FusionEvent],
    features: list[FusionFeature],
    gene_config: GeneConfig,
    *,
    domains: list[KeyDomain] | None = None,
    hit_statuses: frozenset[str] = RETAINED_STATUSES,
) -> list[list[int]]:
    """Build the Fisher table for a configured target-domain flag.

    ``domains`` defaults to ``gene_config.key_domains`` (the retention
    test); pass ``gene_config.disruption_required_domains`` with
    ``hit_statuses=DISRUPTED_STATUSES`` for the inverse disruption test.
    Row 0 is "domain status is in ``hit_statuses``", row 1 is everything
    else with a known status.

    Features with an unknown domain state are excluded: they carry no
    retention/loss observation and assigning them to either cell would create
    artificial evidence.  Events not marked as an in-frame protein fusion
    form the comparison column.
    """
    domains = gene_config.key_domains if domains is None else domains
    target_key = _target_domain_key(domains, gene_config.gene_symbol)
    event_by_id = {event.Event_id: event for event in events}
    counts = [[0, 0], [0, 0]]
    for feature in _target_features(features, gene_config):
        event = event_by_id.get(feature.Event_id)
        status = _domain_status(feature, target_key)
        if event is None or status is None:
            continue
        row = 0 if status in hit_statuses else 1
        column = 0 if _is_in_frame_protein_fusion(event) else 1
        counts[row][column] += 1
    return counts


def permutation_null_test(
    events: list[FusionEvent],
    features: list[FusionFeature],
    gene_config: GeneConfig,
    seed: int = 42,
    n_permutations: int = 10_000,
    *,
    domains: list[KeyDomain] | None = None,
    hit_statuses: frozenset[str] = RETAINED_STATUSES,
    genome_nexus_client: GenomeNexusClient | None = None,
) -> tuple[float, float, tuple[float, ...]]:
    """Estimate a one-sided empirical p-value by breakpoint randomization.

    ``domains`` and ``hit_statuses`` select which test is run, matching
    :func:`build_frame_domain_contingency_table` (retention by default;
    pass ``gene_config.disruption_required_domains`` and
    ``hit_statuses=DISRUPTED_STATUSES`` for the disruption test).

    The returned tuple is ``(empirical_p_value, observed_hit_rate,
    null_hit_rates)``.  If a ``genome_nexus_client`` is supplied, each
    random breakpoint is sampled across the target transcript's real genomic
    span (including its observed intron/exon structure) and mapped back to a
    CDS/protein position with the reviewed Genome Nexus arithmetic.  Offline
    callers can omit the client; then the observed resolved protein positions
    are resampled as a deterministic, conservative fallback.

    Genome Nexus describes exon coordinates but ``FusionFeature`` deliberately
    stores only the final domain call, so the null classifier is fitted from
    the available mapped breakpoint/domain observations by nearest resolved
    protein position.  This keeps the statistic usable for pre-annotated
    fixture data while avoiding another network-only domain lookup.
    """
    if n_permutations <= 0:
        raise ValueError("n_permutations must be positive")

    domains = gene_config.key_domains if domains is None else domains
    target_key = _target_domain_key(domains, gene_config.gene_symbol)
    event_by_id = {event.Event_id: event for event in events}
    records = [
        (event_by_id[feature.Event_id], feature, _domain_status(feature, target_key))
        for feature in _target_features(features, gene_config)
        if feature.Event_id in event_by_id and _domain_status(feature, target_key) is not None
    ]
    in_frame_records = [record for record in records if _is_in_frame_protein_fusion(record[0])]
    if not in_frame_records:
        raise ValueError("no in-frame protein-fusion events have a known domain-retention state")

    observed_hits = sum(status in hit_statuses for _, _, status in in_frame_records)
    observed_rate = observed_hits / len(in_frame_records)
    reference_positions, reference_statuses = _reference_positions(records)
    rng = np.random.default_rng(seed)

    genomic_mapper = _genome_nexus_mapper(gene_config, genome_nexus_client)
    null_rates: list[float] = []
    for _ in range(n_permutations):
        if genomic_mapper is None:
            sampled_positions = rng.choice(
                reference_positions, size=len(in_frame_records), replace=True
            )
        else:
            sampled_positions = genomic_mapper(rng, len(in_frame_records))
        hits = sum(
            _nearest_status(int(position), reference_positions, reference_statuses) in hit_statuses
            for position in sampled_positions
        )
        null_rates.append(hits / len(in_frame_records))

    empirical_p_value = (1 + sum(rate >= observed_rate for rate in null_rates)) / (
        n_permutations + 1
    )
    return float(empirical_p_value), float(observed_rate), tuple(null_rates)


def gene_breakpoint_domain_status_records(
    events: list[FusionEvent], features: list[FusionFeature], gene_config: GeneConfig
) -> list[tuple[int, str]]:
    """Return ``(breakpoint_protein_position, domain_status)`` pairs for a gene.

    One entry per feature belonging to ``gene_config.gene_symbol`` that has
    both a mapped breakpoint protein position and a known status (retained
    or lost/disrupted) for the gene's configured key domain. This is the
    same underlying (breakpoint, domain-status) observation the
    domain-retention algorithm consumes, exposed independent of frame
    status or partner-gene identity so other gene-agnostic algorithms (e.g.
    cutpoint scanning) can reuse it without duplicating the domain-key/
    status-extraction logic.
    """
    target_key = _target_domain_key(gene_config.key_domains, gene_config.gene_symbol)
    records: list[tuple[int, str]] = []
    for feature in _target_features(features, gene_config):
        status = _domain_status(feature, target_key)
        if status is None or feature.Junction_position_aa is None:
            continue
        records.append((feature.Junction_position_aa, status))
    return records


def gene_breakpoint_domain_status_event_records(
    events: list[FusionEvent], features: list[FusionFeature], gene_config: GeneConfig
) -> list[tuple[str, int, str]]:
    """Return ``(event_id, breakpoint_protein_position, domain_status)`` triples.

    Same underlying observation as :func:`gene_breakpoint_domain_status_records`,
    with the owning ``Event_id`` carried alongside so a caller can identify
    *which* events fall inside a candidate region -- needed by window-based
    scans (e.g. ``window_detection``) to de-duplicate candidate windows by
    event membership rather than by numeric position alone.
    """
    target_key = _target_domain_key(gene_config.key_domains, gene_config.gene_symbol)
    records: list[tuple[str, int, str]] = []
    for feature in _target_features(features, gene_config):
        status = _domain_status(feature, target_key)
        if status is None or feature.Junction_position_aa is None:
            continue
        records.append((feature.Event_id, feature.Junction_position_aa, status))
    return records


def domain_retention_descriptive_table(
    features: list[FusionFeature],
    gene_config: GeneConfig,
    *,
    domains: list[KeyDomain] | None = None,
) -> list[dict]:
    """Summarize optional quantitative retention without changing test calls.

    Only features carrying the additive ``Domain_retention_details`` data
    contribute. This keeps pre-existing/externally constructed features valid
    and leaves the contingency-table classification entirely flag-driven.
    """
    domains = gene_config.key_domains if domains is None else domains
    rows: list[dict] = []
    target_features = list(_target_features(features, gene_config))
    for domain in domains:
        domain_key = domain.key or domain.name.lower().replace(" ", "_")
        fractions: list[float] = []
        non_retained_fractions: list[float] = []
        truncated_count = 0
        fully_retained_count = 0
        fully_lost_count = 0
        for feature in target_features:
            detail = (feature.Domain_retention_details or {}).get(domain_key)
            if detail is None or detail.Retained_fraction is None:
                continue
            fraction = detail.Retained_fraction
            fractions.append(fraction)
            if detail.Is_truncated:
                truncated_count += 1
            elif fraction == 1.0:
                fully_retained_count += 1
            elif fraction == 0.0:
                fully_lost_count += 1
            status = (feature.Domain_retention_flags or {}).get(domain_key)
            if status in _NON_RETAINED:
                non_retained_fractions.append(fraction)
        rows.append(
            {
                "Domain_key": domain_key,
                "Domain_name": domain.name,
                "Quantitative_call_count": len(fractions),
                "Fully_retained_count": fully_retained_count,
                "Truncated_count": truncated_count,
                "Fully_lost_count": fully_lost_count,
                "Mean_retained_fraction_among_non_retained_calls": (
                    sum(non_retained_fractions) / len(non_retained_fractions)
                    if non_retained_fractions
                    else None
                ),
            }
        )
    return rows


def _target_domain_key(domains: list[KeyDomain], gene_symbol: str | None) -> str:
    if not domains:
        raise ValueError(f"{gene_symbol} has no configured domains for this test")
    domain = domains[0]
    return domain.key or domain.name.lower().replace(" ", "_")


def _target_features(features: list[FusionFeature], gene_config: GeneConfig):
    return (
        feature for feature in features if feature.Gene.upper() == gene_config.gene_symbol.upper()
    )


def _domain_status(feature: FusionFeature, domain_key: str) -> str | None:
    status = (feature.Domain_retention_flags or {}).get(domain_key)
    return status if status in _KNOWN_STATUSES else None


def _is_in_frame_protein_fusion(event: FusionEvent) -> bool:
    return event.Is_protein_fusion is True and event.Frame_status == "in-frame"


def _reference_positions(
    records: list[tuple[FusionEvent, FusionFeature, str]],
) -> tuple[np.ndarray, list[str]]:
    positions: list[int] = []
    statuses: list[str] = []
    for _, feature, status in records:
        if feature.Junction_position_aa is not None:
            positions.append(feature.Junction_position_aa)
            statuses.append(status)
    if not positions:
        raise ValueError("permutation test requires at least one mapped protein breakpoint")
    return np.asarray(positions, dtype=int), statuses


def _nearest_status(position: int, reference_positions: np.ndarray, statuses: list[str]) -> str:
    return statuses[int(np.argmin(np.abs(reference_positions - position)))]


def _genome_nexus_mapper(gene_config: GeneConfig, client: GenomeNexusClient | None):
    """Return a seeded random genomic-to-protein sampler, if a client is supplied."""
    if client is None:
        return None
    payload = client.fetch_canonical_transcript(gene_config.gene_symbol)
    transcript = parse_canonical_transcript(payload)
    if not transcript.exons:
        return None
    strand = transcript.exons[0].strand
    cds_min, cds_max = cds_bounds_from_utrs(transcript.utrs)
    genomic_start = min(exon.start for exon in transcript.exons)
    genomic_end = max(exon.end for exon in transcript.exons)

    def _sample(rng: np.random.Generator, size: int) -> np.ndarray:
        genomic_positions = rng.integers(genomic_start, genomic_end + 1, size=size)
        return np.asarray(
            [
                map_genomic_breakpoint_to_protein_position(
                    transcript.exons,
                    int(position),
                    strand,
                    cds_min_genomic=cds_min,
                    cds_max_genomic=cds_max,
                ).protein_position
                for position in genomic_positions
            ],
            dtype=int,
        )

    return _sample
