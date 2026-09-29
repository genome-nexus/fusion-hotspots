"""Collapse repeated observations of one fusion in one patient.

Repeat biopsies and re-sequenced samples from the same patient report the same
rearrangement several times. Every domain and spatial test in this package
treats its inputs as exchangeable, independent observations; counting the same
patient's fusion two or three times violates that assumption and makes the
permutation nulls anti-conservative (see ``CALIBRATION.md``: the
repeated-patient null rejected at 60% for cutpoint and 88% for window at a
nominal 5%).

The ``patient`` unit keeps one representative per
``(cohort, patient, partner gene, target role, junction protein position)``.
Distinct rearrangements in one patient remain separate observations, so this
unit removes repeated measurement, not every possible within-patient
dependence. Events without a patient ID fall back to their sample ID, then to
their event ID, and are counted in the report.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from cfh.model.fusion_event import FusionEvent
from cfh.model.fusion_feature import FusionFeature

OBSERVATION_UNITS = ("patient", "event")
DEFAULT_OBSERVATION_UNIT = "patient"

_DETERMINATE_FRAMES = ("in-frame", "out-of-frame")


@dataclass
class CollapseReport:
    observation_unit: str
    input_event_count: int
    output_event_count: int
    collapsed_event_count: int = 0
    patients_with_repeats: int = 0
    identity_fallback_counts: dict[str, int] = field(default_factory=dict)
    frame_conflict_groups: int = 0

    def as_summary(self) -> dict:
        return {
            "observation_unit": self.observation_unit,
            "input_mapped_event_count": self.input_event_count,
            "analyzed_observation_count": self.output_event_count,
            "collapsed_repeat_observation_count": self.collapsed_event_count,
            "patients_with_repeated_observations": self.patients_with_repeats,
            "observation_identity_fallbacks": dict(self.identity_fallback_counts),
            "frame_conflict_groups": self.frame_conflict_groups,
        }


def _identity(event: FusionEvent) -> tuple[str, str]:
    if event.Patient_id:
        return "patient", event.Patient_id
    if event.Sample_id:
        return "sample", event.Sample_id
    return "event", event.Event_id


def _partner(event: FusionEvent, target_gene: str) -> str:
    target = target_gene.upper()
    genes = [gene for gene in (event.Site1_gene, event.Site2_gene) if gene]
    partners = sorted({gene.upper() for gene in genes if gene.upper() != target})
    return partners[0] if partners else "unknown"


def _representative_rank(event: FusionEvent) -> tuple[int, str, str]:
    """Prefer a determinate frame call, then the earliest sample and event ID."""
    frame_rank = 0 if event.Frame_status in _DETERMINATE_FRAMES else 1
    return frame_rank, event.Sample_id or "", event.Event_id


def collapse_repeated_observations(
    events: list[FusionEvent],
    features: list[FusionFeature],
    target_gene: str,
    *,
    observation_unit: str = DEFAULT_OBSERVATION_UNIT,
) -> tuple[list[FusionEvent], list[FusionFeature], CollapseReport]:
    """Return the events and features to test under ``observation_unit``.

    ``events`` and ``features`` are the mapped, index-aligned pairs used by
    the domain and spatial tests. The returned lists preserve input order.
    ``observation_unit="event"`` returns the inputs unchanged.
    """
    if observation_unit not in OBSERVATION_UNITS:
        raise ValueError(f"observation_unit must be one of {OBSERVATION_UNITS}")
    if len(events) != len(features):
        raise ValueError("events and features must be index-aligned")
    if observation_unit == "event":
        return (
            list(events),
            list(features),
            CollapseReport("event", len(events), len(events)),
        )

    groups: dict[tuple, list[int]] = {}
    fallbacks: dict[str, int] = {}
    for index, (event, feature) in enumerate(zip(events, features, strict=True)):
        if event.Event_id != feature.Event_id:
            raise ValueError("events and features must be index-aligned")
        kind, identity = _identity(event)
        if kind != "patient":
            fallbacks[kind] = fallbacks.get(kind, 0) + 1
        key = (
            event.Cohort,
            kind,
            identity,
            _partner(event, target_gene),
            feature.Role,
            feature.Junction_position_aa,
        )
        groups.setdefault(key, []).append(index)

    keep: set[int] = set()
    patients_with_repeats: set[tuple[str, str]] = set()
    frame_conflicts = 0
    for (cohort, kind, identity, *_), indices in groups.items():
        keep.add(min(indices, key=lambda i: _representative_rank(events[i])))
        if len(indices) > 1:
            patients_with_repeats.add((cohort, f"{kind}:{identity}"))
            frames = {events[i].Frame_status for i in indices} & set(_DETERMINATE_FRAMES)
            frame_conflicts += len(frames) > 1

    kept = sorted(keep)
    report = CollapseReport(
        observation_unit="patient",
        input_event_count=len(events),
        output_event_count=len(kept),
        collapsed_event_count=len(events) - len(kept),
        patients_with_repeats=len(patients_with_repeats),
        identity_fallback_counts=fallbacks,
        frame_conflict_groups=frame_conflicts,
    )
    return [events[i] for i in kept], [features[i] for i in kept], report
