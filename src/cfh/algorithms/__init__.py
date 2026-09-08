"""Built-in hotspot-detection algorithm plugins.

Importing this package registers every built-in algorithm plugin via each
module's ``@register`` decorator.
"""

from cfh.algorithms.composite_score import CompositeScoreAlgorithm
from cfh.algorithms.confidence_stats import ConfidenceStatsAlgorithm
from cfh.algorithms.cutpoint_detection import CutpointDetectionAlgorithm
from cfh.algorithms.domain_disruption import DomainDisruptionAlgorithm
from cfh.algorithms.domain_retention import DomainRetentionAlgorithm
from cfh.algorithms.exon_retention import ExonRetentionAnalysis
from cfh.algorithms.frequency import FrequencyAnalysis
from cfh.algorithms.joint_partner import JointPartnerMode
from cfh.algorithms.window_detection import WindowDetectionAlgorithm

__all__ = [
    "CompositeScoreAlgorithm",
    "ConfidenceStatsAlgorithm",
    "CutpointDetectionAlgorithm",
    "DomainDisruptionAlgorithm",
    "DomainRetentionAlgorithm",
    "ExonRetentionAnalysis",
    "FrequencyAnalysis",
    "JointPartnerMode",
    "WindowDetectionAlgorithm",
]
