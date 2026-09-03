"""Built-in hotspot-detection algorithm plugins."""

from cfh.algorithms.exon_retention import ExonRetentionAnalysis
from cfh.algorithms.frequency import FrequencyAnalysis

__all__ = ["ExonRetentionAnalysis", "FrequencyAnalysis"]
