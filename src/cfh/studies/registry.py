"""Configuration registry for study-specific source metadata."""

from __future__ import annotations

from pathlib import Path

import yaml
from pydantic import BaseModel, ConfigDict

CONFIGS_DIR = Path(__file__).parent / "configs"


class StudyConfig(BaseModel):
    """Source settings shared by one or more cBioPortal studies."""

    model_config = ConfigDict(extra="forbid")

    study_ids: list[str]
    structural_variant_profile_template: str = "{study_id}_structural_variants"
    mutation_profile_template: str = "{study_id}_mutations"
    discrete_cna_profile_template: str = "{study_id}_cna"
    all_sample_list_template: str = "{study_id}_all"
    genome_nexus_base_url: str = "https://www.genomenexus.org"
    mrna_expression_profile_template: str | None = None
    """Molecular-profile-id template for an mRNA expression assay (e.g. TCGA
    PanCancer Atlas's per-gene z-score profile), or ``None`` when this
    cohort has no such profile at all -- e.g. a targeted DNA panel like
    ``msk_impact_50k_2026`` (which has no study config entry here in the
    first place, so this never even applies). Opt-in, mirroring
    ``structural_variant_profile_template``; a cohort without it makes
    :meth:`mrna_expression_profile_id` return ``None`` rather than guess."""

    def molecular_profile_id(self, study_id: str) -> str:
        if study_id not in self.study_ids:
            raise ValueError(f"Study {study_id!r} is not covered by this config")
        return self.structural_variant_profile_template.format(study_id=study_id)

    def mutation_profile_id(self, study_id: str) -> str:
        if study_id not in self.study_ids:
            raise ValueError(f"Study {study_id!r} is not covered by this config")
        return self.mutation_profile_template.format(study_id=study_id)

    def discrete_cna_profile_id(self, study_id: str) -> str:
        if study_id not in self.study_ids:
            raise ValueError(f"Study {study_id!r} is not covered by this config")
        return self.discrete_cna_profile_template.format(study_id=study_id)

    def all_sample_list_id(self, study_id: str) -> str:
        if study_id not in self.study_ids:
            raise ValueError(f"Study {study_id!r} is not covered by this config")
        return self.all_sample_list_template.format(study_id=study_id)

    def mrna_expression_profile_id(self, study_id: str) -> str | None:
        if study_id not in self.study_ids:
            raise ValueError(f"Study {study_id!r} is not covered by this config")
        if self.mrna_expression_profile_template is None:
            return None
        return self.mrna_expression_profile_template.format(study_id=study_id)


def load_study_config(study_id: str) -> StudyConfig | None:
    """Return the unique config covering ``study_id``, or ``None`` for defaults."""
    matches = []
    for path in CONFIGS_DIR.glob("*.yaml"):
        with path.open() as handle:
            config = StudyConfig.model_validate(yaml.safe_load(handle))
        if study_id in config.study_ids:
            matches.append(config)
    if len(matches) > 1:
        raise ValueError(f"Multiple study configs cover {study_id!r}")
    return matches[0] if matches else None
