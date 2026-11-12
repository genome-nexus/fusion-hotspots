"""Sequence-level evidence that a chimeric (fusion-reconstituted) domain is functional.

The junction screen (:mod:`cfh.cohort.domain_screen`) says whether a shared Pfam
domain is rebuilt cleanly at a fusion junction. A non-clean junction can still
give a working domain, and a clean one can still break it, so this package
reports independent, established metrics for each distinct chimeric domain:

* :mod:`.record` -- one versioned record per chimeric domain sequence, the
  shared input for later structure-prediction and stability steps;
* :mod:`.pfam_hmm` -- profile-HMM integrity with HMMER (pyhmmer) against the
  family's Pfam HMM and its gathering threshold;
* :mod:`.uniprot` / :mod:`.afdb` -- parent key-residue annotations and parent
  AlphaFold DB models (DSSP secondary structure, pLDDT, residue contacts);
* :mod:`.metrics` -- key-residue retention, junction secondary structure and a
  SCHEMA-like broken-contact count;
* :mod:`.tiers` -- a transparent, documented, HEURISTIC rule-based tier.

Nothing here is a calibrated predictor of activity; see
``docs/domain_function_seq.md``.
"""

SCHEMA_VERSION = "domain_function_seq/2"
"""Version of the committed record/metrics layout (v2: origin-aware record ids,
observations, sample-level read support, structure and annotation statuses)."""

RECORDS_VERSION = 2
"""Bumped whenever record building changes; part of the records-cache fingerprint."""
