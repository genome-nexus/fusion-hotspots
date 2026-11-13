"""Protein-language-model scores for chimeric fused domains (PR 3 of the domain-function stack).

Consumes the records written by :mod:`cfh.domain_function` (PR #155) and the
structure flags written by :mod:`cfh.domain_structure` (PR #160). Every fused
domain and both parent domain windows are scored with one ESM-2 model by masked
pseudo-log-likelihood (one residue masked per forward pass):

* :mod:`.worker` -- standalone ESM-2 process, run with the external PLM env's
  Python (torch is never imported by the ``cfh`` package itself);
* :mod:`.cache` -- per-sequence score cache keyed by model name + version +
  weights hash + sequence hash; every entry is re-validated on read (fail closed);
* :mod:`.scoring` -- positions to score per sequence and the worker client;
* :mod:`.metrics` -- length-normalised chimera-vs-parent comparisons and the
  junction dip;
* :mod:`.compare` -- a descriptive rank label and the disagreement list against
  PR #155 tiers and PR #160 fold flags.

The scores are reported next to PR #155's tier and PR #160's flag; they never
change either. No combined score and no functional cutoff is defined.
"""

SCHEMA_VERSION = "domain_plm/1"
