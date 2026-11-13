"""Structure prediction of shortlisted chimeric domains (PR 2 of the domain-function stack).

Consumes the records and shortlist written by :mod:`cfh.domain_function`
(``domain_function_seq/2``). Each distinct chimeric domain (flanked window) and
both parent domain windows are predicted with identical ColabFold settings:

* :mod:`.jobs` -- one prediction job per distinct sequence;
* :mod:`.predictor` -- configurable backend, batched + cached MSAs, cached and
  validated predictions (fail closed);
* :mod:`.geometry` -- Kabsch, RMSD, TM-score with a fixed residue mapping;
* :mod:`.metrics` -- pLDDT / PAE / pTM per prediction and chimera-vs-parent
  comparisons (ColabFold parents and AlphaFold DB parents);
* :mod:`.flags` -- a HEURISTIC screen-out flag (``fold_supported`` /
  ``fold_doubtful`` / ``not_assessed``).

See ``docs/domain_structure.md``. PR 1's tiers are not changed.
"""

SCHEMA_VERSION = "domain_structure/1"
