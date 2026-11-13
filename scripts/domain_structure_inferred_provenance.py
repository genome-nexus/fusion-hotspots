"""ONE-OFF: mark the 2026-11-12 legacy checkpoints as inferred-provenance (see docs).

Usage::

    python scripts/domain_structure_inferred_provenance.py --work runs/.domain_structure_work

Never run by ``domain_structure.py predict`` or ``analyze``. Uses the default
PredictorConfig and the fixed evidence of
:data:`cfh.domain_structure.provenance.EVIDENCE_2026_11_12`; it has no option to
name another server, and refuses if anything disagrees with the evidence.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from cfh.domain_structure.predictor import PredictionFailure, PredictorConfig
from cfh.domain_structure.provenance import EVIDENCE_2026_11_12, record_inferred_provenance


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--work", type=Path, default=Path("runs/.domain_structure_work"))
    args = parser.parse_args()
    try:
        record_inferred_provenance(args.work, PredictorConfig(), EVIDENCE_2026_11_12)
    except PredictionFailure as error:
        print(f"REFUSED: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
