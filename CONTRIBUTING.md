# Contributing to Cancer Fusion Hotspots

Thanks for your interest in contributing.

## Development setup

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

## Running checks locally

```bash
ruff check .
ruff format --check .
mypy src
pytest -m "not network"
```

The same lint and formatting checks can run automatically before each commit:

```bash
pre-commit install
pre-commit run --all-files
```

Tests marked `@pytest.mark.network` make real calls to external services
(cBioPortal, Ensembl, UniProt, InterPro) and are excluded from the default
test run. A weekly and manually dispatchable workflow runs them against the
live services. To run them locally, opt in explicitly:

```bash
CFH_RUN_NETWORK_TESTS=1 pytest -m network
```

## Adding support for a new gene

Gene biology (canonical transcript, protein accession, key domains,
autoinhibitory domains) belongs in a new YAML file under
`src/cfh/genes/configs/`, loaded via `GeneConfig`. Generic
ingestion/normalization/mapping code must never hardcode a specific gene
symbol, transcript ID, or protein accession — those always come from the
`GeneConfig` passed in by the caller.

## Adding a hotspot-detection algorithm

Implement the `Algorithm` interface in `src/cfh/algorithms/base.py` and
register it via `src/cfh/algorithms/registry.py` so it can be discovered by
name.

## Committing run artifacts

Some tests and benchmarks reuse real, already-computed run artifacts
committed under `runs/<type>_<ISO8601-timestamp>/` (e.g.
`runs/braf_msk-impact-50k-2026_20260904T172738Z/`), where `<type>` is the
gene-or-scan-kind + study_id combination (e.g. `braf_msk-impact-50k-2026`,
`cohort-scan_msk_impact_50k_2026`). Regenerating one of these (e.g. after a
fix round) produces a new timestamped directory alongside the old one.

**Before committing a new run artifact under `runs/<type>_<timestamp>/`,
remove any existing `runs/<type>_*/` directory for the same
gene-or-scan-kind + study_id combination so only the latest run per type is
ever committed.** Run the helper script to do this automatically:

```bash
python scripts/prune_old_runs.py            # dry run: reports what would be removed
python scripts/prune_old_runs.py --apply    # actually deletes the stale directories
```

It groups `runs/` directories by everything before the trailing
`_<timestamp>` suffix and removes every directory in a group except the one
with the latest timestamp; it never touches standalone files directly under
`runs/`. Run it (with `--apply`) before committing new run artifacts, and
include the removals in the same commit/PR.

## Pull requests

- Keep changes focused and covered by tests.
- Run `pre-commit run --all-files`, `mypy src`, and
  `pytest -m "not network"` before opening a PR.
- Describe what changed and how it was verified.
