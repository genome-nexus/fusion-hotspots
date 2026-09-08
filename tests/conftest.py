import socket
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).parent.parent
RUNS_DIR = REPO_ROOT / "runs"
FIXTURES_DIR = Path(__file__).parent / "fixtures"
SV_FIXTURE_DIR = FIXTURES_DIR / "sv"
SV_FIXTURE_FILE = SV_FIXTURE_DIR / "data_sv.txt"
SV_SURPLUS_FIELDS_FIXTURE_FILE = SV_FIXTURE_DIR / "surplus_fields.txt"
CLINICAL_SAMPLE_FIXTURE = FIXTURES_DIR / "clinical" / "data_clinical_sample.txt"
CLINICAL_PATIENT_FIXTURE = FIXTURES_DIR / "clinical" / "data_clinical_patient.txt"
UNIPROT_FIXTURE = FIXTURES_DIR / "uniprot" / "P15056.json"
GENOME_NEXUS_CANONICAL_TRANSCRIPT_FIXTURE = (
    FIXTURES_DIR / "genome_nexus" / "canonical_transcript_braf.json"
)
GENOME_NEXUS_TRANSCRIPT_FIXTURE = FIXTURES_DIR / "genome_nexus" / "transcript_ENST00000288602.json"
GENOME_NEXUS_CANONICAL_TRANSCRIPT_PIK3CA_FIXTURE = (
    FIXTURES_DIR / "genome_nexus" / "canonical_transcript_pik3ca.json"
)


@pytest.fixture
def sv_fixture_dir() -> Path:
    return SV_FIXTURE_DIR


@pytest.fixture
def sv_fixture_file() -> Path:
    return SV_FIXTURE_FILE


@pytest.fixture
def sv_surplus_fields_fixture_file() -> Path:
    return SV_SURPLUS_FIELDS_FIXTURE_FILE


@pytest.fixture
def clinical_sample_fixture() -> Path:
    return CLINICAL_SAMPLE_FIXTURE


@pytest.fixture
def clinical_patient_fixture() -> Path:
    return CLINICAL_PATIENT_FIXTURE


@pytest.fixture
def uniprot_fixture_path() -> Path:
    return UNIPROT_FIXTURE


@pytest.fixture
def genome_nexus_canonical_transcript_fixture_path() -> Path:
    return GENOME_NEXUS_CANONICAL_TRANSCRIPT_FIXTURE


@pytest.fixture
def genome_nexus_transcript_fixture_path() -> Path:
    return GENOME_NEXUS_TRANSCRIPT_FIXTURE


@pytest.fixture
def genome_nexus_canonical_transcript_pik3ca_fixture_path() -> Path:
    return GENOME_NEXUS_CANONICAL_TRANSCRIPT_PIK3CA_FIXTURE


def latest_run_dir(prefix: str) -> Path:
    """Resolve the most recently generated committed run directory under
    ``runs/`` matching ``<prefix>_<ISO8601-timestamp>`` (this repo's
    keep-latest-run convention -- see CONTRIBUTING.md and
    ``scripts/prune_old_runs.py``).

    Tests must never hardcode a specific run timestamp directly: regenerating
    a run (e.g. after a fix round) prunes the old timestamped directory and
    commits a new one, which silently breaks any test pinning the old name
    with a ``FileNotFoundError`` unrelated to the change that prompted the
    regeneration. This glob-based lookup is immune to that by construction.
    """
    candidates = sorted(RUNS_DIR.glob(f"{prefix}_*"))
    if not candidates:
        pytest.skip(f"no committed real run directory found for {prefix!r} under {RUNS_DIR}")
    return candidates[-1]


@pytest.fixture(autouse=True)
def _block_real_network_calls(request, monkeypatch):
    """Repo-wide guard: any test not marked ``@pytest.mark.network`` that
    tries to open a real socket fails immediately, instead of silently
    succeeding against a live external service.
    """
    if request.node.get_closest_marker("network"):
        yield
        return

    def _blocked(*args, **kwargs):
        raise AssertionError(
            "real network access attempted in a test not marked @pytest.mark.network"
        )

    monkeypatch.setattr(socket, "create_connection", _blocked)
    yield
