import json
import re
from pathlib import Path

from cfh.reporting.html_viewer import (
    CBIOPORTAL_URL,
    _cbioportal_case_link,
    _render_events_table,
    build_run_viewer,
)

COMMITTED_ALK = (
    Path(__file__).resolve().parents[1]
    / "runs"
    / "cohort-scan_msk_impact_50k_2026_20261026T200631Z"
    / "cohort_scan"
    / "gene_reports"
    / "alk"
)


def test_links_use_cbioportal_patient_view_urls():
    assert CBIOPORTAL_URL == "https://www.cbioportal.org"
    sample = _cbioportal_case_link("msk_impact_50k_2026", "P-0053737-T01-IM6", kind="sample")
    patient = _cbioportal_case_link("msk_impact_50k_2026", "P-0053737", kind="patient")
    assert (
        'href="https://www.cbioportal.org/patient/summary?studyId=msk_impact_50k_2026'
        '&amp;sampleId=P-0053737-T01-IM6"'
    ) in sample
    assert (
        'href="https://www.cbioportal.org/patient/summary?studyId=msk_impact_50k_2026'
        '&amp;caseId=P-0053737"'
    ) in patient
    for link in (sample, patient):
        assert 'target="_blank" rel="noopener noreferrer"' in link


def test_ids_are_plain_text_without_a_study_or_id_and_always_escaped():
    assert _cbioportal_case_link(None, "S1", kind="sample") == "S1"
    assert _cbioportal_case_link("study", "", kind="patient") == ""
    link = _cbioportal_case_link("my study", "S<1>&2", kind="sample")
    assert "studyId=my+study&amp;sampleId=S%3C1%3E%262" in link
    assert ">S&lt;1&gt;&amp;2</a>" in link


def test_events_table_links_sample_and_patient_columns():
    html = _render_events_table(
        [{"event_id": "E1", "sample_id": "S1", "patient_id": "P1", "domain_status": "lost"}],
        "study_x",
    )
    assert "<th>Sample</th><th>Patient</th>" in html
    assert "studyId=study_x&amp;sampleId=S1" in html
    assert "studyId=study_x&amp;caseId=P1" in html


def test_committed_alk_page_links_every_event_to_cbioportal(tmp_path):
    if not (COMMITTED_ALK / "results.json").exists():
        import pytest

        pytest.skip("committed 26 October run not present")
    payload = json.loads((COMMITTED_ALK / "results.json").read_text())
    run_dir = tmp_path / "alk"
    run_dir.mkdir()
    (run_dir / "results.json").write_text(json.dumps(payload))
    html = build_run_viewer(run_dir).read_text()

    table = html.split('id="events-table"', 1)[1].split("</table>", 1)[0]
    rows = table.count("<tr data-tumor-type=")
    assert rows == len(payload["events"])
    assert len(re.findall(r"&amp;sampleId=", table)) == rows
    assert len(re.findall(r"&amp;caseId=", table)) == rows
    # A link in a clickable counter-intuitive row opens cBioPortal without
    # also jumping within the page.
    assert 'clickEvent.target.closest("a")' in html
