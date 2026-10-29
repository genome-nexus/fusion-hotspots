import json
import re

from cfh.ingestion.cbioportal_api import fusion_viewer_id
from cfh.reporting.html_viewer import (
    CBIOPORTAL_URL,
    _cbioportal_case_link,
    _cbioportal_fusion_link,
    _render_events_table,
    build_run_viewer,
)
from conftest import latest_run_dir


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
    assert ">Sample</button></th>" in html and ">Patient</button></th>" in html
    assert "studyId=study_x&amp;sampleId=S1" in html
    assert "studyId=study_x&amp;caseId=P1" in html


def test_fusion_viewer_id_matches_cbioportal_patient_view():
    # The deep link the cBioPortal fusion viewer itself produces for this
    # GABPA::ERG fusion (site 1 is GABPA, site 2 ERG).
    row = {
        "Sample_Id": "P-0023818-T01-IM6",
        "Site1_Hugo_Symbol": "GABPA",
        "Site2_Hugo_Symbol": "ERG",
        "Site1_Position": 27139889.0,  # pandas upcasts positions with gaps
        "Site2_Position": 39755807,
    }
    assert fusion_viewer_id(row) == "P-0023818-T01-IM6_GABPA_ERG_27139889_39755807"
    # Like the viewer: a missing partner is blank, a missing position 0.
    row.update(Site2_Hugo_Symbol=float("nan"), Site2_Position=-1)
    assert fusion_viewer_id(row) == "P-0023818-T01-IM6_GABPA__27139889_0"
    assert fusion_viewer_id({"Sample_Id": None}) is None


def test_fusion_viewer_link_opens_the_fusion_or_the_sample():
    fusion_id = "P-0023818-T01-IM6_GABPA_ERG_27139889_39755807"
    link = _cbioportal_fusion_link("msk_impact_50k_2026", "P-0023818-T01-IM6", fusion_id)
    assert (
        'href="https://www.cbioportal.org/patient/fusionViewer?studyId=msk_impact_50k_2026'
        "&amp;sampleId=P-0023818-T01-IM6"
        '&amp;fusionId=P-0023818-T01-IM6_GABPA_ERG_27139889_39755807"'
    ) in link
    assert 'target="_blank" rel="noopener noreferrer"' in link
    # Runs that predate the recorded id still reach the sample's fusions.
    fallback = _cbioportal_fusion_link("msk_impact_50k_2026", "S1", None)
    assert 'fusionViewer?studyId=msk_impact_50k_2026&amp;sampleId=S1"' in fallback
    assert _cbioportal_fusion_link(None, "S1", fusion_id) == ""

    html = _render_events_table(
        [{"event_id": "E1", "sample_id": "S1", "cbioportal_fusion_id": "S1_A_B_1_2"}],
        "study_x",
    )
    assert ">cBioPortal</button></th>" in html
    assert "fusionId=S1_A_B_1_2" in html


def test_committed_alk_page_links_every_event_to_cbioportal(tmp_path):
    alk = latest_run_dir("cohort-scan_msk_impact_50k_2026") / "cohort_scan" / "gene_reports" / "alk"
    if not (alk / "results.json").exists():
        import pytest

        pytest.skip("committed cohort scan has no report for this gene")
    payload = json.loads((alk / "results.json").read_text())
    run_dir = tmp_path / "alk"
    run_dir.mkdir()
    (run_dir / "results.json").write_text(json.dumps(payload))
    html = build_run_viewer(run_dir).read_text()

    table = html.split('id="events-table"', 1)[1].split("</table>", 1)[0]
    rows = table.count("<tr data-tumor-type=")
    assert rows == len(payload["events"])
    assert len(re.findall(r"/patient/summary\?studyId=\S+&amp;sampleId=", table)) == rows
    assert table.count("/patient/fusionViewer?") == rows
    assert len(re.findall(r"&amp;caseId=", table)) == rows
    # A link in a clickable counter-intuitive row opens cBioPortal without
    # also jumping within the page.
    assert 'clickEvent.target.closest("a")' in html
