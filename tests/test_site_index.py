import json
import re

import pytest

from cfh.reporting.site_index import collect_site_index, render_site_index, study_label


@pytest.mark.parametrize(
    ("study_id", "label"),
    [
        ("msk_impact_50k_2026", "MSK-IMPACT 50K (2026)"),
        ("msk_impact_2017", "MSK-IMPACT (2017)"),
        ("luad_tcga_pan_can_atlas_2018", "TCGA LUAD PanCancer Atlas"),
        ("some_other_study", "some_other_study"),
    ],
)
def test_study_labels(study_id, label):
    assert study_label(study_id) == label


def _cohort(runs, name, genes, *, reports=(), significant=(), controls=None):
    cohort = runs / name / "cohort_scan"
    cohort.mkdir(parents=True)
    (cohort / "summary.json").write_text(
        json.dumps(
            {
                "study_id": "msk_impact_50k_2026",
                "significant_genes": list(significant),
                "honorable_mentions": [{"gene_symbol": "HMGENE"}],
                "genes": [
                    {
                        "gene_symbol": gene,
                        "status": "failed" if gene == "BROKEN" else "ok",
                        "n_events_analyzed": 10,
                        "distinct_patient_count": 9,
                        "min_fdr_adjusted_q_value": 0.001 if gene in significant else 0.5,
                        "fdr_significant": gene in significant,
                    }
                    for gene in genes
                ],
            }
        )
    )
    for gene in reports:
        (cohort / "gene_reports" / gene.lower()).mkdir(parents=True)
        (cohort / "gene_reports" / gene.lower() / "results.json").write_text("{}")
    if controls:
        rows = ["control\tcheck\tstatus\tdetail\tsource"] + [
            f"G{i}\tselection\t{status}\t-\t-" for i, status in enumerate(controls)
        ]
        (cohort / "positive_controls.tsv").write_text("\n".join(rows) + "\n")


def _single(runs, name, gene, study, *, pair=None):
    (runs / name).mkdir(parents=True)
    summary = {"total_fusions": 5, "known_patient_count": 5}
    if pair:
        summary["gene_pair"] = pair
    (runs / name / "results.json").write_text(
        json.dumps({"gene_symbol": gene, "study_id": study, "summary": summary})
    )


@pytest.fixture
def runs(tmp_path):
    runs = tmp_path / "runs"
    _cohort(
        runs,
        "cohort-scan_msk_impact_50k_2026_20261026T200631Z",
        ["ALK", "TP53", "HMGENE", "BROKEN"],
        reports=["ALK"],
        significant=["ALK"],
        controls=["pass", "pass", "fail"],
    )
    _cohort(runs, "cohort-scan_msk_impact_50k_2026_20260915T163812Z", ["ALK"], reports=["ALK"])
    _single(runs, "alk_msk-impact-50k-2026_20261029T153650Z", "ALK", "msk_impact_50k_2026")
    luad = "luad_tcga_pan_can_atlas_2018"
    _single(runs, "alk_luad-tcga-pan-can-atlas-2018_20261029T153648Z", "ALK", luad)
    _single(runs, "ros1_luad-tcga-pan-can-atlas-2018_20261029T153648Z", "ROS1", luad)
    _single(
        runs,
        "eml4-alk_msk-impact-50k-2026_20261029T153713Z",
        "EML4-ALK",
        "msk_impact_50k_2026",
        pair=["EML4", "ALK"],
    )
    return runs


def _published(runs):
    return sorted(path.name for path in runs.iterdir())


def test_each_gene_goes_to_its_most_informative_page(runs):
    index = collect_site_index(runs, _published(runs))
    genes = {gene["gene"]: gene for gene in index["genes"]}

    alk = genes["ALK"]
    assert alk["best"]["url"] == (
        "runs/cohort-scan_msk_impact_50k_2026_20261026T200631Z/genes/alk.html"
    )
    assert alk["best"]["label"] == "Full gene report"
    assert [(entry["kind"], entry["study"]) for entry in alk["also"]] == [
        ("gene", "MSK-IMPACT 50K (2026)"),  # same study as the newest scan
        ("gene", "TCGA LUAD PanCancer Atlas"),
        ("cohort", "MSK-IMPACT 50K (2026)"),  # the older scan, after replications
        ("pair", "MSK-IMPACT 50K (2026)"),
    ]
    assert alk["also"][2]["date"] == "2026-09-15"
    assert alk["significant"] is True and alk["q"] == 0.001

    # Analyzed without a full report: the cohort summary.
    assert genes["TP53"]["best"]["url"] == (
        "runs/cohort-scan_msk_impact_50k_2026_20261026T200631Z/index.html"
    )
    assert genes["TP53"]["best"]["label"] == "Cohort summary"
    assert genes["HMGENE"]["honorable"] is True
    assert genes["BROKEN"]["analyzed"] is False
    # Only in a replication cohort, and only in a pair.
    assert genes["ROS1"]["best"]["study"] == "TCGA LUAD PanCancer Atlas"
    assert genes["EML4"]["best"]["label"] == "EML4–ALK pair run"


def test_runs_list_newest_cohort_first_with_its_headline_numbers(runs):
    index = collect_site_index(runs, _published(runs))
    first = index["runs"][0]
    assert first["kind"] == "cohort" and first["newest"] is True
    assert first["genes_analyzed"] == 3
    assert first["positive_controls"] == [2, 3]
    assert [run["kind"] for run in index["runs"]][:2] == ["cohort", "cohort"]
    assert index["runs"][1]["newest"] is False


def test_landing_page_embeds_the_index_and_escapes_it(runs):
    index = collect_site_index(runs, _published(runs))
    index["genes"][0]["best"]["label"] = "</script><b>x"
    html = render_site_index(index)

    assert "<title>Fusion Hotspots</title>" in html
    assert '<div class="stat-value">3</div>' in html  # genes analysed
    assert '<div class="stat-value">2 of 3</div>' in html  # textbook fusions recovered
    assert '<button type="button" class="chip" data-gene="ALK">ALK</button>' in html
    assert 'id="gene-search"' in html
    data = re.search(r'id="gene-data">(.*?)</script>', html, re.S).group(1)
    assert "</script>" not in data
    assert json.loads(data)[0]["best"]["label"] == "</script><b>x"
    assert html.count('class="card card-') == len(index["runs"])
