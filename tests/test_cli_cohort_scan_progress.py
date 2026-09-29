"""``cfh cohort-scan`` wires per-gene progress to stderr unless --quiet."""

from __future__ import annotations

from types import SimpleNamespace

import pytest
from click.testing import CliRunner

from cfh import cli


def _fake_result():
    return SimpleNamespace(
        gene_outcomes=[],
        total_genes_before_gating=0,
        genes_after_gating=0,
        curated_gene_count=0,
        auto_config_gene_count=0,
        unresolved_gene_count=0,
        non_coding_gene_count=0,
        significant_genes=[],
        warnings=[],
    )


@pytest.mark.parametrize(("flag", "expect_progress"), [([], True), (["--quiet"], False)])
def test_progress_flag(monkeypatch, tmp_path, flag, expect_progress):
    seen = {}

    def fake_scan(study_id, **kwargs):
        seen["progress"] = kwargs["progress"]
        if kwargs["progress"] is not None:
            kwargs["progress"]("[1/1] BRAF: ok (0.1s; 3 events)")
        return _fake_result()

    monkeypatch.setattr(cli, "run_cohort_scan", fake_scan)
    monkeypatch.setattr(cli, "write_cohort_scan_outputs", lambda *a, **k: {})
    result = CliRunner().invoke(
        cli.main, ["cohort-scan", "study", "--output-dir", str(tmp_path), "--no-pdf", *flag]
    )
    assert result.exit_code == 0, result.output
    assert (seen["progress"] is not None) is expect_progress
    assert ("[1/1] BRAF: ok" in result.stderr) is expect_progress
    assert "[1/1] BRAF" not in result.stdout
