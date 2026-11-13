"""Structure step (cfh.domain_structure): no ColabFold, no network unless marked."""

from __future__ import annotations

import gzip
import json
import shutil
from pathlib import Path

import numpy as np
import pytest

from cfh.domain_function.record import ChimeraRecord, ParentDomain
from cfh.domain_structure.flags import (
    FOLD_DOUBTFUL,
    FOLD_SUPPORTED,
    NOT_ASSESSED,
    RULES,
    assign_flag,
)
from cfh.domain_structure.geometry import kabsch, superposed_rmsd, tm_d0, tm_score_fixed
from cfh.domain_structure.jobs import Job, build_jobs
from cfh.domain_structure.metrics import (
    AfdbModel,
    Window,
    chimera_halves,
    compare,
    junction_window,
    parent_halves,
    record_metrics,
)
from cfh.domain_structure.predictor import (
    DiskSpaceError,
    MsaCache,
    Prediction,
    PredictionFailure,
    Predictor,
    PredictorConfig,
    a3m_query,
    check_disk,
    read_ca,
    validate_outputs,
)
from cfh.domain_structure.run import assess, export_models, flag_counts

DATA = Path(__file__).parent / "data" / "domain_structure"
FLI1_ERG = (
    "ISKNTEQRPQPDPYQILGPTSSRLANPGSGQIQLWQFLLELLSDSANASCITWEGTNGEFKMTDPDEVARRWGERKSKPNMN"
    "YDKLSRALRYYYDKNIYYYDKNIMTKVHGKRYAYKFDFHGIAQALQPHPPESSLYKYPSDLPYMGSYH"
)


# --------------------------------------------------------------------------- geometry
def _helix(n: int) -> np.ndarray:
    t = np.arange(n)
    return np.stack([2.3 * np.cos(t * 1.745), 2.3 * np.sin(t * 1.745), 1.5 * t], axis=1)


def _rotate(coords: np.ndarray) -> np.ndarray:
    angle = 0.7
    rotation = np.array(
        [[np.cos(angle), -np.sin(angle), 0], [np.sin(angle), np.cos(angle), 0], [0, 0, 1]]
    )
    return coords @ rotation.T + np.array([5.0, -3.0, 12.0])


def test_kabsch_recovers_rigid_motion():
    coords = _helix(30)
    moved = _rotate(coords)
    rotation, translation = kabsch(coords, moved)
    assert np.allclose(coords @ rotation.T + translation, moved, atol=1e-6)
    assert superposed_rmsd(coords, moved) < 1e-6


def test_tm_score_identical_is_one_and_normalised_by_given_length():
    coords = _helix(60)
    result = tm_score_fixed(coords, _rotate(coords), 60)
    assert result.tm_score == pytest.approx(1.0, abs=1e-4)
    assert result.rmsd < 1e-3
    half = tm_score_fixed(coords, _rotate(coords), 120)
    assert half.tm_score == pytest.approx(0.5, abs=1e-3)


def test_tm_score_finds_the_conserved_core_despite_a_moved_segment():
    coords = _helix(80)
    target = _rotate(coords)
    target[60:] += np.array([25.0, 0.0, 0.0])  # last quarter displaced
    result = tm_score_fixed(coords, target, 80)
    assert 0.65 < result.tm_score < 0.8  # ~60/80 residues superpose exactly
    assert result.rmsd > 5


def test_tm_d0_floor():
    assert tm_d0(10) == 0.5
    assert tm_d0(100) == pytest.approx(1.24 * 85 ** (1 / 3) - 1.8)


# --------------------------------------------------------------------------- fixtures/parse
def test_fixture_model_parses_and_matches_sequence():
    coords = read_ca(DATA / "fli1erg_rank_001.pdb.gz", FLI1_ERG)
    assert len(coords) == len(FLI1_ERG) == 150


def test_model_with_other_sequence_fails_closed():
    with pytest.raises(PredictionFailure, match="parse_error"):
        read_ca(DATA / "fli1erg_rank_001.pdb.gz", FLI1_ERG[:-1] + "A")


def _colabfold_out(tmp_path: Path, job: str = "job_x") -> Path:
    out = tmp_path / "out"
    out.mkdir()
    with gzip.open(DATA / "fli1erg_rank_001.pdb.gz", "rb") as src:
        (out / f"{job}_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_000.pdb").write_bytes(
            src.read()
        )
    with gzip.open(DATA / "fli1erg_scores_rank_001.json.gz", "rt") as handle:
        scores = json.load(handle)
    path = out / f"{job}_scores_rank_001_alphafold2_ptm_model_1_seed_000.json"
    path.write_text(json.dumps(scores))
    return out


def test_validate_outputs_accepts_a_complete_prediction(tmp_path):
    out = _colabfold_out(tmp_path)
    _, _, scores = validate_outputs(out, "job_x", FLI1_ERG)
    assert len(scores["plddt"]) == 150 and 0 < scores["ptm"] <= 1


@pytest.mark.parametrize(
    "damage, reason",
    [
        ("no_model", "missing_model_file"),
        ("short_plddt", "parse_error"),
        ("nan_pae", "parse_error"),
        ("bad_json", "parse_error"),
    ],
)
def test_validate_outputs_fails_closed(tmp_path, damage, reason):
    out = _colabfold_out(tmp_path)
    score = next(out.glob("*_scores_*.json"))
    if damage == "no_model":
        next(out.glob("*.pdb")).unlink()
    elif damage == "bad_json":
        score.write_text("{")
    else:
        data = json.loads(score.read_text())
        if damage == "short_plddt":
            data["plddt"] = data["plddt"][:-1]
        else:
            data["pae"][3][4] = float("nan")
        score.write_text(json.dumps(data))
    with pytest.raises(PredictionFailure, match=reason):
        validate_outputs(out, "job_x", FLI1_ERG)


def test_a3m_query():
    assert a3m_query("#150\t1\n>101\nACDE\n>hit\nAC-E\n") == "ACDE"
    assert a3m_query("") is None


# --------------------------------------------------------------------------- predictor
def _fake_colabfold(tmp_path: Path, *, exit_code: int = 0, damage: bool = False) -> list[str]:
    """A stand-in colabfold_batch that copies the fixture outputs (no ColabFold needed)."""
    script = tmp_path / "fake_colabfold.py"
    script.write_text(
        "import sys, gzip, json, pathlib\n"
        f"data = pathlib.Path({str(DATA)!r})\n"
        "inp, out = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2])\n"
        "out.mkdir(parents=True, exist_ok=True)\n"
        "job = next(inp.glob('*.a3m')).stem\n"
        f"if {exit_code}: sys.exit({exit_code})\n"
        "pdb = gzip.open(data / 'fli1erg_rank_001.pdb.gz').read()\n"
        "(out / f'{job}_unrelaxed_rank_001_m.pdb').write_bytes(pdb)\n"
        "scores = json.load(gzip.open(data / 'fli1erg_scores_rank_001.json.gz', 'rt'))\n"
        f"if {damage}: scores['plddt'] = scores['plddt'][:10]\n"
        "(out / f'{job}_scores_rank_001_m.json').write_text(json.dumps(scores))\n"
    )
    import sys

    return [sys.executable, str(script)]


def _predictor(tmp_path: Path, **kwargs) -> tuple[Predictor, MsaCache]:
    config = PredictorConfig(
        colabfold_command=_fake_colabfold(tmp_path, **kwargs), min_free_disk_gb=0.0
    )
    msas = MsaCache(tmp_path / "work", config)
    msas.path(FLI1_ERG).write_text(f">101\n{FLI1_ERG}\n")
    return Predictor(tmp_path / "work" / "predictions", config, msas), msas


def test_predictor_caches_only_validated_results(tmp_path):
    predictor, _ = _predictor(tmp_path)
    first = predictor.predict("job_a", FLI1_ERG, log=lambda _: None)
    assert first.ptm > 0 and len(first.pae) == 150
    assert (predictor.job_dir("job_a") / "result.json").exists()
    assert not (predictor.job_dir("job_a") / "colabfold").exists()
    # Second call reuses the checkpoint (the command is never run again).
    predictor.config.colabfold_command = ["false"]
    again = predictor.predict("job_a", FLI1_ERG, log=lambda _: None)
    assert again.plddt == first.plddt


@pytest.mark.parametrize("kwargs, reason", [({"exit_code": 3}, "exit status 3"),
                                            ({"damage": True}, "parse_error")])  # fmt: skip
def test_predictor_failure_is_recorded_but_never_cached(tmp_path, kwargs, reason):
    predictor, _ = _predictor(tmp_path, **kwargs)
    with pytest.raises(PredictionFailure, match=reason):
        predictor.predict("job_b", FLI1_ERG, log=lambda _: None)
    assert not (predictor.job_dir("job_b") / "result.json").exists()
    assert reason in json.loads((predictor.job_dir("job_b") / "failure.json").read_text())["reason"]
    assert predictor.cached("job_b", FLI1_ERG) is None


def test_predictor_without_msa_fails_closed(tmp_path):
    predictor, msas = _predictor(tmp_path)
    msas.path(FLI1_ERG).unlink()
    with pytest.raises(PredictionFailure, match="msa_missing"):
        predictor.predict("job_c", FLI1_ERG, log=lambda _: None)


def test_corrupt_checkpoint_is_discarded(tmp_path):
    predictor, _ = _predictor(tmp_path)
    predictor.predict("job_d", FLI1_ERG, log=lambda _: None)
    (predictor.job_dir("job_d") / "scores.json.gz").write_bytes(b"junk")
    assert predictor.cached("job_d", FLI1_ERG) is None


def test_msa_server_error_is_not_cached(tmp_path):
    config = PredictorConfig(msa_pause_seconds=0, min_free_disk_gb=0.0)
    msas = MsaCache(tmp_path, config)

    def failing_worker(command: list[str]) -> int:
        Path(command[-1]).write_text(json.dumps({"ok": False, "error": "RATELIMIT"}))
        return 1

    errors = msas.fetch(["ACDEFGHIK"], log=lambda _: None, runner=failing_worker)
    assert errors["ACDEFGHIK"].startswith("msa_server_error")
    assert msas.get("ACDEFGHIK") is None


def test_msa_batches_and_caches(tmp_path):
    config = PredictorConfig(msa_pause_seconds=0, msa_batch_size=2, min_free_disk_gb=0.0)
    msas = MsaCache(tmp_path, config)
    calls: list[int] = []

    def worker(command: list[str]) -> int:
        request = json.loads(Path(command[-2]).read_text())
        calls.append(len(request["queries"]))
        a3m = {q["key"]: f">101\n{q['sequence']}\n" for q in request["queries"]}
        Path(command[-1]).write_text(json.dumps({"ok": True, "a3m": a3m}))
        return 0

    sequences = ["AAAA", "CCCCC", "DDDDDD"]
    assert msas.fetch(sequences, log=lambda _: None, runner=worker) == {}
    assert calls == [2, 1]
    assert msas.fetch(sequences, log=lambda _: None, runner=worker) == {}
    assert calls == [2, 1]  # cached: no new ticket


def test_disk_guard(tmp_path):
    with pytest.raises(DiskSpaceError):
        check_disk(tmp_path, 1e9)


def test_settings_id_changes_with_prediction_settings_only():
    base = PredictorConfig()
    assert PredictorConfig(num_recycle=3).settings_id() != base.settings_id()
    assert PredictorConfig(colabfold_command=["x"]).settings_id() == base.settings_id()
    assert "--templates" not in base.batch_args()
    with pytest.raises(ValueError):
        PredictorConfig(backend="esmfold").batch_args()


def test_config_rejects_unknown_fields(tmp_path):
    path = tmp_path / "c.json"
    path.write_text(json.dumps({"num_recycle": 1, "colabfold_command": ["ssh", "gpu", "cf"]}))
    assert PredictorConfig.load(path).num_recycle == 1
    path.write_text(json.dumps({"predictor": "esmfold"}))
    with pytest.raises(ValueError):
        PredictorConfig.load(path)


# --------------------------------------------------------------------------- records/metrics
_PARENT = "MSTEKLLRAGHPWVDEQNTYFKAILGSRDCPEMHWNKAVLTGEFSRQDYPLIAKTEHMNGWSA"
_PARTNER = "MSTDKLLRAGHPWIDEQNTYFRAILGSKDCPEMHWNKSVLTGEFARQDYPLLAKSEHMNGWTA"


def _record(five_last: int = 30, three_first: int = 31) -> ChimeraRecord:
    protein = _PARENT[:five_last] + _PARTNER[three_first - 1 :]
    return ChimeraRecord(
        record_id="chim_test",
        pfam_family="PF00000",
        junction_class="test",
        class_reason="",
        five=ParentDomain("GENEA", "TXA", _PARENT, 5, 60, five_last),
        three=ParentDomain("GENEB", "TXB", _PARTNER, 5, 60, three_first),
        hybrid_residue="",
        chimeric_protein=protein,
        domain_start=5,
        domain_end=five_last + (60 - three_first + 1),
        junction_index=five_last,
        pilot_indel_residues=None,
        pilot_indel_sequence="",
        pilot_alignment_parent="",
        pilot_alignment_identity=None,
        sequence_id="seq_test",
    )


def _write_pdb(path: Path, sequence: str, coords: np.ndarray) -> None:
    three = {v: k for k, v in __import__("cfh.domain_structure.predictor").domain_structure
             .predictor.THREE_TO_ONE.items()}  # fmt: skip
    lines = []
    for i, (letter, xyz) in enumerate(zip(sequence, coords, strict=True), start=1):
        lines.append(
            f"ATOM  {i:5d}  CA  {three[letter]} A{i:4d}    "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00 90.00           C"
        )
    with gzip.open(path, "wt") as handle:
        handle.write("\n".join(lines + ["END"]) + "\n")


def _prediction(
    tmp_path: Path, name: str, sequence: str, coords: np.ndarray, plddt: float, pae: float
) -> Prediction:
    path = tmp_path / f"{name}.pdb.gz"
    _write_pdb(path, sequence, coords)
    n = len(sequence)
    return Prediction(
        job_id=name,
        sequence=sequence,
        plddt=[plddt] * n,
        pae=[[pae] * n for _ in range(n)],
        ptm=0.8,
        model_name="m",
        pdb_path=path,
        wall_seconds=1.0,
    )


def _jobs(record: ChimeraRecord):
    row = {"record_id": record.record_id, "sequence_id": record.sequence_id, "tier": "uncertain"}
    record_jobs, jobs = build_jobs([row], [record])
    return record_jobs[0], jobs


def _predictions(tmp_path, record, chimera_plddt=90.0, chimera_pae=3.0, chimera_coords=None):
    rj, _ = _jobs(record)
    chim_seq = record.chimeric_domain_flanked
    p5 = record.five.protein
    p3 = record.three.protein
    return rj, {
        rj.chimera: _prediction(
            tmp_path,
            "c",
            chim_seq,
            chimera_coords if chimera_coords is not None else _helix(len(chim_seq)),
            chimera_plddt,
            chimera_pae,
        ),
        rj.parent5: _prediction(tmp_path, "p5", p5, _rotate(_helix(len(p5))), 90.0, 3.0),
        rj.parent3: _prediction(tmp_path, "p3", p3, _helix(len(p3)), 88.0, 3.0),
    }


def test_build_jobs_dedupes_sequences_and_keeps_windows():
    record = _record()
    rj, jobs = _jobs(record)
    assert len(jobs) == 3
    assert rj.chimera_start == 1 and rj.parent5_start == 1
    second = _record()
    second.record_id = "chim_other"
    rows = [{"record_id": r.record_id, "sequence_id": r.sequence_id} for r in (record, second)]
    _, jobs = build_jobs(rows, [record, second])
    assert len(jobs) == 3  # same three sequences, predicted once


def test_build_jobs_rejects_sequence_id_mismatch():
    record = _record()
    with pytest.raises(ValueError):
        build_jobs([{"record_id": record.record_id, "sequence_id": "seq_x"}], [record])


def test_windows_and_halves():
    record = _record(five_last=30, three_first=26)  # 5-residue tandem overlap
    assert junction_window(record, [31, 32, 33, 34, 35])[0] == 21
    assert junction_window(record, [31, 32, 33, 34, 35])[-1] == 45
    five, three = chimera_halves(record)
    assert five[-1] == 30 and three[0] == 31
    a, b = parent_halves(record.three, "three")
    assert a[-1] == 25 and b[0] == 26


def test_record_metrics_identical_fold_is_supported(tmp_path):
    record = _record()
    rj, predictions = _predictions(tmp_path, record)
    row = record_metrics(rj, predictions, {}, {})
    assert row["vs_5p_status"] == "ok" and row["vs_5p_tm"] == pytest.approx(1.0, abs=1e-3)
    assert row["vs_3p_tm"] == pytest.approx(1.0, abs=1e-3)
    assert row["vs_5p_local_rmsd"] < 1e-2
    assert row["delta_domain_plddt_vs_better_parent"] == 0.0
    assert row["chimera_interhalf_pae"] == 3.0
    assert row["vs_afdb_5p_status"].startswith("not_assessed")
    flag, rules, _ = assign_flag(row)
    assert (flag, rules) == (FOLD_SUPPORTED, ["S1_all_pass"])


def test_record_metrics_flags_low_confidence_and_unpacked_halves(tmp_path):
    record = _record()
    rj, predictions = _predictions(tmp_path, record, chimera_plddt=55.0, chimera_pae=20.0)
    row = record_metrics(rj, predictions, {}, {})
    flag, rules, _ = assign_flag(row)
    assert flag == FOLD_DOUBTFUL
    assert {"D1_domain_plddt_low", "D2_domain_plddt_drop", "D5_halves_not_packed"} <= set(rules)


def test_record_metrics_flags_a_different_fold(tmp_path):
    record = _record()
    n = len(record.chimeric_domain_flanked)
    rng = np.random.default_rng(0)
    scrambled = rng.normal(scale=15.0, size=(n, 3))
    rj, predictions = _predictions(tmp_path, record, chimera_coords=scrambled)
    row = record_metrics(rj, predictions, {}, {})
    assert row["vs_5p_tm"] < 0.5
    assert "D6_fold_differs" in assign_flag(row)[1]


def test_missing_prediction_is_not_assessed(tmp_path):
    record = _record()
    rj, predictions = _predictions(tmp_path, record)
    del predictions[rj.parent3]
    row = record_metrics(rj, predictions, {rj.parent3: "msa_server_error: 503"}, {})
    flag, rules, missing = assign_flag(row)
    assert flag == NOT_ASSESSED and rules == ["N0_prediction_missing"]
    assert "msa_server_error" in missing[0]


def test_parents_not_confident_is_not_assessed():
    row = {
        "chimera_prediction_status": "ok",
        "parent5_prediction_status": "ok",
        "parent3_prediction_status": "ok",
        "parent5_domain_plddt": 45.0,
        "parent3_domain_plddt": 60.0,
        "chimera_domain_plddt": 40.0,
    }
    assert assign_flag(row)[:2] == (NOT_ASSESSED, ["N1_parents_not_confident"])


def test_missing_metric_blocks_supported_but_not_doubtful():
    row = {
        "chimera_prediction_status": "ok",
        "parent5_prediction_status": "ok",
        "parent3_prediction_status": "ok",
        "parent5_domain_plddt": 90.0,
        "parent3_domain_plddt": 90.0,
        "chimera_domain_plddt": 89.0,
    }
    assert assign_flag(row)[0] == NOT_ASSESSED
    row["chimera_domain_plddt"] = 60.0
    assert assign_flag(row)[0] == FOLD_DOUBTFUL


def test_every_rule_is_documented():
    assert all(rule in RULES for rule in ("D1_domain_plddt_low", "N0_prediction_missing"))
    doc = (Path(__file__).parents[1] / "docs" / "domain_structure.md").read_text()
    for rule in RULES:
        assert rule in doc


def test_afdb_comparison_uses_the_mapped_residues(tmp_path):
    record = _record()
    rj, predictions = _predictions(tmp_path, record)
    model = AfdbModel(
        accession="P00000",
        version="4",
        sequence=_PARENT,
        ca=tuple(tuple(x) for x in _rotate(_helix(len(_PARENT))).tolist()),
        plddt=tuple([80.0] * len(_PARENT)),
    )
    row = record_metrics(rj, predictions, {}, {"five": model, "three": "no AFDB model"})
    assert row["vs_afdb_5p_status"] == "ok"
    assert row["vs_afdb_5p_tm"] == pytest.approx(1.0, abs=1e-3)
    assert row["parent5_vs_afdb_tm"] == pytest.approx(1.0, abs=1e-3)
    assert row["vs_afdb_3p_status"] == "not_assessed: no AFDB model"


def test_assess_and_counts(tmp_path):
    record = _record()
    record.cohorts = ["cohortA"]
    rj, predictions = _predictions(tmp_path, record)
    rows = assess([rj], predictions, {}, {})
    assert rows[0]["pr1_tier"] == "uncertain"
    counts = flag_counts(rows, [record])
    assert counts["by_flag"] == {FOLD_SUPPORTED: 1}
    assert counts["by_cohort"] == {"cohortA": {FOLD_SUPPORTED: 1}}


def test_window_pae_helpers(tmp_path):
    prediction = _prediction(tmp_path, "w", "ACDEF", _helix(5), 70.0, 2.0)
    prediction.pae[0][4] = 10.0
    window = Window(prediction, 11)
    assert window.inside([10, 11, 15, 16]) == [11, 15]
    assert window.mean_pae([11], [15]) == pytest.approx(6.0)


def test_mapped_tm_is_normalised_by_the_parent_domain(tmp_path):
    """An insertion (unmapped chimera residues) does not lower ``tm``; a deletion does."""
    coords = _helix(40)
    window = Window(_prediction(tmp_path, "n", "A" * 40, coords, 90.0, 3.0), 1)
    pairs = [(c, c) for c in range(1, 21)]  # only the first 20 residues map
    insertion = compare(window, _rotate(coords), 1, pairs, [20, 21], 20, 40)
    assert insertion["tm"] == pytest.approx(1.0, abs=1e-3)
    assert insertion["tm_chim_norm"] < 0.55
    deletion = compare(window, _rotate(coords), 1, pairs, [20, 21], 40, 20)
    assert deletion["tm"] < 0.55
    assert deletion["tm_chim_norm"] == pytest.approx(1.0, abs=1e-3)
    assert compare(window, coords, 1, pairs[:5], [3], 20, 20)["status"].startswith("mapping_failed")


def test_export_models_only_selected_jobs_and_rounds_pae(tmp_path):
    prediction = _prediction(tmp_path, "e", "ACDEF", _helix(5), 70.0, 2.4)
    other = _prediction(tmp_path, "f", "GHIKL", _helix(5), 70.0, 2.6)
    jobs = [Job("job_e", "ACDEF", ["chimera:X::Y"]), Job("job_f", "GHIKL", ["parent5:X"])]
    index = export_models(
        jobs, {"job_e": prediction, "job_f": other}, tmp_path / "out", {"job_e"}, 0
    )
    assert index["job_e"]["exported"] and not index["job_f"]["exported"]
    assert sorted(p.name for p in (tmp_path / "out").iterdir()) == [
        "job_e.pdb.gz",
        "job_e.scores.json.gz",
    ]
    with gzip.open(tmp_path / "out" / "job_e.scores.json.gz", "rt") as handle:
        scores = json.load(handle)
    assert scores["pae"][0][1] == 2 and scores["pae_decimals"] == 0


# --------------------------------------------------------------------------- live
@pytest.mark.network
def test_live_colabfold_msa_server(tmp_path):
    config = PredictorConfig(min_free_disk_gb=0.0)
    if not Path(config.msa_python).exists() or shutil.which(config.msa_python) is None:
        pytest.skip("ColabFold env not installed")
    msas = MsaCache(tmp_path, config)
    assert msas.fetch([FLI1_ERG[:60]], log=lambda _: None) == {}
    assert msas.get(FLI1_ERG[:60]) is not None
