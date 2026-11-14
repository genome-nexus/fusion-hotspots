"""PLM scores for chimeric domains: fail-closed rules, length-robust metrics, no tier changes."""

from __future__ import annotations

import json
import math
import sys
from dataclasses import asdict, replace
from pathlib import Path

import pytest

from cfh.domain_function.record import ChimeraRecord, ParentDomain
from cfh.domain_plm.cache import ModelSpec, ScoreCache, sequence_sha, valid_logp
from cfh.domain_plm.compare import add_ranks, disagreements
from cfh.domain_plm.metrics import METRIC_COLUMNS, record_metrics
from cfh.domain_plm.run import collect, metrics_row, score_in_order, scoring_order
from cfh.domain_plm.scoring import (
    ScoringError,
    WorkerScorer,
    check_against_metrics,
    plan_record,
    score_sequences,
)
from cfh.domain_plm.worker import MAX_RESIDUES, batch_rows, context_windows
from cfh.domain_plm.write import COLUMNS, counts, report_tables

SPEC = ModelSpec("esm2_test", "2.0.0", "a" * 64)
PARENT = "MSTEKLLRAGHPWVDEQNTYFKAILGSRDCPEMHWNKAVLTGEFSRQDYPLIAKTEHMNGWSA"
PARTNER = "MSTDKLLRAGHPWIDEQNTYFRAILGSKDCPEMHWNKSVLTGEFARQDYPLLAKSEHMNGWTA"


def _record(
    five: str = PARENT,
    three: str = PARTNER,
    *,
    five_last: int = 30,
    three_first: int = 31,
    hybrid: str = "",
    record_id: str = "chim_test",
) -> ChimeraRecord:
    protein = five[:five_last] + hybrid + three[three_first - 1 :]
    junction = five_last + len(hybrid)
    return ChimeraRecord(
        record_id=record_id,
        pfam_family="PF00000",
        junction_class="test",
        class_reason="",
        five=ParentDomain("GENEA", "TXA", five, 5, 60, five_last),
        three=ParentDomain("GENEB", "TXB", three, 5, 60, three_first),
        hybrid_residue=hybrid,
        chimeric_protein=protein,
        domain_start=5,
        domain_end=junction + (60 - three_first + 1),
        junction_index=junction,
        pilot_indel_residues=None,
        pilot_indel_sequence="",
        pilot_alignment_parent="",
        pilot_alignment_identity=None,
    )


class FakeScorer:
    """Deterministic scores; ``value`` overrides every score (e.g. NaN)."""

    def __init__(self, value: float | None = None, fail: bool = False, drop: int = 0):
        self.spec = SPEC
        self.value = value
        self.fail = fail
        self.drop = drop
        self.calls = 0

    def score(self, sequence: str, positions: list[int]) -> dict[int, float]:
        self.calls += 1
        if self.fail:
            raise ScoringError("model crashed")
        out = {
            p: (self.value if self.value is not None else -1.0 - (ord(sequence[p - 1]) % 7) / 3)
            for p in positions
        }
        for p in positions[: self.drop]:
            del out[p]
        return out


def _scored(tmp_path: Path, record: ChimeraRecord, scorer: FakeScorer | None = None):
    plan = plan_record(record)
    cache = ScoreCache(tmp_path, SPEC)
    work = [(side.sequence, side.needed()) for side in plan.sides().values()]
    failures = score_sequences(work, scorer or FakeScorer(), cache, log=lambda _: None)
    return plan, cache, failures


# ------------------------------------------------------------------------------- cache
def test_cache_round_trip_and_keying(tmp_path):
    cache = ScoreCache(tmp_path, SPEC)
    cache.put("ACD", {1: -0.5, 3: -2.0})
    assert cache.get("ACD") == {1: -0.5, 3: -2.0}
    assert cache.get("ACE") == {}
    sha = sequence_sha("ACD")
    assert cache.path("ACD").name == f"{sha}.json"
    assert cache.path("ACD").parent.parent.name == SPEC.spec_id()


@pytest.mark.parametrize(
    "other",
    [
        ModelSpec("esm2_other", "2.0.0", "a" * 64),
        ModelSpec("esm2_test", "2.1.0", "a" * 64),
        ModelSpec("esm2_test", "2.0.0", "b" * 64),
    ],
)
def test_cache_from_another_model_name_version_or_weights_is_not_reused(tmp_path, other):
    ScoreCache(tmp_path, SPEC).put("ACD", {1: -0.5})
    assert other.spec_id() != SPEC.spec_id()
    assert ScoreCache(tmp_path, other).get("ACD") == {}


def _tamper(cache: ScoreCache, seq: str, **changes):
    path = cache.path(seq)
    entry = json.loads(path.read_text())
    entry.update(changes)
    path.write_text(json.dumps(entry))


@pytest.mark.parametrize(
    "changes",
    [
        {"spec": asdict(ModelSpec("esm2_test", "1.0", "a" * 64))},
        {"sequence": "ACE"},
        {"sequence_sha256": "0" * 64},
        {"logp": {"1": "NaN"}},
        {"logp": {"1": 0.7}},
        {"logp": {"9": -0.1}},
        {"logp": {"x": -0.1}},
        {"logp": {"1": True}},
        {"logp": [None, -0.5]},
        {"logp": "invalid"},
        {"logp": {"01": -0.5}},
    ],
)
def test_stale_or_corrupt_cache_entry_reads_as_absent(tmp_path, changes):
    cache = ScoreCache(tmp_path, SPEC)
    cache.put("ACD", {1: -0.5})
    _tamper(cache, "ACD", **changes)
    assert cache.get("ACD") == {}


def test_cache_with_non_finite_json_reads_as_absent(tmp_path):
    cache = ScoreCache(tmp_path, SPEC)
    cache.put("ACD", {1: -0.5})
    path = cache.path("ACD")
    path.write_text(path.read_text().replace("-0.5", "NaN"))
    assert cache.get("ACD") == {}
    path.write_text("{not json")
    assert cache.get("ACD") == {}


@pytest.mark.parametrize("entry", [None, [], "invalid", 5])
def test_cache_with_wrong_json_shape_reads_as_absent(tmp_path, entry):
    cache = ScoreCache(tmp_path, SPEC)
    cache.put("ACD", {1: -0.5})
    cache.path("ACD").write_text(json.dumps(entry))
    assert cache.get("ACD") == {}


@pytest.mark.parametrize("bad", [math.nan, math.inf, -math.inf, 0.3])
def test_cache_refuses_to_store_invalid_scores(tmp_path, bad):
    cache = ScoreCache(tmp_path, SPEC)
    with pytest.raises(ValueError):
        cache.put("ACD", {1: bad})
    assert not cache.path("ACD").exists()


def test_valid_logp():
    assert valid_logp(-1.0) and valid_logp(0.0) and valid_logp(-3)
    assert not any(valid_logp(v) for v in (math.nan, math.inf, -math.inf, 0.01, "1", None, True))


# ----------------------------------------------------------------------------- scoring
@pytest.mark.parametrize(
    "scorer",
    [
        FakeScorer(value=math.nan),
        FakeScorer(value=math.inf),
        FakeScorer(value=-math.inf),
        FakeScorer(value=0.5),
    ],
)
def test_non_finite_or_positive_scores_fail_the_sequence_and_are_not_cached(tmp_path, scorer):
    plan, cache, failures = _scored(tmp_path, _record(), scorer)
    assert set(failures) == {s.sequence for s in plan.sides().values()}
    assert all(cache.get(s.sequence) == {} for s in plan.sides().values())
    row = metrics_row(plan, cache, failures, "", "")
    assert row["plm_status"] == "not_assessed"
    assert "scoring failed" in row["plm_not_assessed_reason"]
    assert all(row[c] == "" for c in METRIC_COLUMNS[2:])


def test_failed_model_call_is_not_assessed(tmp_path):
    plan, cache, failures = _scored(tmp_path, _record(), FakeScorer(fail=True))
    assert len(failures) == 3
    assert metrics_row(plan, cache, failures, "", "")["plm_status"] == "not_assessed"


def test_reply_missing_positions_is_a_failure(tmp_path):
    plan, cache, failures = _scored(tmp_path, _record(), FakeScorer(drop=1))
    assert len(failures) == 3
    assert metrics_row(plan, cache, failures, "", "")["plm_status"] == "not_assessed"


def test_reply_with_extra_positions_is_not_assessed_or_cached(tmp_path):
    class ExtraScorer(FakeScorer):
        def score(self, sequence, positions):
            return {**super().score(sequence, positions), len(sequence) + 1: -1.0}

    plan, cache, failures = _scored(tmp_path, _record(), ExtraScorer())
    assert len(failures) == 3
    assert all(cache.get(s.sequence) == {} for s in plan.sides().values())
    assert metrics_row(plan, cache, failures, "", "")["plm_status"] == "not_assessed"


@pytest.mark.parametrize(
    "changes",
    [
        {"spec": asdict(replace(SPEC, model="other"))},
        {"spec": asdict(replace(SPEC, weights_sha256="b" * 64))},
        {"sequence": "OTHER"},
        {"logp": {"1": math.nan}},
        {"logp": {"1": math.inf}},
        {"logp": {"1": -math.inf}},
        {"logp": []},
    ],
)
def test_stale_or_invalid_required_cache_makes_record_not_assessed(tmp_path, changes):
    plan, cache, _ = _scored(tmp_path, _record())
    _tamper(cache, plan.chimera.sequence, **changes)
    row = metrics_row(plan, cache, {}, "", "")
    assert row["plm_status"] == "not_assessed"
    assert all(row[c] == "" for c in METRIC_COLUMNS[2:])
    add_ranks([row])
    assert row["plm_rank"] == "not_assessed"
    assert disagreements(row) == []


def test_missing_model_makes_every_record_not_assessed(tmp_path):
    plan, cache, failures = _scored(tmp_path, _record())
    row = metrics_row(plan, cache, failures, "", "PLM model unavailable: OSError")
    assert row["plm_status"] == "not_assessed"
    row = metrics_row(plan, None, {}, "", "PLM model unavailable: OSError")
    assert row["plm_status"] == "not_assessed"


def test_unscored_or_partially_cached_record_is_not_assessed(tmp_path):
    plan = plan_record(_record())
    cache = ScoreCache(tmp_path, SPEC)
    assert metrics_row(plan, cache, {}, "", "")["plm_status"] == "not_assessed"
    chim = plan.chimera
    cache.put(chim.sequence, {chim.needed()[0]: -1.0})
    row = metrics_row(plan, cache, {}, "", "")
    assert row["plm_status"] == "not_assessed"
    assert "not scored" in row["plm_not_assessed_reason"]


def test_cache_is_reused_without_calling_the_model_again(tmp_path):
    record = _record()
    scorer = FakeScorer()
    _scored(tmp_path, record, scorer)
    calls = scorer.calls
    _scored(tmp_path, record, scorer)
    assert scorer.calls == calls


def test_fully_scored_record_has_finite_metrics(tmp_path):
    plan, cache, failures = _scored(tmp_path, _record(hybrid="W"))
    row = metrics_row(plan, cache, failures, "", "")
    assert row["plm_status"] == "scored"
    for column in METRIC_COLUMNS[2:]:
        value = row[column]
        if isinstance(value, float):
            assert math.isfinite(value)
    assert row["n_paired_5p"] + row["n_paired_3p"] == len(_record(hybrid="W").chimeric_domain) - 1


def test_metrics_with_a_non_finite_score_are_not_assessed(tmp_path):
    plan, cache, failures = _scored(tmp_path, _record())
    scores, _ = collect(plan, cache)
    first = next(iter(scores["parent3"]))
    scores["parent3"][first] = math.nan
    assert record_metrics(plan, scores, {})["plm_status"] == "not_assessed"


# ----------------------------------------------------------------------- record checks
def test_non_canonical_or_too_long_sequences_are_not_assessed(tmp_path):
    odd = PARENT[:20] + "X" + PARENT[21:]
    plan = plan_record(_record(five=odd))
    assert "non-canonical" in plan.problem
    record = _record(three=PARTNER + "A" * 4100)
    record.three = replace(record.three, domain_end=4150)
    plan_long = plan_record(record)
    assert "exceeds limit" in plan_long.problem
    assert metrics_row(plan, ScoreCache(tmp_path, SPEC), {}, "", "")["plm_status"] == (
        "not_assessed"
    )


def test_context_windows_cover_every_position_inside_the_model_limit():
    assert context_windows(300, [5, 1, 300]) == {1: [1, 5, 300]}
    groups = context_windows(1500, list(range(1, 1501)))
    assert sorted(p for g in groups.values() for p in g) == list(range(1, 1501))
    for start, group in groups.items():
        assert 1 <= start and start + MAX_RESIDUES - 1 <= 1500
        assert all(start <= p <= start + MAX_RESIDUES - 1 for p in group)
    assert context_windows(1500, [1])[1] == [1]
    assert context_windows(1500, [750]) == {750 - MAX_RESIDUES // 2: [750]}


@pytest.mark.parametrize("role", ["chimera", "parent5", "parent3"])
def test_truncated_domain_sequence_is_not_assessed(tmp_path, role):
    record = _record()
    if role == "chimera":
        record.chimeric_protein = record.chimeric_protein[:40]
    else:
        parent = record.five if role == "parent5" else record.three
        parent.protein = parent.protein[:40]
    plan = plan_record(record)
    assert f"{role}: domain span is outside" in plan.problem
    assert metrics_row(plan, None, {}, "", "")["plm_status"] == "not_assessed"


def test_record_longer_than_the_model_context_is_scored_with_a_note(tmp_path):
    record = _record(three=PARTNER + "A" * 1100)
    record.three = replace(record.three, domain_end=1150)
    plan, cache, failures = _scored(tmp_path, record)
    assert plan.problem == "" and plan.windowed_roles() == ["parent3"]
    row = metrics_row(plan, cache, failures, "", "")
    assert row["plm_status"] == "scored"
    assert row["plm_note"] == f"parent3_windowed_context_{len(plan.parent3.sequence)}aa"


def _metrics_row_for(record: ChimeraRecord) -> dict[str, str]:
    return {
        "sequence_id": record.sequence_id,
        "fusion": record.fusion,
        "chimeric_domain_span": f"{record.domain_start}-{record.domain_end}",
        "five_prime_domain": f"{record.five.domain_start}-{record.five.domain_end}",
        "three_prime_domain": f"{record.three.domain_start}-{record.three.domain_end}",
        "junction_index": str(record.junction_index),
        "chimeric_domain_length": str(len(record.chimeric_domain)),
    }


def test_sequence_mismatch_against_pr155_metrics_is_not_assessed(tmp_path):
    record = _record()
    plan, cache, failures = _scored(tmp_path, record)
    good = _metrics_row_for(record)
    assert check_against_metrics(plan, good) == ""
    assert check_against_metrics(plan, None) != ""
    bad = {**good, "chimeric_domain_length": "999"}
    mismatch = check_against_metrics(plan, bad)
    assert "mismatch" in mismatch
    assert metrics_row(plan, cache, failures, mismatch, "")["plm_status"] == "not_assessed"


# ------------------------------------------------------------------ length robustness
def test_paired_context_shift_ignores_length_change_from_a_deletion(tmp_path):
    # identical parents; the 3' side restarts 6 residues later (a 5-residue deletion)
    record = _record(five=PARENT, three=PARENT, five_last=30, three_first=36)
    assert len(record.chimeric_domain) == 56 - 5  # 56-residue parent domain, 5 deleted

    class ContextFree(FakeScorer):
        def score(self, sequence, positions):  # score depends on the residue only
            return {p: -1.0 - (ord(sequence[p - 1]) % 5) / 4 for p in positions}

    plan, cache, failures = _scored(tmp_path, record, ContextFree())
    row = metrics_row(plan, cache, failures, "", "")
    assert row["plm_status"] == "scored"
    assert row["context_shift_5p"] == pytest.approx(0.0)
    assert row["context_shift_3p"] == pytest.approx(0.0)
    assert row["n_paired_5p"] == 30 - 5 + 1
    assert row["n_paired_3p"] == 60 - 36 + 1


# ------------------------------------------------------- ranks, disagreements, tiers
def _row(source="pr155", status="scored", tier="likely_intact", fold="fold_supported", **kw):
    row = {
        "source": source,
        "plm_status": status,
        "pr155_tier": tier,
        "pr160_fold_flag": fold,
        "delta_vs_better_parent": kw.get("delta", 0.0),
        "dip_excess_worst": kw.get("dip", 0.0),
    }
    return row


def test_not_assessed_never_gets_a_rank_or_takes_part_in_a_plm_disagreement():
    rows = [_row(delta=-float(i), dip=0.0) for i in range(20)]
    rows.append(_row(status="not_assessed", tier="likely_broken", fold="fold_doubtful"))
    rows.append(_row(delta=math.nan))
    rows.append(_row(delta=""))
    add_ranks(rows)
    for row in rows[-3:]:
        assert row["plm_rank"] == "not_assessed"
        assert not any(d.startswith("plm_") for d in disagreements(row))
    assert disagreements(rows[-3]) == []


def test_bottom_decile_and_disagreement_labels():
    rows = [_row(delta=-float(i)) for i in range(20)]
    add_ranks(rows)
    low = [r for r in rows if r["plm_rank"] == "bottom_decile"]
    assert {r["delta_vs_better_parent"] for r in low} == {-19.0, -18.0}
    assert "plm_low_but_tier_likely_intact" in disagreements(low[0])
    assert "plm_low_but_fold_supported" in disagreements(low[0])
    broken = {**rows[0], "pr155_tier": "likely_broken", "pr160_fold_flag": "fold_doubtful"}
    assert disagreements(broken) == [
        "plm_not_low_but_tier_likely_broken",
        "plm_not_low_but_fold_doubtful",
    ]


def test_extra_records_are_ranked_but_do_not_shape_the_reference():
    rows = [_row(delta=-float(i)) for i in range(10)]
    rows.append(_row(source="pcawg_tcga_dna", delta=-100.0))
    add_ranks(rows)
    assert rows[-1]["plm_rank"] == "bottom_decile"
    assert rows[-1]["rank_delta_vs_better_parent"] == 0.0
    assert [r["plm_rank"] for r in rows[:10]].count("bottom_decile") == 1


def test_plm_never_changes_pr155_tier_or_pr160_flag(tmp_path):
    plan, cache, failures = _scored(tmp_path, _record())
    metrics = metrics_row(plan, cache, failures, "", "")
    assert not {"tier", "pr155_tier", "structure_flag", "pr160_fold_flag"} & set(metrics)
    rows = [
        {**_row(tier=t, fold=f, delta=-float(i)), "disagreements": ""}
        for i, (t, f) in enumerate(
            [("likely_intact", "fold_supported"), ("likely_broken", "fold_doubtful")] * 10
        )
    ]
    before = [(r["pr155_tier"], r["pr160_fold_flag"]) for r in rows]
    add_ranks(rows)
    for row in rows:
        disagreements(row)
    assert [(r["pr155_tier"], r["pr160_fold_flag"]) for r in rows] == before


# ------------------------------------------------------------------- order and outputs
def test_scoring_order_puts_priority_records_first_and_is_deterministic():
    plans = [plan_record(_record(record_id=f"r{i}")) for i in range(8)]
    order = scoring_order(plans, ["r5", "r2", "missing"], [], {"r7"})
    assert [p.record.record_id for _, p in order[:3]] == ["r5", "r2", "r7"]
    assert [g for g, _ in order[:3]] == ["shortlist", "shortlist", "pcawg_dna"]
    assert len(order) == 8
    again = scoring_order(plans, ["r5", "r2"], [], {"r7"})
    assert [p.record.record_id for _, p in again] == [p.record.record_id for _, p in order]


def test_score_in_order_stops_at_the_deadline(tmp_path):
    plans = [plan_record(_record(record_id=f"r{i}")) for i in range(3)]
    order = scoring_order(plans, [], [], set())
    done, _ = score_in_order(order, FakeScorer(), ScoreCache(tmp_path, SPEC), lambda _: None, 0)
    assert done == 0


def test_counts_and_tables_render_not_assessed_rows():
    rows = []
    for status in ("scored", "not_assessed"):
        row = {c: "" for c in COLUMNS}
        row.update(
            source="pr155",
            scoring_group="shortlist",
            fusion="FLI1::ERG",
            plm_status=status,
            plm_not_assessed_reason="" if status == "scored" else "chimera: not scored",
            plm_rank="bottom_decile" if status == "scored" else "not_assessed",
            pr155_tier="likely_intact",
            pr160_fold_flag="fold_supported",
            disagreements="plm_low_but_tier_likely_intact" if status == "scored" else "",
        )
        rows.append(row)
    summary = counts(rows)
    assert summary["plm_status"] == {"not_assessed": 1, "scored": 1}
    assert summary["disagreements"] == {"plm_low_but_tier_likely_intact": 1}
    text = report_tables(rows, ("FLI1::ERG",))
    assert "not assessed: chimera: not scored" in text


# --------------------------------------------------------------------- worker client
FAKE_WORKER = """
import json, sys
mode = sys.argv[1]
if mode == "die":
    sys.exit(3)
if mode == "noheader":
    print("garbage", flush=True)
    sys.exit(0)
print(json.dumps({"header": {"model": mode, "esm_version": "x", "weights_sha256": "w"}}),
      flush=True)
for line in sys.stdin:
    req = json.loads(line)
    if mode == "notobject":
        print("[]", flush=True)
        continue
    if mode == "wrongid":
        print(json.dumps({"id": "wrong"}), flush=True)
        continue
    seq = req["sequence"]
    if mode == "wrongseq":
        seq = seq[::-1]
    if mode == "error":
        print(json.dumps({"id": req["id"], "error": "boom"}), flush=True)
        continue
    print(json.dumps({"id": req["id"], "sequence": seq,
                      "logp": {str(p): -1.0 for p in req["positions"]}}), flush=True)
"""


@pytest.fixture
def fake_worker(tmp_path):
    path = tmp_path / "fake_worker.py"
    path.write_text(FAKE_WORKER)
    return path


def test_worker_client_round_trip(fake_worker):
    scorer = WorkerScorer(sys.executable, "ok", script=fake_worker)
    assert scorer.spec == ModelSpec("ok", "x", "w")
    assert scorer.score("ACDE", [1, 3]) == {1: -1.0, 3: -1.0}
    scorer.close()


@pytest.mark.parametrize("mode", ["wrongseq", "error", "notobject", "wrongid"])
def test_worker_client_rejects_bad_replies(fake_worker, mode):
    scorer = WorkerScorer(sys.executable, mode, script=fake_worker)
    with pytest.raises(ScoringError):
        scorer.score("ACDE", [1])
    scorer.close()


@pytest.mark.parametrize("mode", ["die", "noheader"])
def test_worker_that_cannot_start_raises(fake_worker, mode):
    with pytest.raises(ScoringError):
        WorkerScorer(sys.executable, mode, script=fake_worker)


def test_batch_rows_bound_tokens_per_batch():
    assert batch_rows(100, 64) == 64
    assert batch_rows(609, 64) == 8192 // 609
    assert batch_rows(1024, 64) == 8
    assert batch_rows(20000, 64) == 1
