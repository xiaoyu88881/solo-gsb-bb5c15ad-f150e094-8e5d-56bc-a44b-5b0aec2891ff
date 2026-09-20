"""批量盲样分配：候选边筛选、全局一对一最优分配与破同分的单元 + API 端到端测试。"""
from fastapi.testclient import TestClient

from app.main import app
from app.batch_identify import resolve_batch_identify
from app.models import BatchIdentifyRequest

client = TestClient(app)


# ---------------------------------------------------------------- 构造工具

def _req(profiles, lanes, tolerance=10, max_unmatched=100):
    return BatchIdentifyRequest(
        profiles=[
            {"name": n, "predicted_bands": b} for n, b in profiles
        ],
        lanes=[
            {"name": n, "observed_bands": b} for n, b in lanes
        ],
        band_match_tolerance=tolerance,
        max_unmatched_bands=max_unmatched,
    )


def _assignment_map(res):
    return {item.lane: item.profile for item in res.assignments}


# ---------------------------------------------------------------- 纯算法

def test_batch_perfect_assignment_uses_evidence():
    # alpha 完美匹配 lane-1，beta 完美匹配 lane-3；lane-2 两条边都有未匹配
    req = _req(
        [("alpha", [230, 320, 450]), ("beta", [318, 460, 700]),
         ("gamma", [300, 600, 900])],
        [("lane-1", [231, 319, 455]), ("lane-2", [200, 322, 451]),
         ("lane-3", [317, 462, 698])],
        max_unmatched=2,
    )
    res = resolve_batch_identify(req)
    assert res.evaluated_edge_count == 9
    # gamma 三条边均未匹配 6 条 > 2 -> 不可用；alpha-lane-3、beta-lane-1 同理
    assert res.eligible_edge_count < 9
    assert _assignment_map(res) == {"lane-1": "alpha", "lane-3": "beta"}
    assert res.unassigned_lanes == ["lane-2"]
    assert res.unused_profiles == ["gamma"]
    assert res.assigned_lane_count == 2
    assert res.total_unmatched_count == 0
    # 逐项依据与复用的一对一匹配统计自洽
    for item in res.assignments:
        m = item.match
        assert item.unmatched_count == m.missing_count + m.extra_count
        assert m.matched_count == len(m.matches)
        assert item.description


def test_batch_threshold_removes_edges_and_grows_cardinality():
    # 阈值收紧时不可用边剔除、已分配数下降；放宽后更多泳道可分配
    profiles = [("A", [100]), ("B", [200]), ("C", [300])]
    lanes = [("x", [100]), ("y", [200, 500]), ("z", [600])]
    # 阈值 0：仅 A-x 完全配对
    req0 = _req(profiles, lanes, tolerance=0, max_unmatched=0)
    res0 = resolve_batch_identify(req0)
    assert res0.eligible_edge_count == 1
    assert _assignment_map(res0) == {"x": "A"}
    assert res0.unassigned_lanes == ["y", "z"]
    # 阈值 1：B-y（额外 1 条）可用 -> 2 个泳道
    res1 = resolve_batch_identify(req0.model_copy(update={"max_unmatched_bands": 1}))
    assert res1.assigned_lane_count == 2
    assert _assignment_map(res1) == {"x": "A", "y": "B"}
    # 阈值 2：C-z（缺失 + 额外各 1 条）也可用 -> 满配 3 个泳道
    res2 = resolve_batch_identify(req0.model_copy(update={"max_unmatched_bands": 2}))
    assert res2.assigned_lane_count == 3
    assert res2.unassigned_lanes == []
    assert res2.total_unmatched_count == 3


def test_batch_maximizes_cardinality_before_unmatched():
    # 贪心从“未匹配最少的边”出发会把 alpha 锁在 w1（u=0），再无谱型可配 w2；
    # 全局最优放弃 u=0 的边：beta→w1（u=1）、alpha→w2（u=1），配满 2 条
    req = _req(
        [("alpha", [200]), ("beta", [100, 200])],
        [("w1", [200]), ("w2", [200, 300])],
        tolerance=0,
        max_unmatched=1,
    )
    res = resolve_batch_identify(req)
    assert res.assigned_lane_count == 2
    assert _assignment_map(res) == {"w1": "beta", "w2": "alpha"}
    assert res.total_unmatched_count == 2


def test_batch_full_tie_on_all_metrics_falls_back_to_lexicographic():
    # 两个条带完全相同的谱型、两个观测完全相同的泳道：
    # tol=0 时 100↔100 配对，200 缺失、300 额外，对角线两种满配的未匹配总数
    # （各 2，共 4）与总误差（各 0）完全一致，最终只能由泳道序上的谱型名
    # 序列字典序决定：x→A、y→B
    req = _req(
        [("A", [100, 200]), ("B", [100, 200])],
        [("x", [100, 300]), ("y", [100, 300])],
        tolerance=0,
        max_unmatched=100,
    )
    res = resolve_batch_identify(req)
    assert res.assigned_lane_count == 2
    assert res.total_unmatched_count == 4
    assert res.total_error == 0
    assert _assignment_map(res) == {"x": "A", "y": "B"}


def test_batch_tiebreak_same_profile_bands_is_lexicographic():
    # 两个完全相同的谱型、两个完全相同的泳道：必须按谱型名字典序
    req = _req(
        [("zeta", [100, 200]), ("alpha", [100, 200])],
        [("l2", [100, 200]), ("l1", [100, 200])],
        tolerance=0,
        max_unmatched=0,
    )
    res = resolve_batch_identify(req)
    assert [item.lane for item in res.assignments] == ["l1", "l2"]
    assert _assignment_map(res) == {"l1": "alpha", "l2": "zeta"}


def test_batch_tiebreak_unassigned_sorts_after_all_profile_names():
    # 三个相同谱型（a、b、c）、两个相同泳道：满配 2 条，未分配项排在
    # 所有谱型名之后 -> 两个已用谱型必为字典序最小的 a、b
    req = _req(
        [("c", [100]), ("a", [100]), ("b", [100])],
        [("y", [100]), ("x", [100])],
        tolerance=0,
        max_unmatched=0,
    )
    res = resolve_batch_identify(req)
    assert _assignment_map(res) == {"x": "a", "y": "b"}
    assert res.unused_profiles == ["c"]
    assert res.unassigned_lanes == []


def test_batch_result_is_input_order_independent():
    req_a = _req(
        [("A", [100, 300]), ("B", [150, 400])],
        [("l1", [100, 300]), ("l2", [150, 400]), ("l3", [999])],
        tolerance=0,
        max_unmatched=0,
    )
    req_b = _req(
        [("B", [150, 400]), ("A", [100, 300])],
        [("l3", [999]), ("l2", [150, 400]), ("l1", [100, 300])],
        tolerance=0,
        max_unmatched=0,
    )
    assert resolve_batch_identify(req_a).model_dump() == \
        resolve_batch_identify(req_b).model_dump()


def test_batch_no_eligible_edges_leaves_everything_unassigned():
    req = _req(
        [("A", [100]), ("B", [200])],
        [("x", [900]), ("y", [800])],
        tolerance=0,
        max_unmatched=0,
    )
    res = resolve_batch_identify(req)
    assert res.eligible_edge_count == 0
    assert res.assignments == []
    assert res.unassigned_lanes == ["x", "y"]
    assert res.unused_profiles == ["A", "B"]
    assert res.assigned_lane_count == 0
    assert res.total_unmatched_count == 0
    assert res.total_error == 0


def test_batch_one_to_one_same_profile_never_used_twice():
    # 单谱型完美匹配两个泳道，也只能分配一次，另一个泳道未分配
    req = _req(
        [("A", [100, 200])],
        [("w1", [100, 200]), ("w2", [100, 200])],
        tolerance=0,
        max_unmatched=0,
    )
    res = resolve_batch_identify(req)
    assert res.assigned_lane_count == 1
    used = [item.profile for item in res.assignments]
    assert used == ["A"]
    assert res.unassigned_lanes == ["w2"]  # 泳道按名排序，w1 字典序占先
    assert res.unused_profiles == []


def test_batch_total_error_tiebreak_when_unmatched_equal():
    # 两种满配的未匹配总数均为 0，但总误差不同：
    # A=[100], B=[102], x=[101], y=[103], tol=5：
    # 边误差 A-x=1, A-y=3, B-x=1, B-y=1
    # 方案 A→x,B→y: 误差 2；方案 A→y,B→x: 误差 4 -> 取前者
    req = _req(
        [("A", [100]), ("B", [102])],
        [("x", [101]), ("y", [103])],
        tolerance=5,
        max_unmatched=0,
    )
    res = resolve_batch_identify(req)
    assert _assignment_map(res) == {"x": "A", "y": "B"}
    assert res.total_unmatched_count == 0
    assert res.total_error == 2


# ---------------------------------------------------------------- API 端到端

def _payload(**kw):
    payload = {
        "profiles": [
            {"name": "alpha", "predicted_bands": [230, 320, 450]},
            {"name": "beta", "predicted_bands": [318, 460, 700]},
            {"name": "gamma", "predicted_bands": [300, 600, 900]},
        ],
        "lanes": [
            {"name": "lane-1", "observed_bands": [231, 319, 455]},
            {"name": "lane-2", "observed_bands": [200, 322, 451]},
            {"name": "lane-3", "observed_bands": [317, 462, 698]},
        ],
        "band_match_tolerance": 10,
        "max_unmatched_bands": 2,
    }
    payload.update(kw)
    return payload


def test_api_batch_identify_full_response():
    res = client.post("/api/v1/digest/batch-identify", json=_payload())
    assert res.status_code == 200
    body = res.json()

    assert body["profile_count"] == 3
    assert body["lane_count"] == 3
    assert body["evaluated_edge_count"] == 9
    assert body["band_match_tolerance"] == 10
    assert body["max_unmatched_bands"] == 2

    assert [a["lane"] for a in body["assignments"]] == ["lane-1", "lane-3"]
    assert [a["profile"] for a in body["assignments"]] == ["alpha", "beta"]
    assert body["unassigned_lanes"] == ["lane-2"]
    assert body["unused_profiles"] == ["gamma"]
    assert body["assigned_lane_count"] == 2
    assert body["total_unmatched_count"] == 0
    assert body["assignment_reason"]

    first = body["assignments"][0]
    match = first["match"]
    assert match["matched_count"] == 3
    assert match["missing_bands"] == [] and match["extra_bands"] == []
    assert match["total_error"] == 7  # |230-231| + |320-319| + |450-455|
    assert [m["predicted_length"] for m in match["matches"]] == [230.0, 320.0, 450.0]
    for m in match["matches"]:
        assert m["error"] == abs(m["predicted_length"] - m["observed_length"])
        assert m["description"]


def test_api_batch_identify_assignments_sorted_by_lane_name():
    # 调换输入提交顺序，响应仍按泳道名排序
    payload = _payload(lanes=list(reversed(_payload()["lanes"])))
    body = client.post("/api/v1/digest/batch-identify", json=payload).json()
    assert [a["lane"] for a in body["assignments"]] == ["lane-1", "lane-3"]
    assert body["unassigned_lanes"] == ["lane-2"]
    assert body["unused_profiles"] == ["gamma"]


# ---------------------------------------------------------------- 422 校验

def _expect_422(payload, code_part=None):
    res = client.post("/api/v1/digest/batch-identify", json=payload)
    assert res.status_code == 422
    body = res.json()
    assert body["code"] == "validation_failed"
    assert body["errors"]
    if code_part:
        assert any(code_part in e["type"] for e in body["errors"])
    return body


def test_batch_422_duplicate_profile_and_lane_names():
    _expect_422(_payload(profiles=[
        {"name": "A", "predicted_bands": [100]},
        {"name": "A", "predicted_bands": [200]},
    ]), "duplicate_profile_name")
    _expect_422(_payload(lanes=[
        {"name": "x", "observed_bands": [100]},
        {"name": "x", "observed_bands": [200]},
    ]), "duplicate_lane_name")


def test_batch_422_blank_names():
    _expect_422(_payload(profiles=[
        {"name": "  ", "predicted_bands": [100]},
        {"name": "B", "predicted_bands": [200]},
    ]), "invalid_profile_name")


def test_batch_422_empty_lists_and_missing_fields():
    _expect_422(_payload(profiles=[]))
    _expect_422(_payload(lanes=[]))
    _expect_422(_payload(profiles=[{"name": "A", "predicted_bands": []}]))
    _expect_422(_payload(lanes=[{"name": "x", "observed_bands": []}]))
    payload = _payload()
    del payload["max_unmatched_bands"]
    _expect_422(payload)
    del payload["band_match_tolerance"]
    _expect_422(payload)


def test_batch_422_negative_values_and_bad_types():
    _expect_422(_payload(band_match_tolerance=-1))
    _expect_422(_payload(max_unmatched_bands=-1))
    _expect_422(_payload(profiles=[
        {"name": "A", "predicted_bands": [100, -5]},
    ]))
    _expect_422(_payload(max_unmatched_bands="two"))


def test_batch_422_too_many_profiles_or_lanes(monkeypatch):
    from app.config import settings
    monkeypatch.setattr(settings, "MAX_BATCH_PROFILES", 3)
    monkeypatch.setattr(settings, "MAX_BATCH_LANES", 3)
    _expect_422(_payload(profiles=[
        {"name": f"P{i}", "predicted_bands": [100 + i]} for i in range(4)
    ]), "too_many_batch_profiles")
    _expect_422(_payload(lanes=[
        {"name": f"L{i}", "observed_bands": [100 + i]} for i in range(4)
    ]), "too_many_batch_lanes")
