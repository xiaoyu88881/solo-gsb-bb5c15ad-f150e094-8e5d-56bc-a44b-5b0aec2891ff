"""批量盲样分配：逐边匹配筛边与全局一对一最优分配的单元 + API 端到端测试。"""
from fastapi.testclient import TestClient

from app.batch import assign_lanes, evaluate_edges, hungarian
from app.main import app
from app.models import BatchIdentifyRequest

client = TestClient(app)


# ---------------------------------------------------------------- 纯函数

def _req(profiles, lanes, tolerance=10, max_unmatched=2):
    return BatchIdentifyRequest(
        profiles=[{"name": n, "predicted_bands": b} for n, b in profiles],
        lanes=[{"name": n, "observed_bands": b} for n, b in lanes],
        band_match_tolerance=tolerance,
        max_unmatched_bands=max_unmatched,
    )


def test_assign_matches_example_scenario():
    req = _req(
        [
            ("strain-A", [230, 320, 450]),
            ("strain-B", [215, 400, 520]),
            ("strain-C", [228, 318, 452, 700]),
            ("strain-D", [600, 800, 950]),
        ],
        [
            ("lane-1", [231, 318, 455]),
            ("lane-2", [214, 402, 518]),
            ("lane-3", [229, 451, 701]),
        ],
    )
    res = assign_lanes(req)
    pairs = [(a.lane, a.assigned_profile) for a in res.assignments]
    assert pairs == [
        ("lane-1", "strain-A"),
        ("lane-2", "strain-B"),
        ("lane-3", "strain-C"),
    ]
    assert res.assigned_lane_count == 3
    assert res.unassigned_lanes == []
    assert res.unused_profiles == ["strain-D"]
    # lane-3 vs strain-C：229-228=1、452-451=1、701-700=1，318 缺失 -> 未匹配 1、误差 3
    a3 = res.assignments[2]
    assert a3.unmatched_count == 1
    assert a3.total_error == 3.0
    assert res.total_unmatched_count == 1
    assert res.total_error == 16.0  # lane-1:8 + lane-2:5 + lane-3:3
    assert res.assignment_reason


def test_threshold_marks_edges_unusable():
    req = _req(
        [("P", [100, 200, 300])],
        [("L", [100, 500, 600])],  # 未匹配 4 条
        tolerance=5,
        max_unmatched=2,
    )
    res = assign_lanes(req)
    edge = res.edges[0]
    assert edge.usable is False
    assert edge.unmatched_count == 4
    assert edge.match is None
    assert "超过阈值" in edge.reason
    assert res.assigned_lane_count == 0
    assert res.unassigned_lanes == ["L"]
    assert res.unused_profiles == ["P"]
    # 未分配泳道的观测条带全部计入未匹配
    assert res.assignments[0].unmatched_count == 3
    assert res.total_unmatched_count == 3
    assert res.total_error == 0
    # 恰好等于阈值的边可用（边界含等号）
    req_eq = _req([("P", [100])], [("L", [500, 600])], tolerance=0, max_unmatched=3)
    assert assign_lanes(req_eq).edges[0].usable is True


def test_maximizes_assigned_lanes_first():
    # 只有一个谱型：lane-1 完美匹配，lane-2 有 1 条未匹配；二者不能同配一个谱型
    req = _req(
        [("P", [100, 200, 300])],
        [("L1", [100, 200, 300]), ("L2", [100, 200, 900])],
        tolerance=0,
        max_unmatched=2,
    )
    res = assign_lanes(req)
    mapping = {a.lane: a.assigned_profile for a in res.assignments}
    assert mapping == {"L1": "P", "L2": None}
    assert res.assigned_lane_count == 1
    # 未匹配总数 = 已配边 0 + 未分配泳道全部 3 条观测
    assert res.total_unmatched_count == 3
    assert "优先用于其他泳道" in res.assignments[1].description


def test_minimizes_unmatched_before_error_and_lex():
    # 一个谱型两条泳道争：配 lane-A 未匹配 1（+未分配 lane-B 的 2 条 = 3），
    # 配 lane-B 未匹配 0（+未分配 lane-A 的 3 条 = 3），并列 -> 字典序 A 得谱型
    req = _req(
        [("P", [100, 200])],
        [("lane-A", [100, 200, 900]), ("lane-B", [100, 200])],
        tolerance=0,
        max_unmatched=2,
    )
    res = assign_lanes(req)
    mapping = {a.lane: a.assigned_profile for a in res.assignments}
    assert mapping == {"lane-A": "P", "lane-B": None}
    assert res.total_unmatched_count == 3

    # 真正拉开未匹配总数时，字典序让步
    req2 = _req(
        [("P", [100, 200, 300])],
        [("lane-A", [100, 200, 999]), ("lane-B", [100, 200, 300])],
        tolerance=0,
        max_unmatched=3,
    )
    res2 = assign_lanes(req2)
    mapping2 = {a.lane: a.assigned_profile for a in res2.assignments}
    assert mapping2 == {"lane-A": None, "lane-B": "P"}


def test_total_error_breaks_unmatched_tie():
    req = _req(
        [("P1", [105, 200]), ("P2", [101, 200])],
        [("L", [100, 200])],
        tolerance=10,
        max_unmatched=2,
    )
    res = assign_lanes(req)
    assert res.assignments[0].assigned_profile == "P2"
    assert res.total_error == 1


def test_lexicographic_tie_break_unassigned_ranks_after_all_profiles():
    # 3 个完全等价谱型、3 个完全等价泳道：按泳道序依次取最小谱型名
    req = _req(
        [(n, [100, 200]) for n in ("zeta", "alpha", "mid")],
        [(f"L{i}", [100, 200]) for i in (3, 1, 2)],
        tolerance=0,
        max_unmatched=0,
    )
    res = assign_lanes(req)
    pairs = [(a.lane, a.assigned_profile) for a in res.assignments]
    assert pairs == [("L1", "alpha"), ("L2", "mid"), ("L3", "zeta")]

    # 谱型多于泳道：两个泳道取前两个小谱型名，余下为未使用谱型
    req2 = _req(
        [(n, [100]) for n in ("c", "a", "b")],
        [(n, [100]) for n in ("Z", "A")],
        tolerance=0,
        max_unmatched=0,
    )
    res2 = assign_lanes(req2)
    pairs2 = [(a.lane, a.assigned_profile) for a in res2.assignments]
    assert pairs2 == [("A", "a"), ("Z", "b")]
    assert res2.unused_profiles == ["c"]


def test_result_is_order_independent():
    base_profiles = [("P1", [230, 320, 450]), ("P2", [231, 319, 455]), ("P3", [100, 800])]
    base_lanes = [("L1", [229, 321, 448]), ("L2", [102, 799]), ("L3", [232, 318])]
    ref = assign_lanes(_req(base_profiles, base_lanes)).model_dump()
    shuffled = assign_lanes(_req(
        list(reversed(base_profiles)), list(reversed(base_lanes)),
    )).model_dump()
    assert shuffled == ref


def test_assignments_sorted_by_lane_and_edges_sorted():
    req = _req(
        [("b", [100]), ("a", [100])],
        [("Z", [100]), ("A", [100]), ("M", [100])],
        tolerance=0,
        max_unmatched=0,
    )
    res = assign_lanes(req)
    # 仅 2 个谱型：字典序最优谱型名序列为 [a, b, 未分配]（未分配排在所有谱型名后）
    pairs = [(a.lane, a.assigned_profile) for a in res.assignments]
    assert pairs == [("A", "a"), ("M", "b"), ("Z", None)]
    edge_keys = [(e.lane, e.profile) for e in res.edges]
    assert edge_keys == sorted(edge_keys)
    assert res.unassigned_lanes == ["Z"]
    assert res.unused_profiles == []


def test_counts_and_evidence_are_consistent():
    req = _req(
        [("A", [230, 320, 450]), ("B", [100, 900])],
        [("x", [231, 319, 455]), ("y", [102, 898, 500])],
        tolerance=10,
        max_unmatched=2,
    )
    res = assign_lanes(req)
    # 全局汇总 = 逐项之和
    assert sum(a.unmatched_count for a in res.assignments) == res.total_unmatched_count
    assert round(sum(a.total_error for a in res.assignments), 2) == res.total_error
    assigned = [a for a in res.assignments if a.assigned_profile]
    assert len(assigned) == res.assigned_lane_count
    # 已分配带逐条依据，未分配不带；计数自洽
    for a in res.assignments:
        if a.assigned_profile:
            assert a.match is not None
            assert a.match.matched_count + a.match.missing_count == a.match.predicted_band_count
        else:
            assert a.match is None
    # 边依据：可用边带匹配明细，不可用边不带
    for e in res.edges:
        assert (e.match is not None) is e.usable
        assert e.missing_count + e.extra_count == e.unmatched_count


def test_profile_and_lane_names_share_namespace_names_allowed():
    req = _req([("X", [100])], [("X", [100])], tolerance=0, max_unmatched=0)
    res = assign_lanes(req)
    assert res.assignments[0].assigned_profile == "X"


def test_total_error_is_exact_sum_of_rounded_edge_errors():
    # 两条边各 0.1 bp 误差：全局总误差须精确等于 0.2（按百分之一 bp 整数编码）
    req = BatchIdentifyRequest(
        profiles=[{"name": "P1", "predicted_bands": [100.0]},
                  {"name": "P2", "predicted_bands": [200.0]}],
        lanes=[{"name": "L1", "observed_bands": [100.1]},
               {"name": "L2", "observed_bands": [200.1]}],
        band_match_tolerance=1,
        max_unmatched_bands=0,
    )
    res = assign_lanes(req)
    assert res.total_error == 0.2
    assert round(sum(a.total_error for a in res.assignments), 2) == res.total_error


def test_one_profile_many_lanes_full_tie_goes_to_lex_smallest_lane():
    req = BatchIdentifyRequest(
        profiles=[{"name": "P", "predicted_bands": [100]}],
        lanes=[{"name": f"L{i:02d}", "observed_bands": [100]} for i in range(32)],
        band_match_tolerance=0,
        max_unmatched_bands=0,
    )
    res = assign_lanes(req)
    assert res.assigned_lane_count == 1
    assert res.assignments[0].lane == "L00"
    assert res.assignments[0].assigned_profile == "P"
    assert all(a.assigned_profile is None for a in res.assignments[1:])
    assert res.unassigned_lanes == [f"L{i:02d}" for i in range(1, 32)]
    assert res.total_unmatched_count == 31


# ---------------------------------------------------------------- 匈牙利

def test_hungarian_identity_and_permutation():
    assert hungarian([[5]]) == (0,)
    # 唯一最优在置换位置
    cost = [
        [9, 1, 9],
        [1, 9, 9],
        [9, 9, 1],
    ]
    assert hungarian(cost) == (1, 0, 2)


def test_hungarian_rectangular_padded_case():
    # 1 真实行 + 1 虚拟行、1 真实列 + 1 虚拟列：真实行取更优的真实列
    # cost(row,col): 真实-真实=3, 真实-虚拟=5, 其余 0
    cost = [[3, 5], [0, 0]]
    assert hungarian(cost) == (0, 1)
    # 真实列比虚拟列更差时真实行宁可"未分配"
    cost = [[8, 5], [0, 0]]
    assert hungarian(cost) == (1, 0)


# ---------------------------------------------------------------- API 端到端

def _payload(**kw):
    payload = {
        "profiles": [
            {"name": "strain-A", "predicted_bands": [230, 320, 450]},
            {"name": "strain-B", "predicted_bands": [215, 400, 520]},
            {"name": "strain-C", "predicted_bands": [228, 318, 452, 700]},
            {"name": "strain-D", "predicted_bands": [600, 800, 950]},
        ],
        "lanes": [
            {"name": "lane-1", "observed_bands": [231, 318, 455]},
            {"name": "lane-2", "observed_bands": [214, 402, 518]},
            {"name": "lane-3", "observed_bands": [229, 451, 701]},
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

    assert body["band_match_tolerance"] == 10
    assert body["max_unmatched_bands"] == 2
    assert body["profile_count"] == 4 and body["lane_count"] == 3
    assert body["assigned_lane_count"] == 3
    assert body["total_unmatched_count"] == 1
    assert body["total_error"] == 16.0
    assert body["unassigned_lanes"] == []
    assert body["unused_profiles"] == ["strain-D"]

    # 逐项匹配依据
    a1 = body["assignments"][0]
    assert a1["lane"] == "lane-1" and a1["assigned_profile"] == "strain-A"
    assert a1["usable"] is True
    assert a1["unmatched_count"] == 0 and a1["total_error"] == 8
    assert [m["predicted_length"] for m in a1["match"]["matches"]] == [230, 320, 450]
    assert a1["description"]

    # 全部边均给出依据（4 谱型 × 3 泳道 = 12），按泳道、谱型排序
    assert len(body["edges"]) == 12
    blocked = [e for e in body["edges"] if not e["usable"]]
    assert blocked
    for e in body["edges"]:
        assert e["reason"]
        if e["usable"]:
            assert e["match"] is not None
            assert e["unmatched_count"] <= 2
        else:
            assert e["match"] is None


def test_api_batch_identify_order_independent():
    body1 = client.post("/api/v1/digest/batch-identify", json=_payload()).json()
    shuffled = _payload(
        profiles=list(reversed(_payload()["profiles"])),
        lanes=list(reversed(_payload()["lanes"])),
    )
    body2 = client.post("/api/v1/digest/batch-identify", json=shuffled).json()
    assert body1 == body2


def test_api_batch_identify_unassigned_lane_evidence():
    payload = _payload(
        profiles=[{"name": "P", "predicted_bands": [100, 200]}],
        lanes=[
            {"name": "L1", "observed_bands": [100, 200]},
            {"name": "L2", "observed_bands": [500, 600, 700]},
        ],
        max_unmatched_bands=1,
    )
    body = client.post("/api/v1/digest/batch-identify", json=payload).json()
    assert body["assigned_lane_count"] == 1
    assert body["unassigned_lanes"] == ["L2"]
    assert body["unused_profiles"] == []
    a2 = body["assignments"][1]
    assert a2["assigned_profile"] is None
    assert a2["usable"] is False
    assert a2["match"] is None
    assert a2["unmatched_count"] == 3
    assert "无可用边" in a2["description"]
    assert body["total_unmatched_count"] == 3


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


def test_422_missing_fields_and_empty_lists():
    _expect_422({})
    _expect_422(_payload(profiles=[]))
    _expect_422(_payload(lanes=[]))


def test_422_negative_and_invalid_numeric_fields():
    _expect_422(_payload(band_match_tolerance=-1))
    _expect_422(_payload(band_match_tolerance=None))
    _expect_422(_payload(max_unmatched_bands=-1))
    bad = _payload(profiles=[{"name": "P", "predicted_bands": [100, -5]}])
    _expect_422(bad)
    bad_lane = _payload(lanes=[{"name": "L", "observed_bands": [0, -1]}])
    _expect_422(bad_lane)


def test_422_empty_band_lists_and_blank_names():
    _expect_422(_payload(profiles=[{"name": "P", "predicted_bands": []}]))
    _expect_422(_payload(lanes=[{"name": "L", "observed_bands": []}]))
    _expect_422(_payload(profiles=[
        {"name": "  ", "predicted_bands": [100]},
        {"name": "B", "predicted_bands": [200]},
    ]))


def test_422_duplicate_profile_or_lane_names():
    body = _expect_422(_payload(profiles=[
        {"name": "X", "predicted_bands": [100]},
        {"name": "X", "predicted_bands": [200]},
    ]), "duplicate_name")
    assert any("profiles" in e["loc"] for e in body["errors"])
    _expect_422(_payload(lanes=[
        {"name": "Y", "observed_bands": [100]},
        {"name": "Y", "observed_bands": [200]},
    ]), "duplicate_name")


def test_422_too_many_items(monkeypatch):
    from app.config import settings
    monkeypatch.setattr(settings, "MAX_BATCH_ITEMS", 3)
    _expect_422(_payload(profiles=[
        {"name": f"P{i}", "predicted_bands": [100]} for i in range(4)
    ]), "too_many_profiles")
    _expect_422(_payload(lanes=[
        {"name": f"L{i}", "observed_bands": [100]} for i in range(4)
    ]), "too_many_lanes")


def test_422_reports_multiple_errors_at_once():
    # 谱型重名 + 泳道重名：语义校验一次性报告（字段均合法以越过 Pydantic）
    body = _expect_422(_payload(
        profiles=[
            {"name": "X", "predicted_bands": [100]},
            {"name": "X", "predicted_bands": [200]},
        ],
        lanes=[
            {"name": "Y", "observed_bands": [100]},
            {"name": "Y", "observed_bands": [200]},
        ],
    ))
    types = {e["type"] for e in body["errors"]}
    assert types == {"duplicate_name"}
    locs = [e["loc"][0] for e in body["errors"]]
    assert "profiles" in locs and "lanes" in locs


def test_edge_evidence_reuses_identify_matching():
    # 容差边界：误差恰好 10 的边可匹配；未匹配数判定阈值
    edges = evaluate_edges(
        {"P": [100, 200]}, {"L": [110, 200, 300]},
        lanes=["L"], profiles=["P"], tolerance=10, max_unmatched=1,
    )
    match, usable = edges[("L", "P")]
    assert match.matched_count == 2  # 100-110、200-200
    assert match.extra_count == 1    # 300 额外
    assert match.total_error == 10
    assert usable is True           # 未匹配恰好 1，阈值 1（边界含等号）
