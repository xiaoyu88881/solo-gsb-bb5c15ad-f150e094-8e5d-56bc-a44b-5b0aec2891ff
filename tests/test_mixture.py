"""混合样本解析：共迁移聚带、子集枚举与最优子集选定的单元 + API 端到端测试。"""
from fastapi.testclient import TestClient

from app.main import app
from app.mixture import cluster_predicted_bands, resolve_mixture
from app.models import MixtureRequest

client = TestClient(app)


# ---------------------------------------------------------------- 单连接聚带（纯函数）

def test_cluster_single_linkage_chains_bands():
    # 100 -> 108 -> 116 逐级差 8（链尾与链头差 16 > 10），仍为同一簇
    assert cluster_predicted_bands([100, 116, 108], tolerance=10) == [108.0]


def test_cluster_splits_when_gap_exceeds_tolerance():
    assert cluster_predicted_bands([100, 108, 130], tolerance=10) == [104.0, 130.0]
    # 边界含等号：恰好相差容差时不分裂
    assert cluster_predicted_bands([100, 110], tolerance=10) == [105.0]


def test_cluster_comigrating_bands_appear_once_with_mean_length():
    # 两个候选共迁移到同一长度：只出现一次，代表长度取均值
    assert cluster_predicted_bands([100, 300, 104, 298], tolerance=10) == [102.0, 299.0]


def test_cluster_representative_mean_is_rounded_to_two_decimals():
    assert cluster_predicted_bands([100, 101, 102], tolerance=10) == [101.0]
    assert cluster_predicted_bands([100, 101, 101], tolerance=10) == [100.67]


def test_cluster_is_order_independent():
    assert cluster_predicted_bands([100, 108, 200], tolerance=10) == \
        cluster_predicted_bands([200, 100, 108], tolerance=10)


# ---------------------------------------------------------------- 子集选定（纯函数）

def _req(components, observed, tolerance=10, max_components=None):
    return MixtureRequest(
        components=[{"name": n, "predicted_bands": b} for n, b in components],
        observed_bands=observed,
        band_match_tolerance=tolerance,
        max_components=len(components) if max_components is None else max_components,
    )


def test_resolve_mixture_two_components_both_required():
    # 单候选都有未匹配条带，两者合并（318/320 共迁移）后全部匹配
    req = _req(
        [("A", [230, 320, 450]), ("B", [318, 460, 700])],
        observed=[231, 319, 455, 702],
    )
    res = resolve_mixture(req)
    assert res.selected_components == ["A", "B"]
    assert res.evaluated_subset_count == 3  # {A}、{B}、{A,B}
    # 共迁移条带 318/320 合并为 319；450/460 差 10 也合并为 455
    assert res.merged_predicted_bands == [230.0, 319.0, 455.0, 700.0]
    assert res.match.matched_count == 4
    assert res.match.missing_count == 0 and res.match.extra_count == 0
    assert res.match.total_error == 3  # 1 + 0 + 0 + 2


def test_resolve_mixture_prefers_fewer_unmatched_then_lower_error():
    # A 单独：未匹配 1 条（200 观测为额外）；B 单独与 A+B：未匹配 0 条，
    # A 与 B 预测相同但 B 整体更贴近观测 -> 总误差决胜选 B
    req = _req(
        [("A", [100, 300, 500]), ("B", [101, 300, 500])],
        observed=[105, 300, 500],
        tolerance=10,
    )
    res = resolve_mixture(req)
    assert res.selected_components == ["B"]
    match = res.match
    assert match.matched_count == 3
    assert match.total_error == 4


def test_resolve_mixture_prefers_fewer_components_on_full_tie():
    # 候选 C 的条带是 A、B 的超集；任何含 C 的子集与 C 单独结果完全一致
    req = _req(
        [("A", [100]), ("B", [300]), ("C", [100, 300])],
        observed=[100, 300],
        tolerance=0,
    )
    res = resolve_mixture(req)
    assert res.selected_components == ["C"]
    assert res.match.total_error == 0
    assert res.evaluated_subset_count == 7


def test_resolve_mixture_leverages_comigration_one_to_one():
    # A、B 在 ~200 处共迁移（只产生 1 条合并带），观测同样只有 1 条 -> 全部匹配
    req = _req(
        [("A", [200, 400]), ("B", [204, 600])],
        observed=[202, 400, 600],
        tolerance=10,
    )
    res = resolve_mixture(req)
    assert res.selected_components == ["A", "B"]
    assert res.merged_predicted_bands == [202.0, 400.0, 600.0]
    assert res.match.matched_count == 3
    assert res.match.missing_bands == [] and res.match.extra_bands == []


def test_resolve_mixture_max_components_limits_enumeration():
    # 上限 1 时只评估单候选子集：A、B 都有 1 条未匹配，A 总误差 7 < B 的 8 -> 选 A
    req = _req(
        [("A", [230, 320, 450]), ("B", [318, 460, 700])],
        observed=[231, 319, 455, 702],
        max_components=1,
    )
    res = resolve_mixture(req)
    assert res.selected_components == ["A"]
    assert res.evaluated_subset_count == 2
    # 上限放大到 2 后，混合子集（未匹配 0 条）胜出
    full = resolve_mixture(req.model_copy(update={"max_components": 2}))
    assert full.selected_components == ["A", "B"]
    assert full.evaluated_subset_count == 3


def test_resolve_mixture_final_tiebreak_is_lexicographic_name_list():
    # 两个完全相同的候选，单候选结果一致 -> 取字典序更小的名字
    req = _req(
        [("zeta", [100, 200]), ("alpha", [100, 200])],
        observed=[100, 200],
        tolerance=0,
    )
    res = resolve_mixture(req)
    assert res.selected_components == ["alpha"]
    assert res.selection_reason


# ---------------------------------------------------------------- API 端到端

def _mixture_payload(**kw):
    payload = {
        "components": [
            {"name": "strain-A", "predicted_bands": [230, 320, 450]},
            {"name": "strain-B", "predicted_bands": [318, 460, 700]},
        ],
        "observed_bands": [231, 319, 455, 702],
        "band_match_tolerance": 10,
        "max_components": 2,
    }
    payload.update(kw)
    return payload


def test_api_mixture_selects_mix_and_returns_full_evidence():
    res = client.post("/api/v1/digest/mixture", json=_mixture_payload())
    assert res.status_code == 200
    body = res.json()

    assert body["selected_components"] == ["strain-A", "strain-B"]
    assert body["merged_predicted_bands"] == [230.0, 319.0, 455.0, 700.0]
    assert body["evaluated_subset_count"] == 3
    assert body["observed_band_count"] == 4
    assert body["band_match_tolerance"] == 10
    assert body["max_components"] == 2
    assert body["selection_reason"]

    match = body["match"]
    assert match["matched_count"] == 4
    assert match["missing_count"] == 0 and match["extra_count"] == 0
    assert match["missing_bands"] == [] and match["extra_bands"] == []
    assert match["total_error"] == 3
    # 合并预测条带数与匹配统计自洽
    assert match["predicted_band_count"] == len(body["merged_predicted_bands"])
    assert match["matched_count"] == len(match["matches"])
    # 逐条匹配依据与合并条带、观测一一对应
    assert [m["predicted_length"] for m in match["matches"]] == [230.0, 319.0, 455.0, 700.0]
    assert [m["observed_length"] for m in match["matches"]] == [231, 319, 455, 702]
    for m in match["matches"]:
        assert m["error"] == abs(m["predicted_length"] - m["observed_length"])
        assert m["description"]


def test_api_mixture_counts_consistent_for_partial_match():
    body = client.post("/api/v1/digest/mixture", json=_mixture_payload(
        observed_bands=[231, 900],
    )).json()
    match = body["match"]
    assert match["matched_count"] + match["missing_count"] == match["predicted_band_count"]
    assert match["matched_count"] + match["extra_count"] == body["observed_band_count"]
    assert match["missing_count"] == len(match["missing_bands"])
    assert match["extra_count"] == len(match["extra_bands"])


# ---------------------------------------------------------------- 422 校验

def _expect_mixture_422(payload, code_part=None):
    res = client.post("/api/v1/digest/mixture", json=payload)
    assert res.status_code == 422
    body = res.json()
    assert body["code"] == "validation_failed"
    assert body["errors"]
    if code_part:
        assert any(code_part in e["type"] for e in body["errors"])
    return body


def test_mixture_422_max_components_zero_or_negative():
    _expect_mixture_422(_mixture_payload(max_components=0))
    _expect_mixture_422(_mixture_payload(max_components=-1))


def test_mixture_422_max_components_above_candidate_count():
    body = _expect_mixture_422(
        _mixture_payload(max_components=3), "max_components_too_large",
    )
    assert any("max_components" in e["loc"] for e in body["errors"])


def test_mixture_422_max_components_missing_or_wrong_type():
    payload = _mixture_payload()
    del payload["max_components"]
    _expect_mixture_422(payload)
    _expect_mixture_422(_mixture_payload(max_components="two"))


def test_mixture_422_duplicate_component_names():
    _expect_mixture_422(_mixture_payload(components=[
        {"name": "A", "predicted_bands": [100]},
        {"name": "A", "predicted_bands": [200]},
    ]), "duplicate_component_name")


def test_mixture_422_blank_component_name():
    _expect_mixture_422(_mixture_payload(components=[
        {"name": "   ", "predicted_bands": [100]},
        {"name": "B", "predicted_bands": [200]},
    ]), "invalid_component_name")


def test_mixture_422_empty_candidate_and_empty_request():
    # 候选自身预测条带为空
    _expect_mixture_422(_mixture_payload(components=[
        {"name": "A", "predicted_bands": [100]},
        {"name": "B", "predicted_bands": []},
    ]))
    # 候选列表为空
    _expect_mixture_422(_mixture_payload(components=[]))


def test_mixture_422_empty_observed():
    _expect_mixture_422(_mixture_payload(observed_bands=[]))


def test_mixture_422_negative_band_and_bad_tolerance():
    _expect_mixture_422(_mixture_payload(observed_bands=[100, -5]))
    _expect_mixture_422(_mixture_payload(components=[
        {"name": "A", "predicted_bands": [-1]},
        {"name": "B", "predicted_bands": [200]},
    ]))
    _expect_mixture_422(_mixture_payload(band_match_tolerance=-1))
    _expect_mixture_422(_mixture_payload(band_match_tolerance=None))


def test_mixture_422_too_many_candidates(monkeypatch):
    from app.config import settings
    monkeypatch.setattr(settings, "MAX_MIXTURE_CANDIDATES", 3)
    _expect_mixture_422(_mixture_payload(components=[
        {"name": f"C{i}", "predicted_bands": [100 + i]} for i in range(4)
    ], max_components=4), "too_many_mixture_candidates")


def test_mixture_422_reports_multiple_errors_at_once():
    # 重名 + 上限超候选数：一次性报告
    body = _expect_mixture_422(_mixture_payload(
        components=[
            {"name": "A", "predicted_bands": [100]},
            {"name": "A", "predicted_bands": [200]},
        ],
        max_components=5,
    ))
    types = {e["type"] for e in body["errors"]}
    assert "duplicate_component_name" in types
    assert "max_components_too_large" in types
