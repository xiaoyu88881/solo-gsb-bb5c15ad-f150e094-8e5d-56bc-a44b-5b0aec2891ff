"""端到端 API 测试（FastAPI TestClient，不启动真实网络服务）。"""
from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)


def _enzyme(name, ref, var):
    return {"name": name, "reference_sites": ref, "variant_sites": var}


def _payload(enzymes, length=1000, min_visible=50, tolerance=20):
    return {
        "plasmid_length": length,
        "min_visible_fragment": min_visible,
        "band_resolution_tolerance": tolerance,
        "enzymes": enzymes,
    }


# ---------------------------------------------------------------- 健康检查

def test_health():
    res = client.get("/health")
    assert res.status_code == 200
    body = res.json()
    assert body["status"] == "ok"
    assert "version" in body


# ---------------------------------------------------------------- 基础推荐

def test_recommend_basic_scenario():
    payload = _payload([
        _enzyme("EcoRI", [100, 600], [100, 600, 850]),
        _enzyme("BamHI", [300], [300]),
    ])
    res = client.post("/api/v1/digest/recommend", json=payload)
    assert res.status_code == 200
    body = res.json()

    assert body["total_enzyme_pairs"] == 1
    assert body["recommended_count"] == 1
    rec = body["recommendations"][0]
    assert rec["enzyme_pair"] == ["BamHI", "EcoRI"]
    assert rec["distinguishable_band_count"] >= 1

    # 每种样本都应报告一个跨原点片段
    for sample in ("reference", "variant"):
        crosses = [f for f in rec["samples"][sample]["fragments"] if f["crosses_origin"]]
        assert len(crosses) == 1

    # 所有差异条带必须可区分
    assert all(d["reference_count"] != d["variant_count"] for d in rec["differences"])
    assert {d["difference_type"] for d in rec["differences"]} <= {
        "reference_only", "variant_only", "count_mismatch",
    }


def test_fragment_lengths_sum_to_plasmid_length():
    payload = _payload([
        _enzyme("EcoRI", [100, 600], [100, 600, 850]),
        _enzyme("BamHI", [300], [300]),
    ])
    rec = client.post("/api/v1/digest/recommend", json=payload).json()["recommendations"][0]
    for sample in ("reference", "variant"):
        fragments = rec["samples"][sample]["fragments"]
        assert sum(f["length"] for f in fragments) == 1000


# ---------------------------------------------------------------- 排序

def test_recommendations_sorting_distinguishable_desc_then_total_asc():
    # L=100、最小可见=1、容差=0，手工精确推演（片段多集合）：
    #   A+B：参考 {20,80} vs 变异 {20,30,50} -> 3 条可区分 / 共 4 条带
    #   A+C：参考 {20,80} vs 变异 {40,60}     -> 4 条可区分 / 共 4 条带
    #   C+D：参考 {30,70} vs 变异 {10,90}     -> 4 条可区分 / 共 4 条带（键序在后）
    #   A+D、B+D：两样本一致 -> 不可区分；B+C：参考仅 1 个合并切点 -> 拒绝
    payload = _payload(
        [
            _enzyme("A", [0], [0]),
            _enzyme("B", [20], [20, 50]),
            _enzyme("C", [20], [40]),
            _enzyme("D", [50], [50]),
        ],
        length=100,
        min_visible=1,
        tolerance=0,
    )
    body = client.post("/api/v1/digest/recommend", json=payload).json()

    assert body["total_enzyme_pairs"] == 6
    keys = [r["pair_key"] for r in body["recommendations"]]
    assert keys == ["A+C", "C+D", "A+B"]
    assert body["recommendations"][0]["distinguishable_band_count"] == 4
    assert body["recommendations"][1]["distinguishable_band_count"] == 4
    ab = body["recommendations"][2]
    assert (ab["distinguishable_band_count"], ab["total_bands"]) == (3, 4)

    assert [p["pair_key"] for p in body["non_discriminating_pairs"]] == ["A+D", "B+D"]
    assert body["non_discriminating_pair_count"] == 2
    assert body["rejected_pairs"][0]["pair_key"] == "B+C"
    assert body["rejected_pairs"][0]["reason_code"] == "insufficient_cuts"


def test_tie_break_on_total_bands_and_name():
    # 两对都不可区分以外的同分时按总带数升序、酶名升序
    payload = _payload([
        _enzyme("EcoRI", [0], [0, 200]),
        _enzyme("BamHI", [400], [400, 700]),
        _enzyme("HindIII", [100, 900], [100, 500, 900]),
    ])
    body = client.post("/api/v1/digest/recommend", json=payload).json()
    recs = body["recommendations"]
    for a, b in zip(recs, recs[1:]):
        ka = (-a["distinguishable_band_count"], a["total_bands"], a["pair_key"])
        kb = (-b["distinguishable_band_count"], b["total_bands"], b["pair_key"])
        assert ka <= kb


# ---------------------------------------------------------------- 拒绝逻辑

def test_rejected_when_sample_has_fewer_than_two_cuts():
    # HindIII 对参考型无切点 -> EcoRI+HindIII 对参考型仅 1 个切点
    payload = _payload([
        _enzyme("EcoRI", [100], [100, 600]),
        _enzyme("HindIII", [], [300]),
    ])
    body = client.post("/api/v1/digest/recommend", json=payload).json()
    assert body["recommended_count"] == 0
    assert body["rejected_pair_count"] == 1
    rejected = body["rejected_pairs"][0]
    assert rejected["reason_code"] == "insufficient_cuts"
    assert rejected["reference_cut_count"] == 1
    assert rejected["variant_cut_count"] == 3


def test_non_discriminating_pair_listed_separately():
    # 两酶在两种样本上切点完全相同 -> 无法区分
    payload = _payload([
        _enzyme("EcoRI", [100, 600], [100, 600]),
        _enzyme("BamHI", [300], [300]),
    ])
    body = client.post("/api/v1/digest/recommend", json=payload).json()
    assert body["recommended_count"] == 0
    assert body["non_discriminating_pair_count"] == 1
    assert body["non_discriminating_pairs"][0]["pair_key"] == "BamHI+EcoRI"


def test_coincident_sites_from_different_enzymes_merge():
    # 参考型两酶在 100 处同位点（合并，不按重复位点拒绝）；变异型 BamHI 在 350 处
    payload = _payload([
        _enzyme("EcoRI", [100, 600], [100, 600]),
        _enzyme("BamHI", [100, 300], [350]),
    ])
    body = client.post("/api/v1/digest/recommend", json=payload).json()
    assert body["recommended_count"] == 1
    rec = body["recommendations"][0]
    cuts = rec["samples"]["reference"]["cuts"]
    at_100 = [c for c in cuts if c["position"] == 100]
    assert len(at_100) == 1
    assert at_100[0]["enzymes"] == ["BamHI", "EcoRI"]
    var_cuts = rec["samples"]["variant"]["cuts"]
    assert var_cuts[0]["position"] == 100 and var_cuts[0]["enzymes"] == ["EcoRI"]


# ---------------------------------------------------------------- 422 校验

def _expect_422(payload, code_part=None):
    res = client.post("/api/v1/digest/recommend", json=payload)
    assert res.status_code == 422
    body = res.json()
    assert body["code"] == "validation_failed"
    assert body["errors"]
    if code_part:
        assert any(code_part in e["type"] for e in body["errors"])
    return body


def test_422_site_out_of_range():
    payload = _payload([
        _enzyme("EcoRI", [100, 1000], [100]),   # 1000 越界（合法 0..999）
        _enzyme("BamHI", [300], [300]),
    ])
    _expect_422(payload, "site_out_of_range")


def test_422_negative_site():
    payload = _payload([
        _enzyme("EcoRI", [-1], [100]),
        _enzyme("BamHI", [300], [300]),
    ])
    _expect_422(payload, "site_out_of_range")


def test_422_duplicate_site():
    payload = _payload([
        _enzyme("EcoRI", [100, 100, 600], [100]),
        _enzyme("BamHI", [300], [300]),
    ])
    _expect_422(payload, "duplicate_site")


def test_422_duplicate_enzyme_name():
    payload = _payload([
        _enzyme("EcoRI", [100], [100]),
        _enzyme("EcoRI", [300], [300]),
    ])
    _expect_422(payload, "duplicate_enzyme_name")


def test_422_min_visible_not_smaller_than_plasmid():
    payload = _payload(
        [_enzyme("EcoRI", [100], [100]), _enzyme("BamHI", [300], [300])],
        min_visible=1000,
    )
    _expect_422(payload, "min_visible_too_large")


def test_422_basic_field_constraints():
    # plasmid_length <= 0、容差为负、酶不足 2 种
    _expect_422(_payload(
        [_enzyme("EcoRI", [100], [100])],
        length=0,
    ))
    res = client.post("/api/v1/digest/recommend", json=_payload(
        [_enzyme("EcoRI", [100], [100]), _enzyme("BamHI", [300], [300])],
    ))
    # 正常对照
    assert res.status_code == 200


# ---------------------------------------------------------------- 鉴别（identify）

def _identify_payload(**kw):
    payload = {
        "reference_bands": [230, 320, 450],
        "variant_bands": [200, 300, 500],
        "observed_bands": [231, 318, 455],
        "band_match_tolerance": 10,
    }
    payload.update(kw)
    return payload


def test_identify_reference_verdict():
    res = client.post("/api/v1/digest/identify", json=_identify_payload())
    assert res.status_code == 200
    body = res.json()

    assert body["verdict"] == "reference"
    assert body["verdict_reason"]
    assert body["observed_band_count"] == 3
    assert body["band_match_tolerance"] == 10

    ref = body["samples"]["reference"]
    assert ref["predicted_band_count"] == 3
    assert ref["matched_count"] == 3
    assert ref["missing_count"] == 0 and ref["extra_count"] == 0
    assert ref["missing_bands"] == [] and ref["extra_bands"] == []
    assert ref["total_error"] == 8  # 1 + 2 + 5

    # 逐条匹配依据：按预测长度升序，误差与说明齐全
    assert [m["predicted_length"] for m in ref["matches"]] == [230, 320, 450]
    for m in ref["matches"]:
        assert m["error"] == abs(m["predicted_length"] - m["observed_length"])
        assert m["description"]

    var = body["samples"]["variant"]
    assert var["matched_count"] == 0
    assert var["missing_bands"] == [200, 300, 500]
    assert var["extra_bands"] == [231, 318, 455]
    assert var["total_error"] == 0


def test_identify_variant_verdict():
    body = client.post("/api/v1/digest/identify", json=_identify_payload(
        observed_bands=[205, 298, 502],
    )).json()
    assert body["verdict"] == "variant"
    assert body["samples"]["variant"]["matched_count"] == 3
    assert body["samples"]["reference"]["matched_count"] == 0


def test_identify_ambiguous_when_samples_equally_good():
    body = client.post("/api/v1/digest/identify", json=_identify_payload(
        variant_bands=[230, 320, 450],
    )).json()
    assert body["verdict"] == "ambiguous"
    assert "无法区分" in body["verdict_reason"]


def test_identify_one_to_one_and_order_independent():
    # 观测乱序、两条预测争一条观测：一对一匹配，且结果与提交顺序无关
    first = client.post("/api/v1/digest/identify", json=_identify_payload(
        reference_bands=[102, 100],
        variant_bands=[900],
        observed_bands=[101],
        band_match_tolerance=5,
    )).json()
    ref = first["samples"]["reference"]
    assert ref["matched_count"] == 1
    assert ref["matches"][0]["predicted_length"] == 100
    assert ref["missing_bands"] == [102]
    var = first["samples"]["variant"]
    assert var["missing_bands"] == [900]
    assert var["extra_bands"] == [101]
    assert first["verdict"] == "reference"

    second = client.post("/api/v1/digest/identify", json=_identify_payload(
        reference_bands=[100, 102],
        variant_bands=[900],
        observed_bands=[101],
        band_match_tolerance=5,
    )).json()
    assert first == second


def test_identify_counts_are_consistent():
    body = client.post("/api/v1/digest/identify", json=_identify_payload(
        reference_bands=[230, 320, 450, 700],
        observed_bands=[231, 318, 455, 900],
    )).json()
    for sample in ("reference", "variant"):
        s = body["samples"][sample]
        assert s["matched_count"] == len(s["matches"])
        assert s["matched_count"] + s["missing_count"] == s["predicted_band_count"]
        assert s["matched_count"] + s["extra_count"] == body["observed_band_count"]
        assert s["missing_count"] == len(s["missing_bands"])
        assert s["extra_count"] == len(s["extra_bands"])


def test_identify_float_bands():
    body = client.post("/api/v1/digest/identify", json=_identify_payload(
        reference_bands=[100, 200],
        variant_bands=[150],
        observed_bands=[100.4, 199.6],
        band_match_tolerance=1,
    )).json()
    ref = body["samples"]["reference"]
    assert ref["matched_count"] == 2
    assert ref["total_error"] == 0.8
    assert body["verdict"] == "reference"


def test_identify_empty_prediction_lists_allowed():
    # 参考/变异预测列表可为空（全部观测计为额外条带），两者一致则无法区分
    body = client.post("/api/v1/digest/identify", json=_identify_payload(
        reference_bands=[], variant_bands=[],
    )).json()
    assert body["verdict"] == "ambiguous"
    assert body["samples"]["reference"]["extra_bands"] == [231, 318, 455]


# ---------------------------------------------------------------- identify 422 校验

def _expect_identify_422(payload, code_part=None):
    res = client.post("/api/v1/digest/identify", json=payload)
    assert res.status_code == 422
    body = res.json()
    assert body["code"] == "validation_failed"
    assert body["errors"]
    if code_part:
        assert any(code_part in e["type"] for e in body["errors"])
    return body


def test_identify_422_empty_observed():
    body = _expect_identify_422(_identify_payload(observed_bands=[]))
    assert any("observed_bands" in e["loc"] for e in body["errors"])


def test_identify_422_negative_band():
    _expect_identify_422(_identify_payload(observed_bands=[100, -5]))
    _expect_identify_422(_identify_payload(reference_bands=[-1]))
    _expect_identify_422(_identify_payload(variant_bands=[-1]))


def test_identify_422_invalid_tolerance():
    _expect_identify_422(_identify_payload(band_match_tolerance=-1))
    _expect_identify_422(_identify_payload(band_match_tolerance="wide"))
    _expect_identify_422(_identify_payload(band_match_tolerance=None))


def test_identify_422_missing_required_fields():
    # 缺观测列表或缺容差
    _expect_identify_422({"reference_bands": [100], "band_match_tolerance": 10})
    _expect_identify_422({"observed_bands": [100]})
    _expect_identify_422({})


def test_identify_422_too_many_bands(monkeypatch):
    from app.config import settings
    monkeypatch.setattr(settings, "MAX_BANDS", 4)
    _expect_identify_422(
        _identify_payload(observed_bands=[1, 2, 3, 4, 5]), "too_many_bands",
    )
