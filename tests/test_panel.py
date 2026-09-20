"""最小酶对面板：精确集合覆盖单元测试 + API 端到端测试。"""
import random
from itertools import combinations

from fastapi.testclient import TestClient

from app.main import app
from app.panel import exact_min_panel

client = TestClient(app)


def _enzyme(name, ref, var_sites):
    return {"name": name, "reference_sites": ref, "variant_sites": var_sites}


def _payload(enzymes, length=1000, min_visible=50, tolerance=20):
    return {
        "plasmid_length": length,
        "min_visible_fragment": min_visible,
        "band_resolution_tolerance": tolerance,
        "enzymes": enzymes,
    }


# ---------------------------------------------------------------- 精确最小面板（纯函数）

def test_exact_min_panel_singletons():
    result = exact_min_panel([(0b01, 1, "A+B"), (0b10, 1, "C+D")], 2)
    assert result == (2, ("A+B", "C+D"))


def test_exact_min_panel_is_not_greedy():
    # 贪心（先选覆盖最多的 P+C）需要 3 对；精确解为 P+A、P+B 两对
    entries = [
        (0b00001111, 4, "P+A"),
        (0b11110000, 4, "P+B"),
        (0b01110111, 6, "P+C"),
    ]
    assert exact_min_panel(entries, 8) == (8, ("P+A", "P+B"))


def test_exact_min_panel_tie_break_prefers_more_bands():
    entries = [
        (0b01, 3, "B+C"),
        (0b10, 1, "B+D"),
        (0b01, 1, "B+E"),
        (0b10, 4, "B+F"),
    ]
    # 同规模（2 对）方案条带总数：B+C/B+F=7 最高，优于字典序更小的组合
    assert exact_min_panel(entries, 2) == (7, ("B+C", "B+F"))


def test_exact_min_panel_tie_break_lexicographic_keys():
    entries = [
        (0b01, 2, "B+C"),
        (0b10, 2, "F+G"),
        (0b01, 2, "A+C"),
        (0b10, 2, "B+D"),
    ]
    # 四个同规模同条带方案，取排序后 pair_key 列表字典序最小者
    assert exact_min_panel(entries, 2) == (4, ("A+C", "B+D"))


def test_exact_min_panel_dedupes_identical_masks_by_best_pair():
    entries = [
        (0b11, 1, "A+B"),
        (0b11, 5, "A+C"),   # 同掩码的更优代表（条带更多）
        (0b01, 9, "A+D"),
    ]
    assert exact_min_panel(entries, 2) == (5, ("A+C",))


def test_exact_min_panel_ignores_empty_masks():
    assert exact_min_panel([(0, 100, "A+B"), (0b11, 1, "C+D")], 2) == (1, ("C+D",))


def test_exact_min_panel_returns_none_when_uncoverable():
    assert exact_min_panel([(0b01, 1, "A+B")], 2) is None
    assert exact_min_panel([], 2) is None


def test_exact_min_panel_matches_brute_force():
    # 与暴力枚举全部子集交叉验证：规模最小 → 条带总数最多 → 排序 key 列表字典序最小
    def brute_force(entries, v):
        full = (1 << v) - 1
        entries = [(m, b, k) for m, b, k in entries if m]
        best = None
        for r in range(1, len(entries) + 1):
            for combo in combinations(entries, r):
                mask = 0
                for m, _, _ in combo:
                    mask |= m
                if mask != full:
                    continue
                bands = sum(b for _, b, _ in combo)
                keys = tuple(sorted(k for _, _, k in combo))
                cand = (r, -bands, keys)
                if best is None or cand < best:
                    best = cand
            if best is not None:
                return (-best[1], best[2])
        return None

    rng = random.Random(20260920)
    for _ in range(500):
        v = rng.randint(1, 5)
        entries = [
            (rng.randint(0, (1 << v) - 1), rng.randint(1, 5), f"K{i:02d}")
            for i in range(rng.randint(1, 8))
        ]
        assert exact_min_panel(entries, v) == brute_force(entries, v)


# ---------------------------------------------------------------- 基础面板

def test_panel_single_pair_covers_all_variants():
    payload = _payload([
        _enzyme("EcoRI", [100, 600], {"del100": [100, 600, 850], "ins250": [100, 600]}),
        _enzyme("BamHI", [300], {"del100": [300], "ins250": [300, 750]}),
        _enzyme("HindIII", [50], {"del100": [50], "ins250": [50]}),
    ])
    res = client.post("/api/v1/digest/panel", json=payload)
    assert res.status_code == 200
    body = res.json()

    assert body["plasmid_length"] == 1000
    assert body["variant_count"] == 2
    assert body["total_enzyme_pairs"] == 3
    assert body["uncovered_variants"] == []
    assert body["covered_variants"] == ["del100", "ins250"]

    # BamHI+EcoRI 一对即可覆盖两个变异型（del100 2 条、ins250 3 条可区分条带）
    assert [p["pair_key"] for p in body["selected_pairs"]] == ["BamHI+EcoRI"]
    assert body["total_distinguishable_bands"] == 5

    pair = body["selected_pairs"][0]
    assert pair["enzyme_pair"] == ["BamHI", "EcoRI"]
    assert pair["covered_variants"] == ["del100", "ins250"]
    assert pair["distinguishable_band_count"] == 5
    counts = {c["variant"]: c["distinguishable_band_count"] for c in pair["coverages"]}
    assert counts == {"del100": 2, "ins250": 3}
    # 逐变异的差异条带依据齐全且可区分
    for coverage in pair["coverages"]:
        assert len(coverage["differences"]) == coverage["distinguishable_band_count"]
        for diff in coverage["differences"]:
            assert diff["reference_count"] != diff["variant_count"]
            assert diff["description"]


def _marker_enzyme(name, special_variant, variants=("v1", "v2", "v3", "v4")):
    # 仅在 special_variant 上多一个切点，其余变异型与参考型一致
    return _enzyme(
        name,
        [100, 600],
        {v: ([100, 600, 850] if v == special_variant else [100, 600]) for v in variants},
    )


def test_panel_exact_minimum_and_lexicographic_tie_break():
    # 酶 Ei 只在变异型 vi 上与参考型不同；N 在所有样本上切点相同。
    # 最小面板为 2 对：{E1+E2,E3+E4}、{E1+E3,E2+E4}、{E1+E4,E2+E3} 三种同规模
    # 同条带方案，精确解取排序后 pair_key 列表字典序最小者。
    payload = _payload([
        _marker_enzyme("E1", "v1"),
        _marker_enzyme("E2", "v2"),
        _marker_enzyme("E3", "v3"),
        _marker_enzyme("E4", "v4"),
        _enzyme("N", [300], {v: [300] for v in ("v1", "v2", "v3", "v4")}),
    ])
    res = client.post("/api/v1/digest/panel", json=payload)
    assert res.status_code == 200
    body = res.json()

    assert body["variant_count"] == 4
    assert body["total_enzyme_pairs"] == 10
    assert body["uncovered_variants"] == []
    assert body["covered_variants"] == ["v1", "v2", "v3", "v4"]
    assert [p["pair_key"] for p in body["selected_pairs"]] == ["E1+E2", "E3+E4"]
    assert body["total_distinguishable_bands"] == 8
    for pair in body["selected_pairs"]:
        assert len(pair["covered_variants"]) == 2
        assert pair["distinguishable_band_count"] == 4


# ---------------------------------------------------------------- 未覆盖变异型

def test_panel_uncovered_variant_returns_200_without_partial_panel():
    # v2 与参考型切点完全一致 -> 任何酶对都无法区分 v2
    payload = _payload([
        _enzyme("EcoRI", [100, 600], {"v1": [100, 600, 850], "v2": [100, 600]}),
        _enzyme("BamHI", [300], {"v1": [300], "v2": [300]}),
    ])
    res = client.post("/api/v1/digest/panel", json=payload)
    assert res.status_code == 200
    body = res.json()
    assert body["uncovered_variants"] == ["v2"]
    # 不返回部分面板：即使 v1 可被覆盖，selected_pairs 也为空
    assert body["selected_pairs"] == []
    assert body["covered_variants"] == []
    assert body["total_distinguishable_bands"] == 0


def test_panel_uncovered_variant_due_to_insufficient_cuts():
    # v2 上两酶合并后仅 1 个切点，无法完成圆环双酶切 -> 该酶对对 v2 不可用
    payload = _payload([
        _enzyme("EcoRI", [100, 600], {"v1": [100, 600, 850], "v2": [100]}),
        _enzyme("BamHI", [300], {"v1": [300], "v2": []}),
    ])
    res = client.post("/api/v1/digest/panel", json=payload)
    assert res.status_code == 200
    body = res.json()
    assert body["uncovered_variants"] == ["v2"]
    assert body["selected_pairs"] == []


# ---------------------------------------------------------------- 422 校验

def _expect_panel_422(payload, code_part=None):
    res = client.post("/api/v1/digest/panel", json=payload)
    assert res.status_code == 422
    body = res.json()
    assert body["code"] == "validation_failed"
    assert body["errors"]
    if code_part:
        assert any(code_part in e["type"] for e in body["errors"])
    return body


def test_panel_422_variant_site_out_of_range():
    payload = _payload([
        _enzyme("EcoRI", [100], {"v1": [100, 1000]}),   # 1000 越界（合法 0..999）
        _enzyme("BamHI", [300], {"v1": [300]}),
    ])
    body = _expect_panel_422(payload, "site_out_of_range")
    assert any("variant_sites" in e["loc"] for e in body["errors"])


def test_panel_422_negative_variant_site():
    payload = _payload([
        _enzyme("EcoRI", [100], {"v1": [-1]}),
        _enzyme("BamHI", [300], {"v1": [300]}),
    ])
    _expect_panel_422(payload, "site_out_of_range")


def test_panel_422_duplicate_variant_site():
    payload = _payload([
        _enzyme("EcoRI", [100], {"v1": [100, 100, 600]}),
        _enzyme("BamHI", [300], {"v1": [300]}),
    ])
    _expect_panel_422(payload, "duplicate_site")


def test_panel_422_duplicate_reference_site():
    payload = _payload([
        _enzyme("EcoRI", [100, 100], {"v1": [100]}),
        _enzyme("BamHI", [300], {"v1": [300]}),
    ])
    _expect_panel_422(payload, "duplicate_site")


def test_panel_422_duplicate_enzyme_name():
    payload = _payload([
        _enzyme("EcoRI", [100], {"v1": [100]}),
        _enzyme("EcoRI", [300], {"v1": [300]}),
    ])
    _expect_panel_422(payload, "duplicate_enzyme_name")


def test_panel_422_inconsistent_variant_names():
    payload = _payload([
        _enzyme("EcoRI", [100], {"v1": [100]}),
        _enzyme("BamHI", [300], {"v2": [300]}),
    ])
    _expect_panel_422(payload, "inconsistent_variant_names")


def test_panel_422_blank_variant_name():
    payload = _payload([
        _enzyme("EcoRI", [100], {"": [100]}),
        _enzyme("BamHI", [300], {"": [300]}),
    ])
    _expect_panel_422(payload, "invalid_variant_name")


def test_panel_422_no_variants():
    payload = _payload([
        _enzyme("EcoRI", [100], {}),
        _enzyme("BamHI", [300], {}),
    ])
    _expect_panel_422(payload, "no_variants")


def test_panel_422_too_many_variants(monkeypatch):
    from app.config import settings
    monkeypatch.setattr(settings, "MAX_VARIANTS", 1)
    payload = _payload([
        _enzyme("EcoRI", [100], {"v1": [100], "v2": [100]}),
        _enzyme("BamHI", [300], {"v1": [300], "v2": [300]}),
    ])
    _expect_panel_422(payload, "too_many_variants")


def test_panel_422_too_many_enzymes(monkeypatch):
    from app.config import settings
    monkeypatch.setattr(settings, "MAX_ENZYMES", 1)
    payload = _payload([
        _enzyme("EcoRI", [100], {"v1": [100]}),
        _enzyme("BamHI", [300], {"v1": [300]}),
    ])
    _expect_panel_422(payload, "too_many_enzymes")


def test_panel_422_min_visible_not_smaller_than_plasmid():
    payload = _payload(
        [
            _enzyme("EcoRI", [100], {"v1": [100]}),
            _enzyme("BamHI", [300], {"v1": [300]}),
        ],
        min_visible=1000,
    )
    _expect_panel_422(payload, "min_visible_too_large")


def test_panel_422_basic_field_constraints():
    # 质粒长度 <= 0、容差为负、酶不足 2 种、variant_sites 类型错误
    _expect_panel_422(_payload(
        [_enzyme("EcoRI", [100], {"v1": [100]}), _enzyme("BamHI", [300], {"v1": [300]})],
        length=0,
    ))
    _expect_panel_422(_payload(
        [_enzyme("EcoRI", [100], {"v1": [100]}), _enzyme("BamHI", [300], {"v1": [300]})],
        tolerance=-1,
    ))
    _expect_panel_422(_payload([_enzyme("EcoRI", [100], {"v1": [100]})]))
    bad = _payload([
        {"name": "EcoRI", "reference_sites": [100], "variant_sites": [[100]]},
        _enzyme("BamHI", [300], {"v1": [300]}),
    ])
    _expect_panel_422(bad)
