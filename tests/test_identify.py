"""单元测试：条带一对一最小代价匹配与样本判定纯算法。"""
from app.identify import decide, identify, match_sample
from app.models import IdentifyRequest


def test_match_basic_all_matched_within_tolerance():
    result = match_sample([450, 320, 230], [455, 318, 231], tolerance=10)
    assert result.matched_count == 3
    assert result.missing_count == 0 and result.extra_count == 0
    assert result.total_error == 8  # 5 + 2 + 1
    # 逐条匹配依据按预测长度升序
    assert [m.predicted_length for m in result.matches] == [230, 320, 450]
    assert [m.error for m in result.matches] == [1, 2, 5]
    assert all(m.description for m in result.matches)


def test_match_respects_tolerance_cutoff():
    result = match_sample([100], [130], tolerance=20)
    assert result.matched_count == 0
    assert result.missing_bands == [100]
    assert result.extra_bands == [130]
    assert result.total_error == 0
    # 容差放宽到恰好 30 即可配对（边界含等号）
    relaxed = match_sample([100], [130], tolerance=30)
    assert relaxed.matched_count == 1
    assert relaxed.total_error == 30


def test_match_is_one_to_one_and_min_cost():
    # 两条预测争一条观测：只配一条，且取代价最小者
    result = match_sample([100, 104], [101], tolerance=5)
    assert result.matched_count == 1
    assert result.matches[0].predicted_length == 100  # 误差 1 < 3
    assert result.missing_bands == [104]
    # 代价并列时确定性地取较短预测条带
    tie = match_sample([100, 102], [101], tolerance=5)
    assert tie.matches[0].predicted_length == 100
    assert tie.missing_bands == [102]
    # 一条预测争两条观测：同样一对一，取较短观测条带
    mirror = match_sample([101], [100, 102], tolerance=5)
    assert mirror.matched_count == 1
    assert mirror.matches[0].observed_length == 100
    assert mirror.extra_bands == [102]


def test_match_prefers_more_pairs_over_lower_cost():
    # 若只追求单对最小代价会只配 108↔109（误差 1），最优应为两对全配
    result = match_sample([100, 108], [104, 109], tolerance=5)
    assert result.matched_count == 2
    assert result.total_error == 5  # |100-104| + |108-109|


def test_match_input_order_does_not_matter():
    a = match_sample([450, 230, 320], [231, 455, 318], tolerance=10)
    b = match_sample([230, 320, 450], [318, 231, 455], tolerance=10)
    assert a.model_dump() == b.model_dump()


def test_match_float_lengths_and_rounding():
    result = match_sample([100, 200], [100.4, 199.6], tolerance=1)
    assert result.matched_count == 2
    assert [m.error for m in result.matches] == [0.4, 0.4]
    assert result.total_error == 0.8


def test_match_empty_prediction_everything_is_extra():
    result = match_sample([], [100, 200], tolerance=10)
    assert result.matched_count == 0
    assert result.missing_bands == []
    assert result.extra_bands == [100, 200]


def test_decide_prefers_fewer_unmatched():
    ref = match_sample([230, 320, 450], [231, 318, 455], tolerance=10)  # 全配
    var = match_sample([200, 500], [231, 318, 455], tolerance=10)       # 全配不上
    verdict, reason = decide(ref, var)
    assert verdict == "reference"
    assert "未匹配条带" in reason


def test_decide_tie_breaks_on_total_error():
    ref = match_sample([100, 200], [101, 202], tolerance=10)   # 总误差 1+2=3
    var = match_sample([105, 208], [101, 202], tolerance=10)   # 总误差 4+6=10
    verdict, reason = decide(ref, var)
    assert verdict == "reference"
    assert "总误差" in reason
    # 交换后判定随之翻转
    verdict_swapped, _ = decide(var, ref)
    assert verdict_swapped == "variant"


def test_decide_ambiguous_when_identical():
    ref = match_sample([100, 200], [101, 202], tolerance=10)
    var = match_sample([100, 200], [101, 202], tolerance=10)
    verdict, reason = decide(ref, var)
    assert verdict == "ambiguous"
    assert "无法区分" in reason


def test_identify_end_to_end():
    req = IdentifyRequest(
        reference_bands=[230, 320, 450],
        variant_bands=[200, 300, 500],
        observed_bands=[231, 318, 455],
        band_match_tolerance=10,
    )
    res = identify(req)
    assert res.verdict == "reference"
    assert res.observed_band_count == 3
    assert res.samples["reference"].matched_count == 3
    assert res.samples["variant"].matched_count == 0
    assert res.samples["variant"].missing_bands == [200, 300, 500]
    assert res.samples["variant"].extra_bands == [231, 318, 455]
