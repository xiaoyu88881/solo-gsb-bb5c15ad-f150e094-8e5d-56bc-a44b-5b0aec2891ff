"""核心算法：凝胶观测条带与两种样本预测条带的一对一最小代价匹配及判定。

所有函数均为纯函数，不依赖任何外部服务或全局状态。

匹配规则
--------
对每种样本，将其预测条带与观测条带分别按长度升序排列后做一对一匹配：
仅当 |预测长度 - 观测长度| <= 匹配容差 时允许两条带配对；
目标为先最大化匹配对数，再最小化总误差（各配对长度差之和）。

一维长度上的绝对差代价满足 Monge 性质：交叉配对"解交叉"后可行性不变、
代价不增，因此最优匹配必可取为保序（不交叉）匹配，可用序列比对式动态
规划在 O(n·m) 内求得。并列时按"跳过预测条带 → 跳过观测条带 → 配对"
的固定优先级选择，等效于让较短的条带优先配对，保证结果完全确定，
且与输入列表的提交顺序无关。

判定规则
--------
分别统计两种样本的未匹配条带数（缺失 + 额外）与总误差：
未匹配条带数越少者判为对应样本；并列时总误差小者胜出；
两者完全一致时判定为 ambiguous（无法区分）。判定仅使用响应中
输出的计数值与（保留两位小数的）总误差，可由响应直接复算。
"""
from app.models import (
    BandMatchOut,
    IdentifyRequest,
    IdentifyResponse,
    SampleMatchOut,
)

# DP 转移优先级：跳过预测条带 > 跳过观测条带 > 配对（仅在严格更优时替换）
_SKIP_PREDICTED = "skip_predicted"
_SKIP_OBSERVED = "skip_observed"
_MATCH = "match"


def _better(a: tuple[int, float], b: tuple[int, float]) -> bool:
    """状态 a 是否严格优于 b：匹配对数多者优先，其次总误差小者。"""
    return a[0] > b[0] or (a[0] == b[0] and a[1] < b[1])


def match_sample(predicted: list[float], observed: list[float],
                 tolerance: float) -> SampleMatchOut:
    """把一种样本的预测条带与观测条带做一对一最小代价匹配。"""
    p = sorted(predicted)
    o = sorted(observed)
    n, m = len(p), len(o)

    # dp[i][j] = (匹配对数, 总误差)：前 i 条预测与前 j 条观测的最优状态
    dp = [[(0, 0.0)] * (m + 1) for _ in range(n + 1)]
    choice = [[None] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            best, move = dp[i - 1][j], _SKIP_PREDICTED
            candidate = dp[i][j - 1]
            if _better(candidate, best):
                best, move = candidate, _SKIP_OBSERVED
            diff = abs(p[i - 1] - o[j - 1])
            if diff <= tolerance:
                candidate = (dp[i - 1][j - 1][0] + 1, dp[i - 1][j - 1][1] + diff)
                if _better(candidate, best):
                    best, move = candidate, _MATCH
            dp[i][j] = best
            choice[i][j] = move

    # 回溯最优路径，还原逐条配对与缺失/额外条带
    matched_pairs: list[tuple[float, float]] = []
    missing: list[float] = []
    extra: list[float] = []
    i, j = n, m
    while i > 0 and j > 0:
        move = choice[i][j]
        if move == _MATCH:
            matched_pairs.append((p[i - 1], o[j - 1]))
            i, j = i - 1, j - 1
        elif move == _SKIP_OBSERVED:
            extra.append(o[j - 1])
            j -= 1
        else:
            missing.append(p[i - 1])
            i -= 1
    while i > 0:
        missing.append(p[i - 1])
        i -= 1
    while j > 0:
        extra.append(o[j - 1])
        j -= 1

    matched_pairs.reverse()
    matches = [
        BandMatchOut(
            predicted_length=pred,
            observed_length=obs,
            error=round(abs(pred - obs), 2),
            description=(
                f"预测条带 {pred:g} bp 与观测条带 {obs:g} bp 配对，"
                f"长度差 {abs(pred - obs):g} bp，在容差 {tolerance:g} bp 之内"
            ),
        )
        for pred, obs in matched_pairs
    ]
    matched_count, total_error = dp[n][m]
    return SampleMatchOut(
        predicted_band_count=n,
        matched_count=matched_count,
        missing_count=n - matched_count,
        extra_count=m - matched_count,
        total_error=round(total_error, 2),
        matches=matches,
        missing_bands=sorted(missing),
        extra_bands=sorted(extra),
    )


def decide(reference: SampleMatchOut, variant: SampleMatchOut) -> tuple[str, str]:
    """按未匹配条带数 → 总误差的顺序比较两种样本，给出判定与中文依据。"""
    ref_unmatched = reference.missing_count + reference.extra_count
    var_unmatched = variant.missing_count + variant.extra_count
    if ref_unmatched != var_unmatched:
        verdict = "reference" if ref_unmatched < var_unmatched else "variant"
        label = "参考型" if verdict == "reference" else "变异型"
        other = "变异型" if verdict == "reference" else "参考型"
        reason = (
            f"{label}未匹配条带 {min(ref_unmatched, var_unmatched)} 条，"
            f"少于{other} {max(ref_unmatched, var_unmatched)} 条"
        )
        return verdict, reason
    if reference.total_error != variant.total_error:
        verdict = "reference" if reference.total_error < variant.total_error else "variant"
        label = "参考型" if verdict == "reference" else "变异型"
        reason = (
            f"两种样本未匹配条带数相同（各 {ref_unmatched} 条），"
            f"{label}总误差更小（参考型 {reference.total_error:g} bp，"
            f"变异型 {variant.total_error:g} bp）"
        )
        return verdict, reason
    return "ambiguous", (
        f"两种样本的未匹配条带数（各 {ref_unmatched} 条）与总误差"
        f"（各 {reference.total_error:g} bp）完全一致，无法区分"
    )


def identify(req: IdentifyRequest) -> IdentifyResponse:
    observed = list(req.observed_bands)
    reference = match_sample(req.reference_bands, observed, req.band_match_tolerance)
    variant = match_sample(req.variant_bands, observed, req.band_match_tolerance)
    verdict, reason = decide(reference, variant)
    return IdentifyResponse(
        band_match_tolerance=req.band_match_tolerance,
        observed_band_count=len(observed),
        samples={"reference": reference, "variant": variant},
        verdict=verdict,
        verdict_reason=reason,
    )
