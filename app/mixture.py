"""核心算法：混合样本（多成分共迁移）的候选子集枚举、聚带与最优子集选定。

所有函数均为纯函数，不依赖任何外部服务或全局状态。

混合解析规则
------------
请求给出若干**命名候选**及其单独的预测条带、一次凝胶观测条带、匹配容差
与成分数上限 ``max_components``。对每个成分数不超过上限的非空候选子集：

1. **共迁移聚带（单连接）**：把子集内各候选的预测条带合并为多重集后按
   长度升序扫描，相邻条带长度差 ≤ 容差即归入同一簇（等价于以容差为阈值
   的单连接聚类，可传递成链）。共迁移条带只出现一次，簇的代表长度取簇内
   全部条带长度的均值。
2. **一对一最小代价匹配**：把合并预测条带与观测条带交给
   :func:`app.identify.match_sample` 复用鉴别接口的同一套确定性匹配
   （先最大化匹配对数，再最小化总误差），得到逐条配对、缺失/额外条带。

最优子集唯一选定
----------------
全部评估子集依次按以下次序比较，全部相同时结果仍一致（候选名列表字典序
为最终次序），因此选定结果唯一且与候选提交顺序无关：

  1. 未匹配条带总数（缺失 + 额外）升序；
  2. 总误差（各配对长度差之和）升序；
  3. 成分数升序；
  4. 排序后的候选名列表字典序升序。
"""
from itertools import combinations

from app.identify import match_sample
from app.models import MixtureRequest, MixtureResponse, SampleMatchOut


# --------------------------------------------------------------- 共迁移聚带

def cluster_predicted_bands(bands: list[float], tolerance: float) -> list[float]:
    """把多个候选的预测条带按容差做单连接聚带，返回各簇代表长度（簇内均值）。

    共迁移条带（含同一候选内重复或多候选共迁移）只出现一次；相邻条带长度差
    ≤ 容差即同簇，链式可传递。代表长度取簇内均值并保留两位小数，保证响应中
    输出的合并条带与参与匹配的条带完全一致、可由响应直接复算。
    """
    ordered = sorted(bands)
    clusters: list[list[float]] = []
    for length in ordered:
        if clusters and length - clusters[-1][-1] <= tolerance:
            clusters[-1].append(length)
        else:
            clusters.append([length])
    return [round(sum(cluster) / len(cluster), 2) for cluster in clusters]


# --------------------------------------------------------------- 子集评估

def evaluate_subset(names: tuple[str, ...], bands_by_name: dict[str, list[float]],
                    observed: list[float],
                    tolerance: float) -> tuple[list[float], SampleMatchOut]:
    """合并子集内候选的预测条带、聚带后复用一对一最小代价匹配。

    返回 (合并预测条带代表长度列表, 匹配结果)。
    """
    pooled: list[float] = []
    for name in names:
        pooled.extend(bands_by_name[name])
    merged = cluster_predicted_bands(pooled, tolerance)
    return merged, match_sample(merged, observed, tolerance)


def resolve_mixture(req: MixtureRequest) -> MixtureResponse:
    """枚举不超过成分数上限的全部非空候选子集，按固定次序选定唯一最优子集。"""
    names = sorted(c.name for c in req.components)
    bands_by_name = {c.name: list(c.predicted_bands) for c in req.components}
    observed = list(req.observed_bands)
    tolerance = req.band_match_tolerance

    best_key: tuple[int, float, int, tuple[str, ...]] | None = None
    best_names: tuple[str, ...] = ()
    best_match: SampleMatchOut | None = None
    best_merged: list[float] = []
    evaluated = 0

    for size in range(1, req.max_components + 1):
        for subset in combinations(names, size):
            evaluated += 1
            merged, match = evaluate_subset(subset, bands_by_name, observed, tolerance)
            unmatched = match.missing_count + match.extra_count
            key = (unmatched, match.total_error, size, subset)
            if best_key is None or key < best_key:
                best_key = key
                best_names = subset
                best_match = match
                best_merged = merged

    unmatched, total_error, size, _ = best_key
    reason = (
        f"在 {evaluated} 个成分数不超过 {req.max_components} 的候选子集中，"
        f"选定 {list(best_names)}（{size} 个成分）：未匹配条带 {unmatched} 条、"
        f"总误差 {total_error:g} bp、成分数 {size}，依次按未匹配条带数升序、"
        f"总误差升序、成分数升序、候选名列表字典序比较后为唯一最优子集"
    )
    return MixtureResponse(
        band_match_tolerance=tolerance,
        observed_band_count=len(observed),
        max_components=req.max_components,
        selected_components=list(best_names),
        merged_predicted_bands=best_merged,
        match=best_match,
        evaluated_subset_count=evaluated,
        selection_reason=reason,
    )
