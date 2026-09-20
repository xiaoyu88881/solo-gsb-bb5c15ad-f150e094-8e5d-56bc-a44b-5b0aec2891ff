"""核心算法：批量盲样（命名观测泳道 vs 命名预测谱型）的全局一对一分配。

所有函数均为纯函数，不依赖任何外部服务或全局状态。

两步规则
--------
1. **逐边复用一对一条带匹配**：对每个"谱型 × 泳道"组合，调用
   :func:`app.identify.match_sample` 复用鉴别接口的同一套确定性匹配
   （先最大化匹配对数，再最小化总误差），得到该边的未匹配条带数
   （缺失 + 额外）与总误差；未匹配数超过 ``max_unmatched_bands`` 的边
   **不可用**，不参与全局分配。
2. **全局一对一分配（赋权二分图匹配）**：每个谱型至多分一条泳道、
   每条泳道至多配一个谱型，依次
     ① 最大化已分配泳道数；
     ② 最小化未匹配总数（已分配边未匹配数 + 未分配泳道的全部观测条带数）；
     ③ 最小化总误差（全部已分配边匹配对长度差之和）；
     ④ 按泳道名字典序排序后，以对应谱型名序列的字典序破同分
        （未分配项排在所有谱型名之后）。

实现方式
--------
补虚拟行（虚拟泳道，吸收未使用谱型）/ 虚拟列（虚拟谱型，吸收未分配
泳道）后，问题成为方阵最小权完美匹配，用经典 O(k³) 匈牙利算法求解。

四个依次比较的目标被编码为一个**混合进制整数**：低序目标的取值范围
严格小于高序目标的进制，因此整数大小次序与四级目标的字典序完全一致，
低位相加不可能进位影响高位。进制全部按请求数据现场取上界（而非拍
脑袋的常数），保证编码恒为保序映射。最优解唯一，且只依赖排序后的
泳道名 / 谱型名，**与请求中谱型、泳道的提交顺序无关**。
"""
from app.identify import match_sample
from app.models import (
    BatchAssignmentOut,
    BatchEdgeOut,
    BatchIdentifyRequest,
    BatchIdentifyResponse,
    SampleMatchOut,
)


def evaluate_edges(
    profile_bands: dict[str, list[float]],
    lane_bands: dict[str, list[float]],
    lanes: list[str],
    profiles: list[str],
    tolerance: float,
    max_unmatched: int,
) -> dict[tuple[str, str], tuple[SampleMatchOut, bool]]:
    """逐"谱型-泳道"复用一对一条带匹配，标记边是否可用。

    返回 (泳道名, 谱型名) -> (匹配结果, 是否可用)。
    """
    edges: dict[tuple[str, str], tuple[SampleMatchOut, bool]] = {}
    for lane in lanes:
        for profile in profiles:
            match = match_sample(profile_bands[profile], lane_bands[lane], tolerance)
            unmatched = match.missing_count + match.extra_count
            edges[(lane, profile)] = (match, unmatched <= max_unmatched)
    return edges


# --------------------------------------------------------------- 匈牙利算法

def hungarian(cost: list[list[int]]) -> tuple[int, ...]:
    """方阵最小权完美匹配（匈牙利，O(k³)），非负整数代价。

    行/列下标均为 0-based；返回每行匹配的列号。
    """
    n = len(cost)
    inf = 10 ** 60
    u = [0] * (n + 1)
    v = [0] * (n + 1)
    p = [0] * (n + 1)
    way = [0] * (n + 1)

    for i in range(1, n + 1):
        p[0] = i
        j0 = 0
        minv = [inf] * (n + 1)
        used = [False] * (n + 1)
        while True:
            used[j0] = True
            i0 = p[j0]
            delta = inf
            j1 = 0
            for j in range(1, n + 1):
                if not used[j]:
                    cur = cost[i0 - 1][j - 1] - u[i0] - v[j]
                    if cur < minv[j]:
                        minv[j] = cur
                        way[j] = j0
                    if minv[j] < delta:
                        delta = minv[j]
                        j1 = j
            for j in range(n + 1):
                if used[j]:
                    u[p[j]] += delta
                    v[j] -= delta
                else:
                    minv[j] -= delta
            j0 = j1
            if p[j0] == 0:
                break
        while True:
            j1 = way[j0]
            p[j0] = p[j1]
            j0 = j1
            if j0 == 0:
                break
    # p[列] = 行 -> 反查每行匹配的列
    assignment = [0] * n
    for j in range(1, n + 1):
        assignment[p[j] - 1] = j - 1
    return tuple(assignment)


# --------------------------------------------------------------- 代价编码

def _error_hundredths(match: SampleMatchOut) -> int:
    """总误差按百分之一 bp 取整为整数，保证全局比较严格精确、无浮点并列。

    直接采用匹配结果中已保留两位小数的 ``total_error``（即响应中可见的值），
    使响应里的全局总误差恒等于逐项边总误差之和、可由响应直接复算。
    """
    return int(round(match.total_error * 100))


def _mixed_radix_weights(
    lanes: list[str],
    profiles: list[str],
    lane_bands: dict[str, list[float]],
    edges: dict[tuple[str, str], tuple[SampleMatchOut, bool]],
) -> tuple[int, int, int, list[int]]:
    """按请求数据现场计算混合进制权重。

    返回 (w0, w_unmatched, w_error, 每个真实泳道的名次权重 w_rank_i)：

      整数代价 = w0·[禁用惩罚] + w_unmatched·[未匹配总数]
                 + w_error·[总误差/0.01bp] + Σ_i w_rank_i·[泳道 i 的谱型名次]

    各进制严格大于该位在任何完美匹配中的最大可能取值，低位相加不会进位，
    故整数比较与"禁用 → 已分配数 → 未匹配总数 → 总误差 → 名次序列"的
    字典序比较完全一致。
    """
    n, m = len(lanes), len(profiles)

    # 未匹配总数上界：每条真实泳道至多贡献其全部观测条带（未分配时）
    max_unmatched_total = sum(len(lane_bands[lane]) for lane in lanes)
    # 总误差（百分之一 bp）上界：只统计可用边
    max_error_total = 0
    for lane in lanes:
        for profile in profiles:
            match, usable = edges[(lane, profile)]
            if usable:
                max_error_total += _error_hundredths(match)
    # 名次每位取值 1..m+1，故进制取 m+2
    rank_base = m + 2

    # 低位（名次）位宽 = n 位 rank_base 进制
    rank_span = rank_base ** n
    w_error = rank_span
    w_unmatched = (max_error_total + 1) * w_error
    # 禁用惩罚只需大于"全部真实泳道均未分配"的合法最大代价
    max_legal = max_unmatched_total * w_unmatched + max_error_total * w_error + rank_span
    w0 = max_legal + 1
    # 名次权重：泳道 i（按字典序越靠前位权越高）
    w_rank = [rank_base ** (n - 1 - i) for i in range(n)]
    return w0, w_unmatched, w_error, w_rank


def build_cost_matrix(
    lanes: list[str],
    profiles: list[str],
    lane_bands: dict[str, list[float]],
    edges: dict[tuple[str, str], tuple[SampleMatchOut, bool]],
) -> list[list[int]]:
    """构造 (n+m)×(n+m) 补虚拟行/列的整数代价矩阵。

    行 0..n-1 为真实泳道，行 n..n+m-1 为虚拟泳道；
    列 0..m-1 为真实谱型，列 m..m+n-1 为虚拟谱型（未分配）。
    真实泳道行的代价编码四个优化目标；虚拟泳道行代价恒为 0
    （吸收未使用谱型与剩余虚拟列）。
    """
    n, m = len(lanes), len(profiles)
    size = n + m
    w0, w_unmatched, w_error, w_rank = _mixed_radix_weights(
        lanes, profiles, lane_bands, edges,
    )

    rows: list[list[int]] = []
    for i, lane in enumerate(lanes):
        observed_count = len(lane_bands[lane])
        rank_w = w_rank[i]
        row: list[int] = []
        for j, profile in enumerate(profiles):
            match, usable = edges[(lane, profile)]
            if not usable:
                # 禁用边：高位惩罚远超任何合法完美匹配，绝不会被选中
                row.append(w0)
                continue
            unmatched = match.missing_count + match.extra_count
            row.append(
                unmatched * w_unmatched
                + _error_hundredths(match) * w_error
                + (j + 1) * rank_w
            )
        # 虚拟谱型列：该真实泳道未分配，未匹配按全部观测条带计，名次取 m+1
        unassigned_cost = (
            observed_count * w_unmatched + (m + 1) * rank_w
        )
        row.extend([unassigned_cost] * n)
        rows.append(row)

    # 虚拟泳道行：分配真实谱型（未使用谱型）或虚拟谱型，代价全为零
    rows.extend([[0] * size for _ in range(m)])
    return rows


# --------------------------------------------------------------- 分配与组装

def assign_lanes(req: BatchIdentifyRequest) -> BatchIdentifyResponse:
    """计算全部边匹配依据并求全局一对一最优分配。"""
    lanes = sorted(lane.name for lane in req.lanes)
    profiles = sorted(profile.name for profile in req.profiles)
    lane_bands = {lane.name: list(lane.observed_bands) for lane in req.lanes}
    profile_bands = {p.name: list(p.predicted_bands) for p in req.profiles}
    tolerance = req.band_match_tolerance
    max_unmatched = req.max_unmatched_bands

    edges = evaluate_edges(
        profile_bands, lane_bands, lanes, profiles, tolerance, max_unmatched,
    )

    n, m = len(lanes), len(profiles)
    cost = build_cost_matrix(lanes, profiles, lane_bands, edges)
    row_cols = hungarian(cost)

    assignment_by_lane: dict[str, str | None] = {}
    for i, lane in enumerate(lanes):
        col = row_cols[i]
        assignment_by_lane[lane] = profiles[col] if col < m else None

    # 输出全部边的逐项依据（按泳道名、谱型名排序）
    edge_outs: list[BatchEdgeOut] = []
    for lane in lanes:
        for profile in profiles:
            match, usable = edges[(lane, profile)]
            unmatched = match.missing_count + match.extra_count
            if usable:
                reason = (
                    f"谱型 {profile} 与泳道 {lane} 一对一匹配："
                    f"匹配 {match.matched_count} 对、未匹配 {unmatched} 条"
                    f"（缺失 {match.missing_count} 条、额外 {match.extra_count} 条），"
                    f"未超过阈值 {max_unmatched} 条，该边可用"
                )
            else:
                reason = (
                    f"谱型 {profile} 与泳道 {lane} 一对一匹配："
                    f"匹配 {match.matched_count} 对、未匹配 {unmatched} 条"
                    f"（缺失 {match.missing_count} 条、额外 {match.extra_count} 条），"
                    f"超过阈值 {max_unmatched} 条，该边不可用"
                )
            edge_outs.append(BatchEdgeOut(
                lane=lane,
                profile=profile,
                usable=usable,
                unmatched_count=unmatched,
                missing_count=match.missing_count,
                extra_count=match.extra_count,
                total_error=match.total_error,
                match=match if usable else None,
                reason=reason,
            ))

    assignments: list[BatchAssignmentOut] = []
    unassigned_lanes: list[str] = []
    used_profiles: set[str] = set()
    total_unmatched = 0
    total_error_hundredths = 0

    for lane in lanes:
        profile = assignment_by_lane[lane]
        if profile is None:
            unassigned_lanes.append(lane)
            observed_count = len(lane_bands[lane])
            total_unmatched += observed_count
            usable_profiles = [p for p in profiles if edges[(lane, p)][1]]
            if usable_profiles:
                description = (
                    f"泳道 {lane} 存在可用谱型 {usable_profiles}，"
                    f"但全局一对一分配中这些谱型需优先用于其他泳道，"
                    f"故该泳道未分配；其 {observed_count} 条观测条带全部计入未匹配"
                )
            else:
                description = (
                    f"泳道 {lane} 与所有谱型的未匹配条带数均超过阈值 "
                    f"{max_unmatched} 条，无可用边，故未分配；"
                    f"其 {observed_count} 条观测条带全部计入未匹配"
                )
            assignments.append(BatchAssignmentOut(
                lane=lane,
                assigned_profile=None,
                usable=False,
                unmatched_count=observed_count,
                total_error=0.0,
                match=None,
                description=description,
            ))
            continue

        used_profiles.add(profile)
        match, _ = edges[(lane, profile)]
        unmatched = match.missing_count + match.extra_count
        total_unmatched += unmatched
        total_error_hundredths += _error_hundredths(match)
        description = (
            f"泳道 {lane} 分配给谱型 {profile}：匹配 {match.matched_count} 对、"
            f"未匹配 {unmatched} 条（缺失 {match.missing_count} 条、"
            f"额外 {match.extra_count} 条）、总误差 {match.total_error:g} bp"
        )
        assignments.append(BatchAssignmentOut(
            lane=lane,
            assigned_profile=profile,
            usable=True,
            unmatched_count=unmatched,
            total_error=match.total_error,
            match=match,
            description=description,
        ))

    unused_profiles = [p for p in profiles if p not in used_profiles]
    total_error = round(total_error_hundredths / 100, 2)
    assigned_count = n - len(unassigned_lanes)

    reason = (
        f"先对全部 {n * m} 条谱型-泳道边复用一对一条带匹配"
        f"（容差 {tolerance:g} bp），未匹配条带数超过 {max_unmatched} 的边不可用；"
        f"再在可用边上做全局一对一分配，依次以 ① 已分配泳道数最多"
        f"（本次 {assigned_count}/{n} 条）、② 未匹配总数最少"
        f"（本次 {total_unmatched} 条）、③ 总误差最小"
        f"（本次 {total_error:g} bp）为目标；并列时按泳道名字典序排序后的"
        f"谱型名序列字典序破同分（未分配排在所有谱型名之后），"
        f"结果唯一且与输入顺序无关"
    )

    return BatchIdentifyResponse(
        band_match_tolerance=tolerance,
        max_unmatched_bands=max_unmatched,
        profile_count=m,
        lane_count=n,
        assigned_lane_count=assigned_count,
        total_unmatched_count=total_unmatched,
        total_error=total_error,
        assignments=assignments,
        unassigned_lanes=unassigned_lanes,
        unused_profiles=unused_profiles,
        edges=edge_outs,
        assignment_reason=reason,
    )
