"""核心算法：批量盲样（命名谱型 ↔ 命名泳道）的全局一对一最优分配。

所有函数均为纯函数，不依赖任何外部服务或全局状态。

候选边计算
----------
对每一个“谱型 p — 泳道 w”组合，复用 :func:`app.identify.match_sample`
的同一套确定性一对一条带匹配（先最大化匹配对数，再最小化总误差），
得到该边的未匹配条带总数 ``u = 缺失数 + 额外数`` 与总误差 ``e``。
仅当 ``u <= max_unmatched_bands`` 时该边**可用**；超过阈值的边从
全局分配的候选图中剔除。

全局一对一分配
--------------
在以可用边构成的二分图上求分配（每个谱型、每个泳道至多使用一次），
依次最大化 / 最小化：

  1. 已分配泳道数最大化（即最大基数匹配）；
  2. 全部分配边的未匹配条带总数最小化；
  3. 全部分配边的总误差之和最小化。

上述三级字典序目标用**最小费用最大流（逐次最短路增广，SSP）**精确求解：
源点 → 谱型 → 泳道 → 汇点，可用边容量 1、费用为该边的 (未匹配数, 误差)；
元组费用天然实现字典序比较，增广到不存在增广路时得到最大基数，
且该基数下费用字典序最小。

字典序破同分
------------
三级指标全部相同的分配仍可能不唯一（例如若干完全相同的谱型）。
此时先把泳道按名称排序，再把该泳道序上“分配到的谱型名序列”按
字典序比较；为表达“未分配”，把未分配项统一排在**所有谱型名之后**
（其字典序名次取谱型总数 P，即比任何真实谱型的 0..P-1 都大）。

该字典序可编码为基数 ``P + 1`` 的位置数：把全部泳道均未分配的序列
（每项名次都是 P）作为基准常数 C，分配一条 p→w 边使第 w 位由 P 变为
p 的名次，费用增量为 ``(名次(p) - P) * (P + 1) ** (L - 1 - w)``
（恒为非正整数）。把该增量作为边费用的第三分量并入最小费用流，
得到唯一最优解，且因名次只由名称排序决定，**结果与输入提交顺序无关**。
"""
from app.identify import match_sample
from app.models import (
    BatchAssignmentOut,
    BatchIdentifyRequest,
    BatchIdentifyResponse,
    SampleMatchOut,
)


class _Edge:
    """残量网络中的一条有向边（反向边成对存储，容量同步更新）。"""

    __slots__ = ("fr", "to", "rev", "cap", "cost", "is_forward")

    def __init__(self, fr: int, to: int, rev: int, cap: int, cost: tuple,
                 is_forward: bool):
        self.fr = fr
        self.to = to
        self.rev = rev
        self.cap = cap
        self.cost = cost
        self.is_forward = is_forward


class _MinCostFlow:
    """小规模二分图分配用的最小费用流：单位容量边，Bellman-Ford 逐次增广。

    费用为可逐分量比较的元组（元组相加按分量相加），因此允许边费用为负
    （字典序破同分增量恒 ≤ 0）；初始残量网络无环，SSP 全程不存在负环。
    邻接表按节点编号、加边顺序构造，算法行为完全确定。
    """

    def __init__(self, n: int):
        self.n = n
        self.graph: list[list[_Edge]] = [[] for _ in range(n)]

    def add_edge(self, fr: int, to: int, cost: tuple) -> None:
        forward = _Edge(fr, to, len(self.graph[to]), 1, cost, True)
        backward = _Edge(to, fr, len(self.graph[fr]), 0,
                         tuple(-c for c in cost), False)
        self.graph[fr].append(forward)
        self.graph[to].append(backward)

    def flow(self, source: int, sink: int) -> tuple[int, tuple, list[_Edge]]:
        """从 source 向 sink 逐次增广到无法再增广，返回 (流量, 总费用, 满流正向边)。"""
        total_flow = 0
        total_cost: tuple = (0, 0, 0)
        used_edges: list[_Edge] = []

        while True:
            # Bellman-Ford：费用元组允许出现负数分量
            dist: list[tuple | None] = [None] * self.n
            prev_edge: list[_Edge | None] = [None] * self.n
            dist[source] = (0, 0, 0)
            for _ in range(self.n - 1):
                updated = False
                for v in range(self.n):
                    if dist[v] is None:
                        continue
                    for edge in self.graph[v]:
                        if edge.cap == 0:
                            continue
                        candidate = tuple(
                            dist[v][k] + edge.cost[k] for k in range(3)
                        )
                        if dist[edge.to] is None or candidate < dist[edge.to]:
                            dist[edge.to] = candidate
                            prev_edge[edge.to] = edge
                            updated = True
                if not updated:
                    break

            if dist[sink] is None:
                break

            # 所有边容量均为 1，每条增广路瓶颈恒为 1
            v = sink
            while v != source:
                edge = prev_edge[v]
                assert edge is not None
                edge.cap = 0
                self.graph[edge.to][edge.rev].cap = 1
                v = edge.fr
            total_flow += 1
            total_cost = tuple(total_cost[k] + dist[sink][k] for k in range(3))

        # 容量降为 0 的原始正向边（源点边 / 谱型-泳道边 / 泳道边）即被选中
        used_edges = [
            edge
            for edges in self.graph
            for edge in edges
            if edge.is_forward and edge.cap == 0
        ]
        return total_flow, total_cost, used_edges


def resolve_batch_identify(req: BatchIdentifyRequest) -> BatchIdentifyResponse:
    """计算全部候选边、剔除超阈值边，再做全局一对一最优分配并组装响应。"""
    profile_names = sorted(p.name for p in req.profiles)
    lane_names = sorted(w.name for w in req.lanes)
    bands_by_profile = {p.name: list(p.predicted_bands) for p in req.profiles}
    bands_by_lane = {w.name: list(w.observed_bands) for w in req.lanes}
    p_count, w_count = len(profile_names), len(lane_names)
    tolerance = req.band_match_tolerance
    max_unmatched = req.max_unmatched_bands

    p_index = {name: i for i, name in enumerate(profile_names)}
    w_index = {name: j for j, name in enumerate(lane_names)}

    # 1) 复用一对一条带匹配评估每条候选边，超过未匹配阈值的边不可用
    matches: dict[tuple[str, str], SampleMatchOut] = {}
    eligible: dict[tuple[str, str], tuple[int, int]] = {}
    for pname in profile_names:
        for wname in lane_names:
            result = match_sample(
                bands_by_profile[pname], bands_by_lane[wname], tolerance,
            )
            unmatched = result.missing_count + result.extra_count
            matches[(pname, wname)] = result
            if unmatched <= max_unmatched:
                # 误差以百分之一 bp 的整数表示，避免浮点数参与费用比较
                eligible[(pname, wname)] = (unmatched, int(round(result.total_error * 100)))

    # 2) 最小费用最大流：source -> 谱型 -> 泳道 -> sink
    source = 0
    profile_base = 1
    lane_base = profile_base + p_count
    sink = lane_base + w_count
    network = _MinCostFlow(sink + 1)

    for i in range(p_count):
        network.add_edge(source, profile_base + i, (0, 0, 0))
    for j in range(w_count):
        network.add_edge(lane_base + j, sink, (0, 0, 0))

    base = p_count + 1  # 字典序编码基数：未分配名次 P 比任何谱型名次都大
    for (pname, wname), (unmatched, error_cents) in eligible.items():
        i, j = p_index[pname], w_index[wname]
        # 相对“全部泳道未分配”基准序列，该边带来的字典序增量（恒为负整数）
        tie_delta = (i - p_count) * base ** (w_count - 1 - j)
        network.add_edge(
            profile_base + i,
            lane_base + j,
            (unmatched, error_cents, tie_delta),
        )

    flow, cost, used_edges = network.flow(source, sink)

    # 3) 从满流边还原谱型 -> 泳道分配
    assignment_by_lane: dict[str, tuple[str, SampleMatchOut, int]] = {}
    for edge in used_edges:
        if profile_base <= edge.fr < lane_base and lane_base <= edge.to < sink:
            pname = profile_names[edge.fr - profile_base]
            wname = lane_names[edge.to - lane_base]
            result = matches[(pname, wname)]
            unmatched = result.missing_count + result.extra_count
            assignment_by_lane[wname] = (pname, result, unmatched)

    assigned_profiles = {pname for pname, _, _ in assignment_by_lane.values()}
    unassigned_lanes = [name for name in lane_names if name not in assignment_by_lane]
    unused_profiles = [name for name in profile_names if name not in assigned_profiles]

    assignments: list[BatchAssignmentOut] = []
    total_unmatched = 0
    total_error = 0.0
    for wname in lane_names:
        if wname not in assignment_by_lane:
            continue
        pname, result, unmatched = assignment_by_lane[wname]
        total_unmatched += unmatched
        total_error += result.total_error
        assignments.append(BatchAssignmentOut(
            lane=wname,
            profile=pname,
            unmatched_count=unmatched,
            match=result,
            description=(
                f"泳道 {wname} 分配给谱型 {pname}：未匹配条带 {unmatched} 条"
                f"（缺失 {result.missing_count} 条、额外 {result.extra_count} 条），"
                f"总误差 {result.total_error:g} bp，未超过上限 {max_unmatched} 条，"
                f"该边可参与全局一对一分配"
            ),
        ))
    total_error = round(total_error, 2)

    reason = (
        f"共评估 {p_count * w_count} 条谱型-泳道候选边（复用一对一最小代价条带匹配），"
        f"其中 {len(eligible)} 条未匹配条带数不超过上限 {max_unmatched} 而可用，"
        f"其余 {p_count * w_count - len(eligible)} 条超过阈值被剔除；"
        f"在可用边构成的二分图上做全局一对一分配，依次以已分配泳道数最大"
        f"（{flow} 个泳道）、未匹配条带总数最小（{total_unmatched} 条）、"
        f"总误差最小（{total_error:g} bp）为目标；指标仍并列时按泳道名排序后"
        f"谱型名序列字典序破同分，未分配项排在所有谱型名之后，"
        f"故结果唯一且与输入顺序无关"
    )

    return BatchIdentifyResponse(
        band_match_tolerance=tolerance,
        max_unmatched_bands=max_unmatched,
        profile_count=p_count,
        lane_count=w_count,
        evaluated_edge_count=p_count * w_count,
        eligible_edge_count=len(eligible),
        assignments=assignments,
        unassigned_lanes=unassigned_lanes,
        unused_profiles=unused_profiles,
        assigned_lane_count=flow,
        total_unmatched_count=total_unmatched,
        total_error=total_error,
        assignment_reason=reason,
    )
