"""核心算法：覆盖全部变异型的最小酶对面板（精确最小集合覆盖，非贪心近似）。

所有函数均为纯函数，不依赖任何外部服务或全局状态。

面板规则
--------
对每一对不同酶，逐变异型复用 :mod:`app.digest` 的圆环双酶切与条带聚类规则：
某变异型相对参考型存在至少一条可区分条带时，称该酶对**覆盖**此变异型
（任一样本有效切点不足 2 个的酶对对该变异型不可用，视为不覆盖）。

精确最小面板
------------
把每个酶对的覆盖变异型集合表示为位掩码，按位或 BFS 逐层扩展：第 k 层为
恰好 k 个酶对可达的全部覆盖掩码，首次到达全掩码的层即为最小面板规模，
等价于最小集合覆盖的精确解，不使用贪心近似。

同一掩码的支配规则（保证 BFS 精确且完备）：
  * 最小规模面板中任意两个酶对的覆盖集合互不相同（否则移除其一仍全覆盖，
    与最小性矛盾），因此每个掩码只需保留一个代表酶对：可区分条带总数
    高者优先，并列取 pair_key 小者；
  * 同一掩码若已在更浅层到达，则深层再次出现时被支配（扩展只取决于掩码，
    浅层方案的后续选择完全相同且规模更小）；
  * 同一掩码在同一层只保留最优取值（条带总数高者优先，其次排序后的
    pair_key 列表字典序小者）：等长有序列表与同一集合的并集保持字典序，
    因此被支配取值不可能经扩展反超。

同规模方案依次按可区分条带总数降序、排序后的 pair_key 列表字典序选定。
"""
from dataclasses import dataclass, field
from itertools import combinations

from app.digest import build_differences, cluster_bands, digest_sample, pair_key
from app.models import (
    DifferenceOut,
    PanelRequest,
    PanelResponse,
    PanelSelectedPairOut,
    PanelVariantCoverageOut,
)


@dataclass
class _PairEval:
    """一个酶对在所有变异型上的评估结果。"""

    enzyme_pair: list[str]
    pair_key: str
    mask: int = 0                                   # 覆盖变异型的位掩码
    bands_total: int = 0                            # 全部变异型上的可区分条带总数
    differences: dict[str, list[DifferenceOut]] = field(default_factory=dict)


# --------------------------------------------------------------- 酶对评估

def evaluate_pairs(req: PanelRequest, variants: list[str]) -> list[_PairEval]:
    """逐酶对、逐变异型计算覆盖掩码与差异条带依据。"""
    vindex = {v: i for i, v in enumerate(variants)}
    ref_sites = {e.name: sorted(e.reference_sites) for e in req.enzymes}
    var_sites = {e.name: dict(e.variant_sites) for e in req.enzymes}
    names = [e.name for e in req.enzymes]

    evaluations: list[_PairEval] = []
    for a, b in combinations(names, 2):
        pair, key = pair_key(a, b)
        ev = _PairEval(enzyme_pair=pair, pair_key=key)
        for v in variants:
            site_map = {
                a: {"reference": ref_sites[a], "variant": sorted(var_sites[a][v])},
                b: {"reference": ref_sites[b], "variant": sorted(var_sites[b][v])},
            }
            samples = {
                "reference": digest_sample(
                    site_map, a, b, "reference",
                    req.plasmid_length, req.min_visible_fragment,
                ),
                "variant": digest_sample(
                    site_map, a, b, "variant",
                    req.plasmid_length, req.min_visible_fragment,
                ),
            }
            # 任一样本有效切点不足 2 个：无法完成圆环双酶切，对该变异型不可用
            if len(samples["reference"].cuts) < 2 or len(samples["variant"].cuts) < 2:
                continue
            bands = cluster_bands(
                samples, req.band_resolution_tolerance, req.min_visible_fragment,
            )
            diffs = build_differences(bands)
            if not diffs:
                continue
            ev.mask |= 1 << vindex[v]
            ev.bands_total += len(diffs)
            ev.differences[v] = diffs
        evaluations.append(ev)
    return evaluations


# --------------------------------------------------------------- 精确最小面板

def _better_value(a: tuple[int, tuple[str, ...]], b: tuple[int, tuple[str, ...]]) -> bool:
    """面板取值 a 是否优于 b：可区分条带总数多者优先，其次 pair_key 列表字典序小者。"""
    return a[0] > b[0] or (a[0] == b[0] and a[1] < b[1])


def exact_min_panel(entries: list[tuple[int, int, str]],
                    variant_count: int) -> tuple[int, tuple[str, ...]] | None:
    """精确求覆盖全部变异型的最小酶对面板。

    entries 为 (覆盖掩码, 可区分条带总数, pair_key)，掩码为 0 的酶对不参与；
    同一掩码的多个酶对在此去重，只保留最优代表（条带总数高者优先，并列取
    pair_key 小者）。返回 (可区分条带总数, 排序后的 pair_key 元组)；
    无法覆盖全部变异型时返回 None。
    """
    best_by_mask: dict[int, tuple[int, str]] = {}
    for mask, bands, key in entries:
        if mask == 0:
            continue
        cur = best_by_mask.get(mask)
        if cur is None or bands > cur[0] or (bands == cur[0] and key < cur[1]):
            best_by_mask[mask] = (bands, key)
    representatives = [(m, bands, key) for m, (bands, key) in best_by_mask.items()]

    full = (1 << variant_count) - 1
    visited = {0}
    # 每层的掩码 -> (可区分条带总数, 排序后的 pair_key 元组)
    current: dict[int, tuple[int, tuple[str, ...]]] = {0: (0, ())}
    while current:
        nxt: dict[int, tuple[int, tuple[str, ...]]] = {}
        best_full: tuple[int, tuple[str, ...]] | None = None
        for mask, (bands, keys) in current.items():
            for emask, ebands, ekey in representatives:
                new = mask | emask
                if new == mask or new in visited:
                    continue
                candidate = (bands + ebands, tuple(sorted(keys + (ekey,))))
                if new == full:
                    if best_full is None or _better_value(candidate, best_full):
                        best_full = candidate
                    continue
                old = nxt.get(new)
                if old is None or _better_value(candidate, old):
                    nxt[new] = candidate
        if best_full is not None:
            return best_full
        visited.update(nxt)
        current = nxt
    return None


# --------------------------------------------------------------- 面板组装

def build_panel(req: PanelRequest) -> PanelResponse:
    variants = sorted({v for e in req.enzymes for v in e.variant_sites})
    evaluations = evaluate_pairs(req, variants)
    by_key = {ev.pair_key: ev for ev in evaluations}

    union = 0
    for ev in evaluations:
        union |= ev.mask
    full = (1 << len(variants)) - 1

    def _response(selected: list[PanelSelectedPairOut], covered: list[str],
                  uncovered: list[str], total_bands: int) -> PanelResponse:
        return PanelResponse(
            plasmid_length=req.plasmid_length,
            min_visible_fragment=req.min_visible_fragment,
            band_resolution_tolerance=req.band_resolution_tolerance,
            variant_count=len(variants),
            total_enzyme_pairs=len(evaluations),
            selected_pairs=selected,
            covered_variants=covered,
            uncovered_variants=uncovered,
            total_distinguishable_bands=total_bands,
        )

    if union != full:
        # 存在任何酶对都无法覆盖的变异型：不返回部分面板
        uncovered = [v for i, v in enumerate(variants) if not (union >> i) & 1]
        return _response([], [], uncovered, 0)

    result = exact_min_panel(
        [(ev.mask, ev.bands_total, ev.pair_key) for ev in evaluations],
        len(variants),
    )
    total_bands, keys = result  # union == full，精确解必然存在

    selected: list[PanelSelectedPairOut] = []
    for key in keys:  # keys 已按字典序排序
        ev = by_key[key]
        covered = sorted(ev.differences)
        selected.append(PanelSelectedPairOut(
            enzyme_pair=ev.enzyme_pair,
            pair_key=ev.pair_key,
            covered_variants=covered,
            distinguishable_band_count=ev.bands_total,
            coverages=[
                PanelVariantCoverageOut(
                    variant=v,
                    distinguishable_band_count=len(ev.differences[v]),
                    differences=ev.differences[v],
                )
                for v in covered
            ],
        ))
    return _response(selected, variants, [], total_bands)
