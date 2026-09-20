"""核心算法：环状质粒双酶切片段计算、条带聚类与酶对推荐。

所有函数均为纯函数，不依赖任何外部服务或全局状态。

圆环片段规则
------------
对圆环上按升序排列的切点 p0 < p1 < ... < p_{n-1}：
  * 片段 (p_i, p_{i+1})，长度 p_{i+1} - p_i
  * 末段 (p_{n-1}, p0) 跨过坐标原点，长度 p0 + L - p_{n-1}
至少需要 2 个切点才能在圆环上切出至少 2 个片段。

条带聚类规则
------------
两种样本的可见片段（长度 >= 最小可见片段）合并后按长度升序排列，
依次扫描：相邻片段长度差 > 分辨容差时断开，否则归入同一条带
（等价于以容差为阈值的单连接聚类）。一条带中两样本的片段数
不同即为"可区分条带"。
"""
from itertools import combinations

from app.models import (
    BandMemberOut,
    BandOut,
    CutPointOut,
    DifferenceOut,
    DigestRequest,
    DigestResponse,
    FragmentOut,
    NonDiscriminatingPairOut,
    PairRecommendationOut,
    RejectedPairOut,
    SampleDigestOut,
)


# --------------------------------------------------------------- 切点合并

def merge_cut_points(site_map: dict[str, dict[str, list[int]]],
                     enzyme_a: str, enzyme_b: str, sample: str) -> list[CutPointOut]:
    """合并两种酶在某样本上的切点；相同位置合并，并记录来源酶。"""
    merged: dict[int, list[str]] = {}
    for enzyme in (enzyme_a, enzyme_b):
        for pos in site_map[enzyme][sample]:
            merged.setdefault(pos, [])
            if enzyme not in merged[pos]:
                merged[pos].append(enzyme)
    return [
        CutPointOut(position=pos, enzymes=sorted(merged[pos]))
        for pos in sorted(merged)
    ]


# --------------------------------------------------------------- 环状片段

def circular_fragments(cuts: list[CutPointOut], plasmid_length: int) -> list[FragmentOut]:
    """按圆环相邻切点计算片段，末段跨过原点。"""
    fragments: list[FragmentOut] = []
    positions = [c.position for c in cuts]
    for i in range(len(positions) - 1):
        start, end = positions[i], positions[i + 1]
        fragments.append(FragmentOut(
            start=start, end=end, length=end - start, crosses_origin=False,
        ))
    # 末段：最后一个切点环绕到第一个切点
    last, first = positions[-1], positions[0]
    fragments.append(FragmentOut(
        start=last,
        end=first,
        length=first + plasmid_length - last,
        crosses_origin=True,
    ))
    return fragments


def digest_sample(site_map: dict[str, dict[str, list[int]]],
                  enzyme_a: str, enzyme_b: str, sample: str,
                  plasmid_length: int, min_visible: int) -> SampleDigestOut:
    cuts = merge_cut_points(site_map, enzyme_a, enzyme_b, sample)
    fragments = circular_fragments(cuts, plasmid_length) if len(cuts) >= 2 else []
    visible = [f for f in fragments if f.length >= min_visible]
    return SampleDigestOut(
        cuts=cuts,
        fragments=fragments,
        visible_fragment_count=len(visible),
        hidden_fragment_count=len(fragments) - len(visible),
    )


# --------------------------------------------------------------- 条带聚类

def cluster_bands(samples: dict[str, SampleDigestOut], tolerance: int,
                  min_visible: int) -> list[BandOut]:
    """把两样本的可见片段按容差聚成条带（单连接聚类）。"""
    members: list[BandMemberOut] = []
    for sample_name, digest in samples.items():
        for frag in digest.fragments:
            if frag.length >= min_visible:
                members.append(BandMemberOut(**frag.model_dump(), sample=sample_name))
    members.sort(key=lambda m: m.length)

    clusters: list[list[BandMemberOut]] = []
    for member in members:
        if clusters and member.length - clusters[-1][-1].length <= tolerance:
            clusters[-1].append(member)
        else:
            clusters.append([member])

    bands: list[BandOut] = []
    for band_id, cluster in enumerate(clusters, start=1):
        lengths = [m.length for m in cluster]
        ref_count = sum(1 for m in cluster if m.sample == "reference")
        var_count = sum(1 for m in cluster if m.sample == "variant")
        avg = round(sum(lengths) / len(lengths), 2)
        bands.append(BandOut(
            band_id=band_id,
            representative_length=avg,
            min_length=min(lengths),
            max_length=max(lengths),
            reference_count=ref_count,
            variant_count=var_count,
            distinguishable=ref_count != var_count,
            members=cluster,
        ))
    return bands


# --------------------------------------------------------------- 差异依据

_DIFF_LABEL = {"reference": "参考型", "variant": "变异型"}


def build_differences(bands: list[BandOut]) -> list[DifferenceOut]:
    differences: list[DifferenceOut] = []
    for band in bands:
        if not band.distinguishable:
            continue
        if band.variant_count == 0:
            diff_type = "reference_only"
            desc = (f"约 {band.representative_length:g} bp 条带仅参考型出现"
                    f"（参考型 {band.reference_count} 条，变异型 0 条）")
        elif band.reference_count == 0:
            diff_type = "variant_only"
            desc = (f"约 {band.representative_length:g} bp 条带仅变异型出现"
                    f"（变异型 {band.variant_count} 条，参考型 0 条）")
        else:
            diff_type = "count_mismatch"
            desc = (f"约 {band.representative_length:g} bp 条带片段数不同："
                    f"参考型 {band.reference_count} 条，变异型 {band.variant_count} 条")
        differences.append(DifferenceOut(
            band_id=band.band_id,
            representative_length=band.representative_length,
            difference_type=diff_type,
            reference_count=band.reference_count,
            variant_count=band.variant_count,
            description=desc,
        ))
    return differences


# --------------------------------------------------------------- 酶对枚举

def pair_key(a: str, b: str) -> tuple[list[str], str]:
    pair = sorted((a, b))
    return pair, f"{pair[0]}+{pair[1]}"


def recommend(req: DigestRequest) -> DigestResponse:
    site_map = {
        enzyme.name: {
            "reference": sorted(enzyme.reference_sites),
            "variant": sorted(enzyme.variant_sites),
        }
        for enzyme in req.enzymes
    }

    recommendations: list[PairRecommendationOut] = []
    rejected: list[RejectedPairOut] = []
    indistinguishable: list[NonDiscriminatingPairOut] = []

    names = [enzyme.name for enzyme in req.enzymes]
    for a, b in combinations(names, 2):
        pair, key = pair_key(a, b)
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

        ref_cuts = len(samples["reference"].cuts)
        var_cuts = len(samples["variant"].cuts)
        if ref_cuts < 2 or var_cuts < 2:
            rejected.append(RejectedPairOut(
                enzyme_pair=pair,
                pair_key=key,
                reason_code="insufficient_cuts",
                reason=f"有效切点不足：参考型 {ref_cuts} 个、变异型 {var_cuts} 个，"
                       f"任一样本少于 2 个有效切点，无法完成圆环双酶切",
                reference_cut_count=ref_cuts,
                variant_cut_count=var_cuts,
            ))
            continue

        bands = cluster_bands(
            samples, req.band_resolution_tolerance, req.min_visible_fragment,
        )
        diffs = build_differences(bands)
        if not diffs:
            indistinguishable.append(NonDiscriminatingPairOut(
                enzyme_pair=pair, pair_key=key, total_bands=len(bands),
            ))
            continue

        recommendations.append(PairRecommendationOut(
            enzyme_pair=pair,
            pair_key=key,
            samples=samples,
            bands=bands,
            total_bands=len(bands),
            distinguishable_band_count=len(diffs),
            differences=diffs,
        ))

    # 可区分条带数降序 → 总条带数升序 → 酶名（pair_key）升序
    recommendations.sort(
        key=lambda r: (-r.distinguishable_band_count, r.total_bands, r.pair_key)
    )

    total_pairs = len(list(combinations(sorted(set(names)), 2)))
    return DigestResponse(
        plasmid_length=req.plasmid_length,
        min_visible_fragment=req.min_visible_fragment,
        band_resolution_tolerance=req.band_resolution_tolerance,
        total_enzyme_pairs=total_pairs,
        recommended_count=len(recommendations),
        rejected_pair_count=len(rejected),
        non_discriminating_pairs=indistinguishable,
        non_discriminating_pair_count=len(indistinguishable),
        recommendations=recommendations,
        rejected_pairs=rejected,
    )
