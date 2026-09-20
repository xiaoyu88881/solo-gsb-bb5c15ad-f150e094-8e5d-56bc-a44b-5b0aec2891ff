"""请求 / 响应数据模型。

输入使用 Pydantic 做基础类型约束；跨字段的语义校验
（位点越界、重复位点、切点数量等）在 :mod:`app.validation` 中完成，
两者最终都输出统一结构的 422 响应。
"""
from typing import Annotated, Literal

from pydantic import BaseModel, Field

SampleName = Literal["reference", "variant"]

# 条带长度 / 匹配容差：非负有限浮点数（JSON 整数自动转换）
NonNegativeBand = Annotated[float, Field(ge=0, allow_inf_nan=False)]


# ---------------------------------------------------------------- 输入模型

class EnzymeInput(BaseModel):
    name: str = Field(..., min_length=1, description="酶名称，如 EcoRI")
    reference_sites: list[int] = Field(
        default_factory=list, description="参考型上该酶的候选切点位置（bp，0-based）"
    )
    variant_sites: list[int] = Field(
        default_factory=list, description="变异型上该酶的候选切点位置（bp，0-based）"
    )


class DigestRequest(BaseModel):
    plasmid_length: int = Field(..., gt=0, description="环状质粒总长度（bp）")
    min_visible_fragment: int = Field(
        ..., ge=0, description="最小可见片段长度（bp），短于此值的片段不参与条带聚类"
    )
    band_resolution_tolerance: int = Field(
        ..., ge=0, description="条带分辨容差（bp），长度差不超过该值的相邻片段聚为同一条带"
    )
    enzymes: list[EnzymeInput] = Field(..., min_length=2, description="参与枚举的酶，至少 2 种")


# ---------------------------------------------------------------- 输出模型

class CutPointOut(BaseModel):
    position: int
    enzymes: list[str] = Field(description="在该位置产生切点的酶（可能两酶同位点）")


class FragmentOut(BaseModel):
    start: int = Field(description="片段起始切点位置")
    end: int = Field(description="片段终止切点位置（跨原点时为环绕后的位置）")
    length: int
    crosses_origin: bool = Field(description="是否为跨过坐标原点（0 点）的片段")


class SampleDigestOut(BaseModel):
    cuts: list[CutPointOut]
    fragments: list[FragmentOut]
    visible_fragment_count: int
    hidden_fragment_count: int = Field(description="短于最小可见片段、未参与成带的片段数")


class BandMemberOut(FragmentOut):
    sample: SampleName


class BandOut(BaseModel):
    band_id: int = Field(description="条带编号，按代表长度升序，从 1 开始")
    representative_length: float = Field(description="条带内片段长度平均值（bp）")
    min_length: int
    max_length: int
    reference_count: int
    variant_count: int
    distinguishable: bool = Field(description="两种样本在该条带中的片段数是否不同")
    members: list[BandMemberOut]


class DifferenceOut(BaseModel):
    band_id: int
    representative_length: float
    difference_type: Literal["reference_only", "variant_only", "count_mismatch"]
    reference_count: int
    variant_count: int
    description: str = Field(description="中文差异依据说明")


class PairRecommendationOut(BaseModel):
    enzyme_pair: list[str]
    pair_key: str
    samples: dict[SampleName, SampleDigestOut]
    bands: list[BandOut]
    total_bands: int
    distinguishable_band_count: int
    differences: list[DifferenceOut]


class RejectedPairOut(BaseModel):
    enzyme_pair: list[str]
    pair_key: str
    reason_code: Literal["insufficient_cuts"]
    reason: str
    reference_cut_count: int
    variant_cut_count: int


class NonDiscriminatingPairOut(BaseModel):
    enzyme_pair: list[str]
    pair_key: str
    total_bands: int = Field(description="两种样本条带完全一致时的条带总数")


class DigestResponse(BaseModel):
    plasmid_length: int
    min_visible_fragment: int
    band_resolution_tolerance: int
    total_enzyme_pairs: int
    recommended_count: int
    rejected_pair_count: int
    non_discriminating_pair_count: int
    recommendations: list[PairRecommendationOut]
    rejected_pairs: list[RejectedPairOut]
    non_discriminating_pairs: list[NonDiscriminatingPairOut]


class HealthResponse(BaseModel):
    status: Literal["ok"]
    service: str
    version: str


# ---------------------------------------------------------------- 鉴别（identify）模型

class IdentifyRequest(BaseModel):
    reference_bands: list[NonNegativeBand] = Field(
        default_factory=list, description="参考型预测条带长度列表（bp）"
    )
    variant_bands: list[NonNegativeBand] = Field(
        default_factory=list, description="变异型预测条带长度列表（bp）"
    )
    observed_bands: list[NonNegativeBand] = Field(
        ..., min_length=1, description="一次凝胶观测到的条带长度列表（bp），至少 1 条"
    )
    band_match_tolerance: NonNegativeBand = Field(
        ..., description="条带匹配容差（bp）：长度差不超过该值的预测/观测条带才允许配对"
    )


class BandMatchOut(BaseModel):
    predicted_length: float = Field(description="预测条带长度（bp）")
    observed_length: float = Field(description="与之配对的观测条带长度（bp）")
    error: float = Field(description="预测与观测长度差的绝对值（bp）")
    description: str = Field(description="中文匹配依据说明")


class SampleMatchOut(BaseModel):
    predicted_band_count: int
    matched_count: int
    missing_count: int = Field(description="缺失条带数：未在观测中找到配对的预测条带")
    extra_count: int = Field(description="额外条带数：未匹配到任何预测条带的观测条带")
    total_error: float = Field(description="全部匹配对的长度差之和（bp）")
    matches: list[BandMatchOut] = Field(description="逐条匹配依据，按预测长度升序")
    missing_bands: list[float] = Field(description="缺失条带：未匹配到观测的预测长度，升序")
    extra_bands: list[float] = Field(description="额外条带：未匹配到预测的观测长度，升序")


class IdentifyResponse(BaseModel):
    band_match_tolerance: float
    observed_band_count: int
    samples: dict[SampleName, SampleMatchOut]
    verdict: Literal["reference", "variant", "ambiguous"]
    verdict_reason: str = Field(description="中文判定依据说明")


# ---------------------------------------------------------------- 面板（panel）模型

class PanelEnzymeInput(BaseModel):
    name: str = Field(..., min_length=1, description="酶名称，如 EcoRI")
    reference_sites: list[int] = Field(
        default_factory=list, description="参考型上该酶的候选切点位置（bp，0-based）"
    )
    variant_sites: dict[str, list[int]] = Field(
        default_factory=dict,
        description="按变异名映射的候选切点：变异名 -> 该变异型上的切点位置列表（bp，0-based）",
    )


class PanelRequest(BaseModel):
    plasmid_length: int = Field(..., gt=0, description="环状质粒总长度（bp）")
    min_visible_fragment: int = Field(
        ..., ge=0, description="最小可见片段长度（bp），短于此值的片段不参与条带聚类"
    )
    band_resolution_tolerance: int = Field(
        ..., ge=0, description="条带分辨容差（bp），长度差不超过该值的相邻片段聚为同一条带"
    )
    enzymes: list[PanelEnzymeInput] = Field(..., min_length=2, description="参与枚举的酶，至少 2 种")


class PanelVariantCoverageOut(BaseModel):
    variant: str = Field(description="变异名")
    distinguishable_band_count: int = Field(description="该酶对在此变异型上的可区分条带数")
    differences: list[DifferenceOut] = Field(description="逐条差异条带依据（与推荐接口同结构）")


class PanelSelectedPairOut(BaseModel):
    enzyme_pair: list[str]
    pair_key: str
    covered_variants: list[str] = Field(description="该酶对覆盖的变异名，按字典序")
    distinguishable_band_count: int = Field(description="该酶对在全部变异型上的可区分条带总数")
    coverages: list[PanelVariantCoverageOut] = Field(description="逐变异型的差异条带依据")


class PanelResponse(BaseModel):
    plasmid_length: int
    min_visible_fragment: int
    band_resolution_tolerance: int
    variant_count: int = Field(description="请求中全部变异型的数量")
    total_enzyme_pairs: int
    selected_pairs: list[PanelSelectedPairOut] = Field(
        description="最小面板中的酶对，按 pair_key 字典序；存在未覆盖变异型时为空（不返回部分面板）"
    )
    covered_variants: list[str] = Field(description="面板覆盖的变异名；未给出面板时为空")
    uncovered_variants: list[str] = Field(
        description="无任何可用酶对覆盖的变异名；非空时不返回部分面板"
    )
    total_distinguishable_bands: int = Field(description="面板全部酶对的可区分条带总数")


# ---------------------------------------------------------------- 混合（mixture）模型

class MixtureComponentInput(BaseModel):
    name: str = Field(..., min_length=1, description="候选成分名（如某变异型 / 单菌编号）")
    predicted_bands: list[NonNegativeBand] = Field(
        ..., min_length=1, description="该候选成分单独的预测条带长度列表（bp），至少 1 条"
    )


class MixtureRequest(BaseModel):
    components: list[MixtureComponentInput] = Field(
        ..., min_length=1, description="命名候选成分列表，名称不可重复，至少 1 个候选"
    )
    observed_bands: list[NonNegativeBand] = Field(
        ..., min_length=1, description="一次凝胶观测到的混合样本条带长度列表（bp），至少 1 条"
    )
    band_match_tolerance: NonNegativeBand = Field(
        ...,
        description="条带容差（bp）：既用于候选间共迁移条带的单连接聚带，"
                    "也用于合并预测条带与观测条带的一对一配对",
    )
    max_components: int = Field(
        ..., ge=1, description="成分数上限：枚举成分数不超过该值的全部候选子集，且不得超过候选总数"
    )


class MixtureResponse(BaseModel):
    band_match_tolerance: float
    observed_band_count: int
    max_components: int
    selected_components: list[str] = Field(
        description="选中的候选成分名（成分数不超过上限的最优子集），按字典序"
    )
    merged_predicted_bands: list[float] = Field(
        description="选中子集全部预测条带按容差单连接聚带后的代表长度（簇内均值），升序"
    )
    match: SampleMatchOut = Field(
        description="合并预测条带与观测条带的一对一最小代价匹配结果（与鉴别接口同结构）"
    )
    evaluated_subset_count: int = Field(
        description="实际评估的候选子集数量（成分数 1..max_components 的全部非空子集）"
    )
    selection_reason: str = Field(
        description="中文选定依据：未匹配条带数、总误差、成分数及候选名列表如何决定唯一结果"
    )


# ------------------------------------------------- 批量盲样分配（batch-identify）模型

class BatchProfileInput(BaseModel):
    name: str = Field(..., min_length=1, description="命名预测谱型名称（如菌株 / 单菌编号）")
    predicted_bands: list[NonNegativeBand] = Field(
        ..., min_length=1, description="该谱型的预测条带长度列表（bp），至少 1 条"
    )


class BatchLaneInput(BaseModel):
    name: str = Field(..., min_length=1, description="命名观测泳道名称（如凝胶道编号）")
    observed_bands: list[NonNegativeBand] = Field(
        ..., min_length=1, description="该泳道的凝胶观测条带长度列表（bp），至少 1 条"
    )


class BatchIdentifyRequest(BaseModel):
    profiles: list[BatchProfileInput] = Field(
        ..., min_length=1, description="命名预测谱型列表，名称不可重复，至少 1 个"
    )
    lanes: list[BatchLaneInput] = Field(
        ..., min_length=1, description="命名观测泳道列表，名称不可重复，至少 1 个"
    )
    band_match_tolerance: NonNegativeBand = Field(
        ..., description="条带匹配容差（bp）：长度差不超过该值的预测/观测条带才允许配对"
    )
    max_unmatched_bands: int = Field(
        ..., ge=0,
        description="单条分配允许的最大未匹配条带数：未匹配数超过该值的谱型-泳道边不可用",
    )


class BatchAssignmentOut(BaseModel):
    lane: str = Field(description="观测泳道名；assignments 按泳道名字典序排列")
    profile: str = Field(description="全局最优分配中分配给该泳道的预测谱型名")
    unmatched_count: int = Field(
        description="该谱型-泳道边的未匹配条带总数（缺失 + 额外），不超过请求上限"
    )
    match: SampleMatchOut = Field(
        description="该谱型预测条带与该泳道观测条带的一对一最小代价匹配结果"
        "（与鉴别接口同结构，含逐条配对、缺失/额外条带依据）"
    )
    description: str = Field(description="中文分配依据说明")


class BatchIdentifyResponse(BaseModel):
    band_match_tolerance: float
    max_unmatched_bands: int
    profile_count: int
    lane_count: int
    evaluated_edge_count: int = Field(
        description="评估的谱型-泳道候选边总数，等于谱型数 × 泳道数"
    )
    eligible_edge_count: int = Field(
        description="未匹配条带数不超过上限、可参与全局一对一分配的候选边数"
    )
    assignments: list[BatchAssignmentOut] = Field(
        description="逐项匹配依据，按泳道名字典序排列；每个谱型与每个泳道至多出现一次"
    )
    unassigned_lanes: list[str] = Field(
        description="未分配到任何谱型的泳道名，按字典序"
    )
    unused_profiles: list[str] = Field(
        description="未分配给任何泳道的谱型名，按字典序"
    )
    assigned_lane_count: int = Field(description="已分配泳道数（全局最大可分配数）")
    total_unmatched_count: int = Field(description="全部分配边的未匹配条带总数（同级最优中最小）")
    total_error: float = Field(description="全部分配边的总误差之和（bp，同级最优中最小）")
    assignment_reason: str = Field(
        description="中文分配依据：可用边阈值、三级优化目标与字典序破同分如何决定唯一结果"
    )
