"""跨字段语义校验：位点越界、重复位点、酶名称重复等。

校验失败统一抛出 :class:`DigestValidationError`，由主应用转换为
统一结构的 HTTP 422 响应（与 Pydantic 字段校验同一响应格式）。
"""
from fastapi import HTTPException

from app.config import settings
from app.models import (
    BatchIdentifyRequest,
    DigestRequest,
    IdentifyRequest,
    MixtureRequest,
    PanelRequest,
)


class DigestValidationError(HTTPException):
    def __init__(self, errors: list[dict]):
        super().__init__(status_code=422, detail={"code": "validation_failed", "errors": errors})


def _err(loc: list[str], message: str, code: str = "invalid") -> dict:
    return {"loc": loc, "msg": message, "type": code}


def validate_request(req: DigestRequest) -> DigestRequest:
    errors: list[dict] = []
    length = req.plasmid_length

    if req.min_visible_fragment >= length:
        errors.append(_err(
            ["min_visible_fragment"],
            f"最小可见片段长度（{req.min_visible_fragment}）必须小于质粒长度（{length}）",
            "min_visible_too_large",
        ))

    if len(req.enzymes) > settings.MAX_ENZYMES:
        errors.append(_err(
            ["enzymes"],
            f"酶数量（{len(req.enzymes)}）超过允许上限 {settings.MAX_ENZYMES}",
            "too_many_enzymes",
        ))

    seen_names: set[str] = set()
    for i, enzyme in enumerate(req.enzymes):
        name = enzyme.name
        if name in seen_names:
            errors.append(_err(
                ["enzymes", i, "name"],
                f"酶名称重复：{name}",
                "duplicate_enzyme_name",
            ))
        seen_names.add(name)

        for field in ("reference_sites", "variant_sites"):
            sites = getattr(enzyme, field)
            sample = "reference" if field == "reference_sites" else "variant"

            out_of_range = [s for s in sites if s < 0 or s >= length]
            if out_of_range:
                errors.append(_err(
                    ["enzymes", i, field],
                    f"酶 {name} 在{('参考型' if sample == 'reference' else '变异型')}"
                    f"中存在越界切点 {out_of_range}：合法范围为 0 ~ {length - 1}",
                    "site_out_of_range",
                ))

            unique_sites = set(sites)
            if len(unique_sites) != len(sites):
                dupes = sorted({s for s in sites if sites.count(s) > 1})
                errors.append(_err(
                    ["enzymes", i, field],
                    f"酶 {name} 在{('参考型' if sample == 'reference' else '变异型')}"
                    f"中存在重复切点 {dupes}",
                    "duplicate_site",
                ))

    if errors:
        raise DigestValidationError(errors)
    return req


def validate_identify_request(req: IdentifyRequest) -> IdentifyRequest:
    """identify 接口的语义校验：条带列表规模上限。

    观测为空、条带为负、容差非法等基础字段约束由 Pydantic 模型直接拦截，
    与此处一样输出统一结构的 422 响应。
    """
    errors: list[dict] = []
    for field in ("reference_bands", "variant_bands", "observed_bands"):
        bands = getattr(req, field)
        if len(bands) > settings.MAX_BANDS:
            errors.append(_err(
                [field],
                f"条带数量（{len(bands)}）超过允许上限 {settings.MAX_BANDS}",
                "too_many_bands",
            ))
    if errors:
        raise DigestValidationError(errors)
    return req


def validate_panel_request(req: PanelRequest) -> PanelRequest:
    """panel 接口的语义校验：位点越界/重复、变异名一致性与规模上限。

    基础字段约束（质粒长度、可见阈值、容差、酶数量下限等）由 Pydantic
    模型直接拦截，与此处一样输出统一结构的 422 响应。
    """
    errors: list[dict] = []
    length = req.plasmid_length

    if req.min_visible_fragment >= length:
        errors.append(_err(
            ["min_visible_fragment"],
            f"最小可见片段长度（{req.min_visible_fragment}）必须小于质粒长度（{length}）",
            "min_visible_too_large",
        ))

    if len(req.enzymes) > settings.MAX_ENZYMES:
        errors.append(_err(
            ["enzymes"],
            f"酶数量（{len(req.enzymes)}）超过允许上限 {settings.MAX_ENZYMES}",
            "too_many_enzymes",
        ))

    all_variants: set[str] = set()
    for enzyme in req.enzymes:
        all_variants.update(enzyme.variant_sites)

    seen_names: set[str] = set()
    for i, enzyme in enumerate(req.enzymes):
        name = enzyme.name
        if name in seen_names:
            errors.append(_err(
                ["enzymes", i, "name"],
                f"酶名称重复：{name}",
                "duplicate_enzyme_name",
            ))
        seen_names.add(name)

        sites = enzyme.reference_sites
        out_of_range = [s for s in sites if s < 0 or s >= length]
        if out_of_range:
            errors.append(_err(
                ["enzymes", i, "reference_sites"],
                f"酶 {name} 在参考型中存在越界切点 {out_of_range}：合法范围为 0 ~ {length - 1}",
                "site_out_of_range",
            ))
        unique_sites = set(sites)
        if len(unique_sites) != len(sites):
            dupes = sorted({s for s in sites if sites.count(s) > 1})
            errors.append(_err(
                ["enzymes", i, "reference_sites"],
                f"酶 {name} 在参考型中存在重复切点 {dupes}",
                "duplicate_site",
            ))

        if set(enzyme.variant_sites) != all_variants:
            errors.append(_err(
                ["enzymes", i, "variant_sites"],
                f"酶 {name} 的变异名 {sorted(enzyme.variant_sites)} 与请求中其他酶不一致："
                f"所有酶必须定义完全相同的变异名集合 {sorted(all_variants)}",
                "inconsistent_variant_names",
            ))

        for vname, vsites in enzyme.variant_sites.items():
            if not vname.strip():
                errors.append(_err(
                    ["enzymes", i, "variant_sites"],
                    f"酶 {name} 存在空白变异名",
                    "invalid_variant_name",
                ))
                continue

            out_of_range = [s for s in vsites if s < 0 or s >= length]
            if out_of_range:
                errors.append(_err(
                    ["enzymes", i, "variant_sites", vname],
                    f"酶 {name} 在变异型 {vname} 中存在越界切点 {out_of_range}："
                    f"合法范围为 0 ~ {length - 1}",
                    "site_out_of_range",
                ))

            unique_sites = set(vsites)
            if len(unique_sites) != len(vsites):
                dupes = sorted({s for s in vsites if vsites.count(s) > 1})
                errors.append(_err(
                    ["enzymes", i, "variant_sites", vname],
                    f"酶 {name} 在变异型 {vname} 中存在重复切点 {dupes}",
                    "duplicate_site",
                ))

    if not all_variants:
        errors.append(_err(
            ["enzymes"],
            "至少需要定义 1 个变异型：所有酶的 variant_sites 均为空",
            "no_variants",
        ))
    elif len(all_variants) > settings.MAX_VARIANTS:
        errors.append(_err(
            ["enzymes"],
            f"变异型数量（{len(all_variants)}）超过允许上限 {settings.MAX_VARIANTS}",
            "too_many_variants",
        ))

    if errors:
        raise DigestValidationError(errors)
    return req


def validate_batch_identify_request(req: BatchIdentifyRequest) -> BatchIdentifyRequest:
    """batch-identify 接口的语义校验：谱型/泳道规模、重名、空名与条带规模。

    谱型/泳道列表为空、条带为负、容差非法或缺失、``max_unmatched_bands``
    为负 / 缺失等基础字段约束由 Pydantic 模型直接拦截，与此处一样输出
    统一结构的 422 响应。谱型名与泳道名分属不同命名空间，允许同名。
    """
    errors: list[dict] = []

    if len(req.profiles) > settings.MAX_BATCH_ITEMS:
        errors.append(_err(
            ["profiles"],
            f"预测谱型数量（{len(req.profiles)}）超过允许上限 {settings.MAX_BATCH_ITEMS}",
            "too_many_profiles",
        ))
    if len(req.lanes) > settings.MAX_BATCH_ITEMS:
        errors.append(_err(
            ["lanes"],
            f"观测泳道数量（{len(req.lanes)}）超过允许上限 {settings.MAX_BATCH_ITEMS}",
            "too_many_lanes",
        ))

    for field_name, items in (("profiles", req.profiles), ("lanes", req.lanes)):
        seen_names: set[str] = set()
        for i, item in enumerate(items):
            name = item.name
            if not name.strip():
                errors.append(_err(
                    [field_name, i, "name"],
                    f"{'谱型' if field_name == 'profiles' else '泳道'}名不能为空白",
                    "invalid_name",
                ))
            if name in seen_names:
                errors.append(_err(
                    [field_name, i, "name"],
                    f"{'谱型' if field_name == 'profiles' else '泳道'}名重复：{name}",
                    "duplicate_name",
                ))
            seen_names.add(name)

            bands_attr = "predicted_bands" if field_name == "profiles" else "observed_bands"
            if len(getattr(item, bands_attr)) > settings.MAX_BANDS:
                errors.append(_err(
                    [field_name, i, bands_attr],
                    f"{'谱型' if field_name == 'profiles' else '泳道'} {name} 的"
                    f"条带数量（{len(getattr(item, bands_attr))}）超过允许上限 "
                    f"{settings.MAX_BANDS}",
                    "too_many_bands",
                ))

    if errors:
        raise DigestValidationError(errors)
    return req


def validate_mixture_request(req: MixtureRequest) -> MixtureRequest:
    """mixture 接口的语义校验：候选规模/重名/空候选与成分数上限。

    观测或候选条带为空、条带为负、容差非法或缺失、``max_components`` 缺失 /
    小于 1 等基础字段约束由 Pydantic 模型直接拦截，与此处一样输出统一结构的
    422 响应。
    """
    errors: list[dict] = []
    candidate_count = len(req.components)

    if candidate_count > settings.MAX_MIXTURE_CANDIDATES:
        errors.append(_err(
            ["components"],
            f"候选成分数量（{candidate_count}）超过允许上限 "
            f"{settings.MAX_MIXTURE_CANDIDATES}",
            "too_many_mixture_candidates",
        ))

    seen_names: set[str] = set()
    for i, component in enumerate(req.components):
        name = component.name
        if not name.strip():
            errors.append(_err(
                ["components", i, "name"],
                "候选成分名不能为空白",
                "invalid_component_name",
            ))
        if name in seen_names:
            errors.append(_err(
                ["components", i, "name"],
                f"候选成分名重复：{name}",
                "duplicate_component_name",
            ))
        seen_names.add(name)
        if len(component.predicted_bands) > settings.MAX_BANDS:
            errors.append(_err(
                ["components", i, "predicted_bands"],
                f"候选 {name} 的预测条带数量（{len(component.predicted_bands)}）"
                f"超过允许上限 {settings.MAX_BANDS}",
                "too_many_bands",
            ))

    if len(req.observed_bands) > settings.MAX_BANDS:
        errors.append(_err(
            ["observed_bands"],
            f"条带数量（{len(req.observed_bands)}）超过允许上限 {settings.MAX_BANDS}",
            "too_many_bands",
        ))

    # 成分数上限必须落在 [1, 候选总数]：为 0 / 负数 / 非整数由 Pydantic 拦截
    if req.max_components > candidate_count:
        errors.append(_err(
            ["max_components"],
            f"成分数上限（{req.max_components}）不能超过候选成分数量"
            f"（{candidate_count}）",
            "max_components_too_large",
        ))

    if errors:
        raise DigestValidationError(errors)
    return req
