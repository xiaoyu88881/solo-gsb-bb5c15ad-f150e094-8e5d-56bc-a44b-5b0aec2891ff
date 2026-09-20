"""FastAPI 应用入口：环状 DNA 双酶切鉴别服务。

核心功能为纯本地计算，不依赖数据库、缓存或任何外部服务。
"""
from fastapi import FastAPI, Request
from fastapi.exceptions import RequestValidationError
from fastapi.responses import JSONResponse

from app.config import settings
from app.digest import recommend
from app.identify import identify
from app.mixture import resolve_mixture
from app.models import (
    DigestRequest,
    DigestResponse,
    HealthResponse,
    IdentifyRequest,
    IdentifyResponse,
    MixtureRequest,
    MixtureResponse,
    PanelRequest,
    PanelResponse,
)
from app.panel import build_panel
from app.validation import (
    DigestValidationError,
    validate_identify_request,
    validate_mixture_request,
    validate_panel_request,
    validate_request,
)

VERSION = "1.3.0"

app = FastAPI(
    title="Circular DNA Double-Digest Discrimination Service",
    description=(
        "圆形（环状）DNA 双酶切鉴别服务：枚举不同酶对，按圆环相邻切点计算"
        "参考型与变异型的片段长度，按分辨容差聚类条带，并推荐可区分两种样本的酶对；"
        "也可将一次凝胶观测条带与两种样本的预测条带做一对一最小代价匹配，"
        "判定观测样本为参考型、变异型或无法区分；还可针对多个变异型精确求出"
        "覆盖全部变异型的最小酶对面板（非贪心近似）；并可枚举命名候选的"
        "成分子集，先按容差单连接合并共迁移预测条带再做一对一匹配，"
        "唯一选定与一次混合样本观测最吻合的成分组合。"
    ),
    version=VERSION,
)


def _validation_response(errors: list[dict]) -> JSONResponse:
    return JSONResponse(
        status_code=422,
        content={"code": "validation_failed", "errors": errors},
    )


@app.exception_handler(DigestValidationError)
async def digest_validation_handler(request: Request, exc: DigestValidationError) -> JSONResponse:
    return _validation_response(exc.detail["errors"])


@app.exception_handler(RequestValidationError)
async def request_validation_handler(request: Request, exc: RequestValidationError) -> JSONResponse:
    errors = [
        {"loc": [str(p) for p in err["loc"]], "msg": err["msg"], "type": err["type"]}
        for err in exc.errors()
    ]
    return _validation_response(errors)


@app.get("/health", response_model=HealthResponse, tags=["health"], summary="健康检查")
async def health() -> HealthResponse:
    return HealthResponse(status="ok", service=settings.APP_NAME, version=VERSION)


@app.post(
    "/api/v1/digest/recommend",
    response_model=DigestResponse,
    tags=["digest"],
    summary="枚举酶对并推荐可区分参考型与变异型的双酶切方案",
)
async def digest_recommend(req: DigestRequest) -> DigestResponse:
    validate_request(req)
    return recommend(req)


@app.post(
    "/api/v1/digest/identify",
    response_model=IdentifyResponse,
    tags=["digest"],
    summary="将凝胶观测条带与两种样本的预测条带匹配，判定样本为参考型 / 变异型 / 无法区分",
)
async def digest_identify(req: IdentifyRequest) -> IdentifyResponse:
    validate_identify_request(req)
    return identify(req)


@app.post(
    "/api/v1/digest/panel",
    response_model=PanelResponse,
    tags=["digest"],
    summary="精确求出覆盖全部变异型的最小酶对面板（非贪心近似）",
)
async def digest_panel(req: PanelRequest) -> PanelResponse:
    validate_panel_request(req)
    return build_panel(req)


@app.post(
    "/api/v1/digest/mixture",
    response_model=MixtureResponse,
    tags=["digest"],
    summary="枚举命名候选的成分子集，合并共迁移条带后匹配一次混合样本观测",
)
async def digest_mixture(req: MixtureRequest) -> MixtureResponse:
    validate_mixture_request(req)
    return resolve_mixture(req)


if __name__ == "__main__":
    import uvicorn

    uvicorn.run("app.main:app", host=settings.HOST, port=settings.PORT)
