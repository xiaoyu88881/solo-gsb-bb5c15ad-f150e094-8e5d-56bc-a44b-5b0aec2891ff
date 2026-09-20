# 环状 DNA 双酶切鉴别服务

基于 FastAPI 的圆形（环状）DNA 双酶切鉴别服务，提供五个互补接口：

- **方案推荐**：客户端提交质粒长度、参考型与变异型的候选酶切位点、最小可见片段长度和条带分辨容差，服务枚举所有不同酶对，按圆环相邻切点计算两种样本的片段长度，将容差内的片段聚为同一条电泳条带，并推荐能够区分两种样本的酶对。
- **样本鉴别**：客户端提交同一质粒参考型、变异型的预测条带长度和一次凝胶观测长度列表，服务在给定容差内做确定性的一对一最小代价匹配，分别计算两种样本的缺失条带、额外条带与总误差，判定观测样本为参考型、变异型或无法区分。
- **最小酶对面板**：客户端提交质粒长度、可见阈值、分辨容差及酶列表（每种酶给出参考型切点和按变异名映射的变异型切点），服务复用同一套圆环酶切与聚带规则评估每个酶对对每个变异型的可区分条带，并**精确**求出覆盖全部变异型的最小酶对面板（非贪心近似）。
- **混合样本解析**：客户端提交若干**命名候选成分**各自的预测条带、一次混合样本观测条带、匹配容差与成分数上限，服务枚举不超过上限的全部候选子集，先把每个子集的预测条带按容差做单连接聚带（共迁移条带只出现一次、代表长度取均值），再复用鉴别接口的一对一最小代价匹配，按未匹配条带数、总误差、成分数、候选名列表字典序唯一选定最吻合的成分组合。
- **批量盲样分配**：客户端提交若干**命名预测谱型**与若干**命名观测泳道**（盲样）、条带匹配容差与单条分配允许的最大未匹配条带数；服务先对每条"谱型-泳道"边复用鉴别接口的一对一条带匹配计算未匹配数与总误差（超阈值的边不可用），再做全局一对一最优分配——依次最大化已分配泳道数、最小化未匹配总数、最小化总误差，并按泳道名排序后的谱型名序列字典序破同分（未分配排在所有谱型名之后），结果与输入顺序无关。

**核心功能为纯本地数学计算，不依赖数据库、缓存或任何外部服务。**

---

## 单条命令启动

需要 Docker（含 Compose 插件，Docker 20.10+ / Docker Desktop 自带）：

```bash
docker compose up --build -d
```

服务启动后监听 **`http://localhost:8000`**。

停止：

```bash
docker compose down
```

### 不使用 Docker（本地 Python 3.11+）

```bash
pip install -r requirements.txt
uvicorn app.main:app --host 0.0.0.0 --port 8000
```

---

## API 地址与端口

| 用途 | 方法 | 地址 |
| --- | --- | --- |
| 健康检查 | GET | `http://localhost:8000/health` |
| 双酶切方案推荐 | POST | `http://localhost:8000/api/v1/digest/recommend` |
| 凝胶条带样本鉴别 | POST | `http://localhost:8000/api/v1/digest/identify` |
| 最小酶对面板 | POST | `http://localhost:8000/api/v1/digest/panel` |
| 混合样本解析 | POST | `http://localhost:8000/api/v1/digest/mixture` |
| 批量盲样分配 | POST | `http://localhost:8000/api/v1/digest/batch-identify` |
| 交互式 API 文档（Swagger） | GET | `http://localhost:8000/docs` |
| OpenAPI Schema | GET | `http://localhost:8000/openapi.json` |

健康检查示例：

```bash
curl http://localhost:8000/health
# {"status":"ok","service":"circular-dna-digest","version":"1.4.0"}
```

---

## 必要配置

所有配置均有默认值，**开箱即用，无需任何外部服务**。可通过环境变量覆盖：

| 环境变量 | 默认值 | 说明 | Compose 中如何改 |
| --- | --- | --- | --- |
| `APP_HOST` | `0.0.0.0` | 服务监听地址 | 编辑 `docker-compose.yml` 的 `environment` |
| `APP_PORT` | `8000` | 容器内监听端口 | 编辑 `docker-compose.yml` 的 `environment` |
| `HOST_PORT` | `8000` | 宿主机映射端口 | `HOST_PORT=9000 docker compose up --build -d` |
| `MAX_ENZYMES` | `64` | 单次请求允许的最大酶数量（酶对数为 C(n,2)） | 编辑 `environment` |
| `MAX_BANDS` | `512` | identify / mixture 接口每个条带列表允许的最大条带数 | 编辑 `environment` |
| `MAX_VARIANTS` | `16` | panel 接口允许的最大变异型数量（最小面板为精确集合覆盖，规模随变异数指数增长） | 编辑 `environment` |
| `MAX_MIXTURE_CANDIDATES` | `12` | mixture 接口允许的最大候选成分数（需枚举至多该规模的全部子集） | 编辑 `environment` |
| `MAX_BATCH_ITEMS` | `32` | batch-identify 接口允许的最大谱型数 / 泳道数（全局一对一分配为 O(k³)） | 编辑 `environment` |

---

## 请求与响应

### 双酶切方案推荐（POST `/api/v1/digest/recommend`）

#### 请求体（JSON）

| 字段 | 类型 | 约束 |
| --- | --- | --- |
| `plasmid_length` | int | > 0，环状质粒总长度（bp），切点合法范围为 `0 ~ plasmid_length-1` |
| `min_visible_fragment` | int | ≥ 0 且 < 质粒长度；短于此长度的片段不参与条带聚类 |
| `band_resolution_tolerance` | int | ≥ 0，条带分辨容差（bp） |
| `enzymes` | array | 至少 2 种酶，酶名称不可重复 |
| `enzymes[].name` | string | 酶名称，如 `EcoRI` |
| `enzymes[].reference_sites` | int[] | 参考型候选切点（0-based，bp） |
| `enzymes[].variant_sites` | int[] | 变异型候选切点（0-based，bp） |

请求示例见 [`examples/request.json`](examples/request.json)：

```bash
curl -X POST http://localhost:8000/api/v1/digest/recommend \
  -H 'Content-Type: application/json' \
  -d @examples/request.json
```

#### 算法说明（推荐）

1. **枚举酶对**：对所有不同酶的两两组合（酶名按字典序规范化）。
2. **合并切点**：每种样本合并两酶切点；两种酶切在同一位置时合并为一个切点并记录来源酶（不视为重复位点）。
3. **圆环片段**：切点升序排列为 p0 < p1 < … < p_{n-1}，相邻切点间形成片段；末段 `(p_{n-1}, p0)` 跨过坐标原点，长度为 `p0 + L - p_{n-1}`，响应中以 `crosses_origin: true` 标记。
4. **可见性过滤**：长度 < `min_visible_fragment` 的片段不参与成带（计数为 `hidden_fragment_count`）。
5. **条带聚类**：两样本可见片段合并并按长度排序，相邻片段长度差 ≤ `band_resolution_tolerance` 归入同一条带（单连接聚类），超过则断开成新带。
6. **可区分条带**：条带内参考型与变异型的片段数量不同（`reference_only` / `variant_only` / `count_mismatch`）。
7. **推荐排序**：可区分条带数降序 → 总条带数升序 → 酶对名称字典序升序。

#### 响应要点（推荐）

- `recommendations`：推荐酶对，包含
  - `samples.reference` / `samples.variant`：合并后的**切点**（含来源酶）、全部**片段**及跨原点标记、可见/不可见片段计数；
  - `bands`：聚类后的**条带**（代表长度、长度范围、两样本各自片段数、成员片段、是否可区分）；
  - `differences`：**差异依据**（中文说明 + 机器可读类型）。
- `rejected_pairs`：**任一样本有效切点不足 2 个**的酶对（含原因与各自切点数量）。
- `non_discriminating_pairs`：切点足够但两样本条带完全一致、无法鉴别的酶对。
- `total_enzyme_pairs` / `recommended_count` / `rejected_pair_count` / `non_discriminating_pair_count`：汇总计数。

### 凝胶条带样本鉴别（POST `/api/v1/digest/identify`）

#### 请求体（JSON）

| 字段 | 类型 | 约束 |
| --- | --- | --- |
| `reference_bands` | number[] | 参考型预测条带长度（bp），每条 ≥ 0，可为空列表 |
| `variant_bands` | number[] | 变异型预测条带长度（bp），每条 ≥ 0，可为空列表 |
| `observed_bands` | number[] | 一次凝胶观测到的条带长度（bp），每条 ≥ 0，**至少 1 条** |
| `band_match_tolerance` | number | ≥ 0，条带匹配容差（bp） |

请求示例见 [`examples/identify-request.json`](examples/identify-request.json)：

```bash
curl -X POST http://localhost:8000/api/v1/digest/identify \
  -H 'Content-Type: application/json' \
  -d @examples/identify-request.json
```

#### 匹配与判定规则

1. **一对一最小代价匹配**：对每种样本，将其预测条带与观测条带分别按长度升序排列后做一对一配对；仅当 `|预测长度 − 观测长度| ≤ band_match_tolerance` 时允许配对。目标为先**最大化匹配对数**，再**最小化总误差**（各配对长度差之和）。一维绝对差代价满足 Monge 性质，最优解必可取不交叉匹配，服务用序列比对式动态规划在 O(n·m) 内求得；并列时按固定优先级选择（等效于较短条带优先配对），结果完全确定且与列表提交顺序无关。
2. **缺失 / 额外条带**：未匹配到观测条带的预测条带为**缺失条带**（`missing_bands`），未匹配到任何预测条带的观测条带为**额外条带**（`extra_bands`）。
3. **判定**：比较两种样本的未匹配条带数（缺失 + 额外），少者判为对应样本（`reference` / `variant`）；并列时总误差小者胜出；两者完全一致时判定为 `ambiguous`（无法区分）。

#### 响应要点（鉴别）

- `samples.reference` / `samples.variant`：每种样本的
  - `matched_count` / `missing_count` / `extra_count` / `total_error`：匹配对数、缺失数、额外数与总误差（bp，保留两位小数）；
  - `matches`：**逐条匹配依据**（预测长度、观测长度、长度差及中文说明），按预测长度升序；
  - `missing_bands` / `extra_bands`：缺失、额外条带长度列表（升序）。
- `verdict`：`reference` / `variant` / `ambiguous`；`verdict_reason` 为中文判定依据，可完全由响应中的计数与总误差复算。

### 最小酶对面板（POST `/api/v1/digest/panel`）

#### 请求体（JSON）

| 字段 | 类型 | 约束 |
| --- | --- | --- |
| `plasmid_length` | int | > 0，环状质粒总长度（bp），切点合法范围为 `0 ~ plasmid_length-1` |
| `min_visible_fragment` | int | ≥ 0 且 < 质粒长度；短于此长度的片段不参与条带聚类 |
| `band_resolution_tolerance` | int | ≥ 0，条带分辨容差（bp） |
| `enzymes` | array | 至少 2 种酶，酶名称不可重复；所有酶的变异名集合必须完全一致 |
| `enzymes[].name` | string | 酶名称，如 `EcoRI` |
| `enzymes[].reference_sites` | int[] | 参考型候选切点（0-based，bp） |
| `enzymes[].variant_sites` | object | 按变异名映射的候选切点：`{"变异名": [切点...]}`，至少定义 1 个变异名 |

请求示例见 [`examples/panel-request.json`](examples/panel-request.json)：

```bash
curl -X POST http://localhost:8000/api/v1/digest/panel \
  -H 'Content-Type: application/json' \
  -d @examples/panel-request.json
```

#### 算法说明（面板）

1. **逐对逐变异评估**：对每一对不同酶、每一个变异型，复用推荐接口的圆环酶切与聚带规则计算参考型 vs 该变异型的条带；任一样本有效切点不足 2 个时该酶对对此变异型不可用。存在至少一条可区分条带时，称该酶对**覆盖**此变异型。
2. **精确最小面板**：把每个酶对的覆盖变异型集合表示为位掩码，按位或 BFS 逐层扩展（第 k 层为恰好 k 个酶对可达的覆盖掩码），首次到达全掩码的层即为最小面板规模——这是最小集合覆盖的**精确解，不使用贪心近似**。同一掩码只保留最优代表酶对（可区分条带总数高者优先，并列取 `pair_key` 小者），同层同掩码只保留最优取值，保证结果完备且确定。
3. **同规模选定**：最小规模相同的方案依次按**可区分条带总数降序**、**排序后的 pair_key 列表字典序**选定。

#### 响应要点（面板）

- `selected_pairs`：最小面板中的酶对（按 `pair_key` 字典序），每个酶对包含
  - `covered_variants`：该酶对覆盖的变异名；
  - `coverages`：**逐变异的差异条带依据**（可区分条带数 + 与推荐接口同结构的 `differences` 明细）；
  - `distinguishable_band_count`：该酶对在全部变异型上的可区分条带总数。
- `covered_variants` / `uncovered_variants`：面板覆盖的变异名 / 无任何可用酶对的变异名。
- `total_distinguishable_bands`：面板全部酶对的可区分条带总数。
- **任一变异型无可用酶对时返回 HTTP 200**：该变异列入 `uncovered_variants`，且 `selected_pairs` 为空——**不返回部分面板**。

### 混合样本解析（POST `/api/v1/digest/mixture`）

#### 请求体（JSON）

| 字段 | 类型 | 约束 |
| --- | --- | --- |
| `components` | array | 命名候选成分列表，**至少 1 个**候选，候选名不可重复，总数不超过 `MAX_MIXTURE_CANDIDATES`（默认 12） |
| `components[].name` | string | 候选成分名（非空白），如 `strain-A` |
| `components[].predicted_bands` | number[] | 该候选成分单独的预测条带长度（bp），每条 ≥ 0，**至少 1 条** |
| `observed_bands` | number[] | 一次凝胶观测到的混合样本条带长度（bp），每条 ≥ 0，**至少 1 条** |
| `band_match_tolerance` | number | ≥ 0，条带容差（bp）：同时用于共迁移聚带与预测/观测配对 |
| `max_components` | int | ≥ 1 且 ≤ 候选总数；枚举成分数 1..该值的全部非空候选子集 |

请求示例见 [`examples/mixture-request.json`](examples/mixture-request.json)：

```bash
curl -X POST http://localhost:8000/api/v1/digest/mixture \
  -H 'Content-Type: application/json' \
  -d @examples/mixture-request.json
```

#### 算法说明（混合）

1. **枚举候选子集**：对成分数 1..`max_components` 的全部非空候选子集逐一评估（候选名按字典序规范化，结果与提交顺序无关）。
2. **共迁移聚带（单连接）**：合并子集内各候选的预测条带（多重集），按长度升序扫描，相邻条带长度差 ≤ `band_match_tolerance` 即归入同一簇（链式可传递的单连接聚类）。**共迁移条带只出现一次**，无论它来自多个候选还是同一候选内的重复预测；每簇的**代表长度取簇内全部条带长度的均值**（保留两位小数）。
3. **一对一最小代价匹配**：把合并预测条带与观测条带复用鉴别接口的同一套确定性匹配（先最大化匹配对数，再最小化总误差），得到逐条配对、缺失条带与额外条带。
4. **唯一选定**：全部评估子集依次按 ① **未匹配条带总数**（缺失 + 额外）升序、② **总误差**升序、③ **成分数**升序、④ **排序后的候选名列表字典序**升序比较，次序④保证即使指标完全相同也有唯一结果。

#### 响应示例（混合）

```json
{
  "band_match_tolerance": 10.0,
  "observed_band_count": 4,
  "max_components": 2,
  "selected_components": ["strain-A", "strain-B"],
  "merged_predicted_bands": [230.0, 319.0, 455.0, 700.0],
  "match": {
    "predicted_band_count": 4,
    "matched_count": 4,
    "missing_count": 0,
    "extra_count": 0,
    "total_error": 3.0,
    "matches": [
      {"predicted_length": 230.0, "observed_length": 231.0, "error": 1.0, "description": "..."},
      {"predicted_length": 319.0, "observed_length": 319.0, "error": 0.0, "description": "..."},
      {"predicted_length": 455.0, "observed_length": 455.0, "error": 0.0, "description": "..."},
      {"predicted_length": 700.0, "observed_length": 702.0, "error": 2.0, "description": "..."}
    ],
    "missing_bands": [],
    "extra_bands": []
  },
  "evaluated_subset_count": 3,
  "selection_reason": "在 3 个成分数不超过 2 的候选子集中，选定 ['strain-A', 'strain-B']（2 个成分）……"
}
```

- `selected_components`：选中的候选成分名（字典序）；`merged_predicted_bands`：选中子集聚带后的代表长度（升序）。
- `match`：与鉴别接口**完全同结构**——逐条匹配依据（按预测长度升序）、`missing_bands` / `extra_bands` 及各项计数与总误差。
- `evaluated_subset_count`：实际评估的子集数，等于 `Σ C(候选数, k)`（k = 1..`max_components`）。
- `selection_reason`：中文选定依据，说明四级比较次序如何决定唯一结果。

### 批量盲样分配（POST `/api/v1/digest/batch-identify`）

一次把多条**盲样泳道**分配到多个**命名预测谱型**，每个谱型至多使用一次、每条泳道至多分配一个谱型（全局一对一）。

#### 请求体（JSON）

| 字段 | 类型 | 约束 |
| --- | --- | --- |
| `profiles` | array | 命名预测谱型列表，**至少 1 个**，谱型名不可重复，总数不超过 `MAX_BATCH_ITEMS`（默认 32） |
| `profiles[].name` | string | 谱型名（非空白），如 `strain-A` |
| `profiles[].predicted_bands` | number[] | 该谱型的预测条带长度（bp），每条 ≥ 0，**至少 1 条** |
| `lanes` | array | 命名观测泳道列表，**至少 1 条**，泳道名不可重复，总数不超过 `MAX_BATCH_ITEMS`（默认 32）；泳道名与谱型名分属不同命名空间，允许同名 |
| `lanes[].name` | string | 泳道名（非空白），如 `lane-1` |
| `lanes[].observed_bands` | number[] | 该泳道观测到的条带长度（bp），每条 ≥ 0，**至少 1 条** |
| `band_match_tolerance` | number | ≥ 0，条带匹配容差（bp）：长度差不超过该值的预测/观测条带才允许配对 |
| `max_unmatched_bands` | int | ≥ 0，**单条分配允许的最大未匹配条带数**（缺失 + 额外）；某谱型-泳道边的未匹配数超过该值时该边不可用 |

请求示例见 [`examples/batch-identify-request.json`](examples/batch-identify-request.json)：

```bash
curl -X POST http://localhost:8000/api/v1/digest/batch-identify \
  -H 'Content-Type: application/json' \
  -d @examples/batch-identify-request.json
```

#### 算法说明（批量分配）

1. **逐边一对一条带匹配**：对全部 `谱型数 × 泳道数` 条边，复用鉴别接口的同一套确定性匹配（条带各自按长度升序、`|预测 − 观测| ≤ band_match_tolerance` 才配对，先最大化匹配对数再最小化总误差），得到每条边的匹配对数、缺失 / 额外条带、未匹配条带数与总误差。**未匹配条带数（缺失 + 额外）超过 `max_unmatched_bands` 的边标记为不可用**，不参与后续分配（边界含等号：恰好等于阈值可用）。
2. **全局一对一最优分配**：补虚拟泳道（吸收未使用谱型）与虚拟谱型（吸收未分配泳道）后，在可用边上做赋权二分图最小权完美匹配（O(k³) 匈牙利算法）。四个依次比较的目标被编码为一个**混合进制整数**（低序目标的取值范围严格小于高序目标的进制，低位相加不会进位），保证整数大小次序与下列字典序完全一致：
   1. **已分配泳道数最多**（首要目标；宁可留谱型未使用，也要尽量多分配泳道）；
   2. **未匹配总数最少** = 全部已分配边的未匹配数之和 + 未分配泳道的全部观测条数（未分配泳道的观测条带全部计为未匹配）；
   3. **总误差最小** = 全部已分配边匹配对长度差之和；
   4. 并列时，泳道先按泳道名字典序排序，再以其**对应谱型名序列的字典序**破同分，**未分配项排在所有谱型名之后**。
3. **确定性**：全部决策只依赖排序后的泳道名 / 谱型名与边匹配结果，**响应不受谱型、泳道及条带在请求中的提交顺序影响**。

#### 响应要点（批量分配）

- `assignments`：**逐项匹配依据**，按泳道名字典序排序，每项包含
  - `lane` / `assigned_profile`（未分配时为 `null`）/ `usable`；
  - `unmatched_count` / `total_error`（已分配为所用边的值；未分配时分别为该泳道全部观测条数与 0）；
  - `match`：已分配泳道的逐条配对依据（与鉴别接口**完全同结构**，含 `matches` / `missing_bands` / `extra_bands`），未分配泳道为 `null`；
  - `description`：中文分配依据（说明无可用边、还是可用谱型在全局分配中被其他泳道占用）。
- `unassigned_lanes`：**未分配泳道**名（字典序）；`unused_profiles`：**未使用谱型**名（字典序）。
- `edges`：**全部谱型-泳道边**的逐项依据（按泳道名、谱型名字典序排序）：`usable`、`unmatched_count` / `missing_count` / `extra_count` / `total_error`，可用边附完整 `match` 明细、不可用边 `match` 为 `null`，并给出中文 `reason`。
- `assigned_lane_count` / `total_unmatched_count` / `total_error`：全局汇总，恒等于 `assignments` 逐项之和，可由响应直接复算；`assignment_reason` 说明阈值与四级比较次序如何决定唯一结果。

### 拒绝条件（HTTP 422）

以下输入问题整体拒绝，返回统一结构 `{"code":"validation_failed","errors":[...]}`，一次性报告全部错误：

推荐接口（`/api/v1/digest/recommend`）：

- 切点越界（`< 0` 或 `≥ plasmid_length`）；
- 同一样本、同一酶的候选切点重复；
- 酶名称重复；酶少于 2 种；
- 质粒长度 ≤ 0、最小可见片段 ≥ 质粒长度、容差为负等字段约束错误。

鉴别接口（`/api/v1/digest/identify`）：

- 观测条带列表为空（`observed_bands` 至少 1 条）；
- 任一条带长度（参考 / 变异 / 观测）为负数，或匹配容差为负、缺失、非数值；
- 任一条带列表超过 `MAX_BANDS` 上限。

面板接口（`/api/v1/digest/panel`）：

- 参考型或任一变异型的切点越界（`< 0` 或 `≥ plasmid_length`）、同一酶同一样本位点重复；
- 酶名称重复、酶少于 2 种、酶数量超过 `MAX_ENZYMES` 上限；
- 各酶的变异名集合不一致、变异名为空白、未定义任何变异型、变异型数量超过 `MAX_VARIANTS` 上限；
- 质粒长度 ≤ 0、最小可见片段 ≥ 质粒长度、容差为负等字段约束错误。

混合接口（`/api/v1/digest/mixture`）：

- `max_components` 非法（缺失、非整数、< 1）或**超过候选成分总数**；
- 候选成分名**重复**或为空白；候选列表为空；
- 任一候选的 `predicted_bands` 为空（**空候选**），或观测列表为空；
- 候选成分数超过 `MAX_MIXTURE_CANDIDATES`（默认 12）上限；
- 任一条带长度为负数、容差为负 / 缺失 / 非数值等字段约束错误；
- 任一条带列表超过 `MAX_BANDS` 上限。

批量分配接口（`/api/v1/digest/batch-identify`）：

- 谱型名或泳道名**重复**或为空白；谱型 / 泳道列表为空；
- 谱型数或泳道数超过 `MAX_BATCH_ITEMS`（默认 32）上限；
- 任一谱型的 `predicted_bands` 或任一泳道的 `observed_bands` 为空；
- `band_match_tolerance` 或 `max_unmatched_bands` 为负 / 缺失 / 非数值；
- 任一条带长度为负数等字段约束错误；
- 任一条带列表超过 `MAX_BANDS` 上限。

注意：推荐接口中切点不足两个属于**酶对级别**的拒绝（记录在 `rejected_pairs`，HTTP 200），而非请求级别错误；面板接口中某变异型无可用酶对时同样返回 HTTP 200（列入 `uncovered_variants`，不返回部分面板）；批量分配接口中某条边未匹配数超阈值属于**边级别**的不可用（记录在 `edges[].usable=false` 与未分配结果中，HTTP 200），所有谱型-泳道边都不可用时返回的是"全部泳道未分配"的成功响应，而非 422。

---

## 项目结构

```
.
├── app/
│   ├── main.py          # FastAPI 入口：健康检查、推荐/鉴别/面板/混合/批量分配路由、统一 422 处理
│   ├── models.py        # 请求/响应 Pydantic 模型
│   ├── validation.py    # 越界、重复位点、变异名一致性、条带/候选规模等语义校验
│   ├── digest.py        # 圆环片段、条带聚类、酶对枚举与排序（纯函数）
│   ├── identify.py      # 条带一对一最小代价匹配与样本判定（纯函数）
│   ├── panel.py         # 逐变异评估与精确最小酶对面板（纯函数）
│   ├── mixture.py       # 混合样本共迁移聚带、子集枚举与最优子集选定（纯函数）
│   ├── batch.py         # 批量盲样逐边匹配筛边与全局一对一最优分配（纯函数）
│   └── config.py        # 环境变量配置（均有默认值）
├── tests/               # pytest 单元 + API 端到端测试
├── examples/request.json
├── examples/identify-request.json
├── examples/panel-request.json
├── examples/mixture-request.json
├── examples/batch-identify-request.json
├── Dockerfile
├── docker-compose.yml
└── requirements.txt
```

## 运行测试

```bash
pip install -r requirements.txt pytest httpx
python3 -m pytest -q
```
