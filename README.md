# MS Precompute API

用于质谱预计算的 API 服务：客户端提交含元素计数的分子式（可带同位素标记），服务按各元素天然同位素丰度将分布表示为质量-强度多项式并逐元素卷积，合并质量差小于 0.0001 Da 的峰，归一化基峰强度（=100）后返回前 N 个峰及分子单同位素质量。

技术栈：Python 3.12 · FastAPI · SQLite（标准库 `sqlite3`）。

## 启动（单条命令）

```bash
docker compose up --build
```

依赖安装（镜像构建时）与数据库初始化（建表 + 同位素参考数据播种）均随启动自动完成，无需任何手工步骤。

## 访问地址

| 地址 | 说明 |
| --- | --- |
| http://localhost:8000 | API 服务（端口 `8000`） |
| http://localhost:8000/docs | Swagger 交互式文档 |
| http://localhost:8000/health | 健康检查 |

## 配置

通过环境变量配置（均有默认值，非必需）：

| 变量 | 默认值 | 说明 |
| --- | --- | --- |
| `DATABASE_PATH` | `/data/msprecompute.db`（容器）/ `./msprecompute.db`（本地） | SQLite 数据库文件路径 |

Compose 中数据库文件持久化在命名卷 `msdata`；端口映射可在 `docker-compose.yml` 的 `ports` 中调整。

## API 用法

### 计算同位素分布

`POST /api/v1/isotope-distribution`

```json
{
  "formula": "C6H12O6",
  "labels": [{"isotope": "13C", "count": 2}],
  "top_n": 5
}
```

- `formula`（必填）：分子式，如 `C6H12O6`；标记也可内嵌写作 `C6H12O6[13C]2`。标记表示该元素中有多少个原子被标记（`C6[13C]2` = 共 6 个碳，其中 2 个为 ¹³C），标记数不能超过该元素原子总数。
- `labels`（可选）：结构化标记列表，与内嵌标记合并校验；同位素须为参考表中该元素的已知同位素（如 `13C`、`15N`、`2H`、`18O`）。
- `top_n`（可选，默认 10）：返回强度最高的前 N 个峰。

响应：

```json
{
  "formula": "C6H12O6",
  "monoisotopic_mass": 182.070098,
  "peaks": [
    {"mass": 182.070098, "intensity": 100.0},
    {"mass": 183.073453, "intensity": 4.3263},
    {"mass": 184.074343, "intensity": 1.233},
    {"mass": 183.074315, "intensity": 0.2286},
    {"mass": 183.076375, "intensity": 0.138}
  ],
  "total_peaks": 46,
  "cached": false
}
```

计算结果按规范化组成缓存进 SQLite（`calculations` 表），重复请求直接命中预计算结果（`cached: true`）。

### 查询支持的元素与同位素

`GET /api/v1/elements` — 返回元素、质量数、同位素质量与天然丰度（启动时由内置数据播种至 SQLite 的 `isotopes` 表）。

### 错误报告

非法分子式或标记返回 HTTP 422，并指出出错位置（`position` 为 `formula` 中的字符下标；若错误在请求体 `labels` 中，则为列表下标）：

```json
{
  "error": "invalid_formula",
  "message": "unknown element symbol 'Xy'",
  "position": 2
}
```

## 计算说明

- 每种元素的同位素分布视为质量-强度多项式 Σ(丰度·x^质量），按原子数做二项式幂运算后逐元素卷积；每步合并质量差 < 0.0001 Da 的峰（强度加权平均质量、强度求和），并裁剪相对强度 < 1e-9 的峰以控制规模。
- 单同位素质量 = 各元素未标记原子 × 最轻天然同位素质量 + 标记原子 × 标记同位素质量。
- 基峰（最强峰）强度归一化为 100。
- 分子式原子总数上限为 10000。

## 本地开发与测试

核心算法（解析器、卷积引擎）不依赖 FastAPI，可直接用标准库自测：

```bash
python3 tests/test_core.py
```

本地运行服务（需 Python 3.11+）：

```bash
pip install -r requirements.txt
uvicorn app.main:app --host 0.0.0.0 --port 8000
```
