# Spec: API 服务与 AI 解析功能测试

> **状态**: 待审批
> **优先级**: 高
> **影响范围**: `tests/` 新增测试文件, `allenricher/ai/interpreter.py` 代码修改

---

## 1. 背景与问题

### 1.1 当前测试覆盖

项目已有 7 个测试文件，约 100+ 测试用例，覆盖了核心富集分析算法、CLI 参数解析、数据库解析器等。但以下关键模块**缺少测试**：

| 模块 | 已有测试 | 缺口 |
|------|----------|------|
| **KEGG REST API** (`kegg_fetcher.py`) | ✅ 已完成 | - |
| **FastAPI 服务** (`api/server.py`) | ⚠️ 仅 4 个基础测试 | POST/DELETE 端点、后台任务、文件上传 |
| **AI 解读** (`ai/interpreter.py`) | ⚠️ 仅 MockInterpreter | OpenAI/Claude/Ollama 的 `_call_api()`、提示词构建 |
| **AI 集成** (CLI → AI → 报告) | ❌ 无 | 端到端流程、JSON 输出、HTML 嵌入 |

### 1.2 AI 解读逻辑变更

**当前行为**：各后端提取前 10 条（Ollama 前 5 条）显著富集条目发送给 AI。

**变更后行为**：
- 统一所有后端提取 **Top 20** 显著富集条目
- 如果某数据库显著条目 < 20 条，则按实际数量发送
- 如果某数据库显著条目 = 0 条，直接报告"没有显著富集结果"，**不调用 AI API**

### 1.3 测试策略

- **不使用真实 API Key**：所有外部 HTTP 调用使用 `unittest.mock` 模拟
- **不依赖网络**：OpenAI API、Claude API 均通过 mock 测试
- **不启动服务器**：FastAPI 使用 `TestClient`（基于 httpx）进行内存测试
- **不写入磁盘**：文件输出使用 `tmp_path` pytest fixture

---

## 2. 测试范围

### 2.1 FastAPI 服务测试

**文件**: `tests/test_api_server.py`

| 测试点 | 端点 | 验证内容 |
|--------|------|----------|
| 服务信息 | `GET /` | 返回服务名称和版本 |
| 物种列表 | `GET /api/species` | 返回支持的物种列表 |
| 数据库列表 | `GET /api/databases` | 返回可用数据库列表 |
| 提交分析任务 | `POST /api/analyze` | 请求体验证、job_id 返回、后台任务启动 |
| 文件上传分析 | `POST /api/upload` | 文件解析、基因列表提取、任务创建 |
| 上传大小限制 | `POST /api/upload` | 超过 10MB 返回 413 |
| 查询任务状态 | `GET /api/status/{job_id}` | 正常状态返回、不存在返回 404 |
| 获取分析结果 | `GET /api/results/{job_id}` | JSON/TSV 格式、未完成返回 400 |
| 获取图表 | `GET /api/results/{job_id}/plot` | PDF 文件返回、无效路径返回 404 |
| 获取报告 | `GET /api/results/{job_id}/report` | HTML 文件返回 |
| 删除任务 | `DELETE /api/jobs/{job_id}` | 任务删除、文件清理 |
| 路径遍历防护 | `GET /api/results/{job_id}/plot` | `../../etc/passwd` 被清理 |
| 无效请求体 | `POST /api/analyze` | 缺少字段返回 422 |
| CORS 配置 | 任意端点 | `Access-Control-Allow-Origin` 头 |

### 2.2 AI 解读模块测试

**文件**: `tests/test_ai_interpreter.py`

| 测试点 | 方法 | 验证内容 |
|--------|------|----------|
| OpenAI 调用 | `OpenAIInterpreter._call_api()` | mock openai.ChatCompletion、提示词格式、响应解析 |
| Claude 调用 | `ClaudeInterpreter._call_api()` | mock anthropic.Messages、提示词格式 |
| Ollama 调用 | `OllamaInterpreter._call_api()` | mock HTTP POST、本地 API 格式 |
| DeepSeek 调用 | `DeepSeekInterpreter._call_api()` | mock openai 兼容接口 |
| GLM 调用 | `GLMInterpreter._call_api()` | mock openai 兼容接口 |
| MiniMax 调用 | `MiniMaxInterpreter._call_api()` | mock openai 兼容接口 |
| 依赖缺失降级 | 各 Interpreter | openai/anthropic 未安装时抛出 ImportError |
| **Top 20 条目提取** | `interpret()` | 提取前 20 条显著条目（非 10 条） |
| **不足 20 条时按实际** | `interpret()` | 仅有 5 条结果时发送 5 条 |
| **0 条时跳过 AI** | `interpret()` | 空结果时不调用 API，返回"无显著富集" |
| 条目总结 | `summarize_term()` | 单条目提示词格式 |
| 门面类解读 | `AIInterpreter.interpret_results()` | 调用后端、返回格式 |
| 条目总结开关 | `include_term_summaries` | True 时生成 `{db}_term_summaries` |
| HTML 报告段 | `generate_report_section()` | HTML 格式、免责声明 |
| 工厂函数 | `create_interpreter()` | 后端名称 → 正确实例 |
| 无效后端 | `create_interpreter("invalid")` | 抛出 ValueError |
| 可用后端列表 | `get_available_backends()` | 返回所有已安装依赖的后端 |

### 2.3 AI 集成测试

**文件**: `tests/test_ai_integration.py`

| 测试点 | 验证内容 |
|--------|----------|
| MockInterpreter 端到端 | analyze → AI 解读 → JSON 文件 → HTML 报告 |
| **0 条显著结果** | 直接报告"没有显著富集结果"，不调用 AI |
| **不足 20 条** | 按实际条目数发送给 AI |
| **恰好 20 条** | 发送全部 20 条 |
| AI JSON 输出格式 | `ai_interpretation.json` 结构验证 |
| HTML 报告嵌入 | AI 解读段落正确嵌入 HTML |

---

## 3. 技术方案

### 3.1 AI 解读逻辑变更

**修改文件**: `allenricher/ai/interpreter.py`

**变更点**（所有后端的 `interpret()` 方法统一修改）：

```python
# === 变更前 ===
top_results = df.head(10)

# === 变更后 ===
if len(df) == 0:
    interpretation[db_name] = "No significant enrichment results found for this database."
    continue
top_results = df.head(20)  # 统一 Top 20
```

**涉及的后端**（6 处修改）：
1. `OpenAIInterpreter.interpret()` — `head(10)` → `head(20)` + 空结果检查
2. `ClaudeInterpreter.interpret()` — `head(10)` → `head(20)` + 空结果检查
3. `OllamaInterpreter.interpret()` — `head(5)` → `head(20)` + 空结果检查
4. `DeepSeekInterpreter.interpret()` — `head(10)` → `head(20)` + 空结果检查
5. `GLMInterpreter.interpret()` — `head(10)` → `head(20)` + 空结果检查
6. `MiniMaxInterpreter.interpret()` — `head(10)` → `head(20)` + 空结果检查

**提示词中的文本变更**：
```python
# === 变更前 ===
"Top 10 enriched terms:\n{summary}"

# === 变更后 ===
f"Top {len(top_results)} enriched terms:\n{summary}"
```

**门面类 `AIInterpreter.interpret_results()` 变更**：
```python
# 条目总结也统一为 Top 20（原为 Top 5）
for _, row in df.head(20).iterrows():  # 原 df.head(5)
```

### 3.2 Mock 策略

```python
# OpenAI mock
@patch('allenricher.ai.interpreter.openai')
def test_openai_interpret(mock_openai):
    mock_openai.ChatCompletion.create.return_value = {
        "choices": [{"message": {"content": "AI response"}}]
    }

# FastAPI mock
from fastapi.testclient import TestClient
client = TestClient(app)
response = client.post("/api/analyze", json={...})
```

### 3.3 测试数据

使用 `tmp_path` fixture 创建临时测试数据：

```python
import pandas as pd

# 创建有 25 条结果的 DataFrame（测试 Top 20 截断）
def create_results_df(n=25):
    return pd.DataFrame({
        'Term_Name': [f'term_{i}' for i in range(n)],
        'P_Value': [10**(-i) for i in range(1, n+1)],
        'Gene_Count': [i+1 for i in range(n)],
    })

# 创建空 DataFrame（测试 0 条跳过）
def create_empty_results_df():
    return pd.DataFrame(columns=['Term_Name', 'P_Value', 'Gene_Count'])
```

---

## 4. 实施任务清单

### Task 1: 修改 AI 解读逻辑 — Top 20 + 空结果处理
- 修改 6 个后端的 `interpret()` 方法：`head(10/5)` → `head(20)` + 空结果检查
- 修改提示词文本：`"Top 10"` → `f"Top {len(top_results)}"`
- 修改门面类 `AIInterpreter.interpret_results()`：`head(5)` → `head(20)`
- 运行现有测试确保不破坏

### Task 2: FastAPI 服务测试 (test_api_server.py)
- TestClient 测试所有 10 个端点
- 测试后台任务执行（mock `run_analysis`）
- 测试安全防护（路径遍历、文件大小限制）
- 测试 CORS 配置

### Task 3: AI 解读模块测试 (test_ai_interpreter.py)
- mock 各后端 `_call_api()` 方法
- 测试 Top 20 条目提取逻辑
- 测试不足 20 条时按实际数量
- 测试 0 条时空结果处理（不调用 API）
- 测试门面类 `AIInterpreter`
- 测试工厂函数和后端列表
- 测试依赖缺失降级

### Task 4: AI 集成测试 (test_ai_integration.py)
- MockInterpreter 端到端流程
- 0 条 / 5 条 / 20 条 / 25 条四种场景
- AI JSON 输出验证
- HTML 报告嵌入验证

### Task 5: 运行全部测试并验证
- `pytest tests/ -v --tb=short`
- 确保所有新测试通过
- 确保现有测试不受影响

---

## 5. 预期产出

### 代码变更

| 文件 | 变更 |
|------|------|
| `allenricher/ai/interpreter.py` | 6 个后端 `head()` 统一为 20 + 空结果检查 |

### 新增测试

| 测试文件 | 预计测试数 | 覆盖模块 |
|----------|-----------|----------|
| `test_api_server.py` | ~15 | FastAPI 服务 |
| `test_ai_interpreter.py` | ~20 | AI 解读模块 |
| `test_ai_integration.py` | ~6 | AI 集成流程 |
| **合计** | **~41** | |

加上现有 ~100 个测试，总计 ~141 个测试用例。
