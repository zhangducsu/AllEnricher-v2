# AllEnricher v2.0 — Claude Code Context

> 本文件为 Claude Code / Trae Solo 提供项目上下文，用于本地继续开发和测试。

---

## 1. 项目概述

**AllEnricher** 是一个综合性的基因集功能富集分析工具，支持多种数据库、算法和物种，并提供 AI 驱动的结果解读。

- **仓库**: https://github.com/zhangducsu/AllEnricher-v2
- **版本**: 2.0.0 (Beta)
- **许可证**: MIT
- **Python**: >=3.8
- **代码量**: ~12,000 行 Python 代码（22 个源文件 + 12 个测试文件）
- **测试**: 679 个测试全部通过

---

## 2. 技术选型

| 层面 | 技术 | 说明 |
|------|------|------|
| 语言 | Python 3.8+ | 核心分析引擎 |
| 统计计算 | scipy, numpy, statsmodels | Fisher精确检验、超几何检验、多重检验校正 |
| 数据处理 | pandas | DataFrame 操作、TSV/Excel 读写 |
| CLI | argparse | 命令行接口 |
| Web API | FastAPI + uvicorn | REST API 服务 |
| 可视化 | matplotlib, seaborn, cutecharts | Python 原生实现，支持 PNG/PDF 双输出 |
| 报告 | Jinja2 风格 HTML 模板 | 交互式 HTML 报告 |
| AI | OpenAI SDK >=1.0, Anthropic SDK | GPT-4 / Claude 结果解读 |
| 配置 | PyYAML + dataclass | YAML 配置文件 + Python 类型安全 |
| 测试 | pytest | 单元测试 + 集成测试 |
| 构建 | setuptools + pyproject.toml | PEP 518 标准构建 |

### 核心依赖

```
numpy>=1.20.0, pandas>=1.3.0, scipy>=1.7.0, statsmodels>=0.12.0,
pyyaml>=5.4.0, requests>=2.25.0, tqdm>=4.60.0, openpyxl>=3.0.0
```

### 可选依赖

```
[api]      fastapi>=0.68.0, uvicorn>=0.15.0, python-multipart>=0.0.5
[ai]       openai>=1.0.0, anthropic>=0.18.0
[vis]      matplotlib>=3.4.0, seaborn>=0.11.0
[dev]      pytest>=6.0.0, pytest-cov>=2.12.0, black>=21.0
```

---

## 3. 项目架构

```
AllEnricher-v2/
├── pyproject.toml              # 项目配置、依赖、工具配置
├── README.md                   # 项目文档
├── .gitignore
├── allenricher/
│   ├── __init__.py             # 包入口，导出核心类
│   ├── cli.py                  # CLI 命令行接口 (argparse)
│   ├── core/
│   │   ├── config.py           # Config dataclass + YAML 模板 + 校验
│   │   └── enrichment.py       # 核心分析引擎（算法、两遍扫描、校正、过滤）
│   ├── database/
│   │   ├── manager.py          # DatabaseManager + 7个数据库子类
│   │   └── species_lookup.py   # SpeciesLookup 物种检索（内置16种 + KEGG API）
│   ├── visualization/
│   │   ├── __init__.py         # 可视化模块入口
│   │   ├── barplot.py          # 柱状图（matplotlib实现）
│   │   ├── bubble.py           # 气泡图（matplotlib实现）
│   │   ├── color_config.py     # 颜色配置（36组风格+暗色模式）
│   │   ├── common_plots.py     # 通用图表（网络图/火山图等）
│   │   ├── gsea_plots.py       # GSEA专属图表（富集曲线/NES条形图）
│   │   ├── gsva_plots.py       # GSVA/ssGSEA专属图表（热图/组间比较）
│   │   ├── plot_config.py      # 图表配置参数
│   │   ├── plot_theme.py       # 主题系统（10种预置主题+自定义）
│   │   └── plotter.py          # Plotter 统一调度类（Python原生实现）
│   ├── report/
│   │   └── generator.py        # ReportGenerator HTML 报告生成
│   ├── ai/
│   │   └── interpreter.py      # AIInterpreter AI 解读（OpenAI/Claude）
│   └── api/
│       └── server.py           # FastAPI REST API 服务器
└── tests/
    ├── test_enrichment.py           # 基础单元测试（14个）
    ├── test_enrichment_extended.py  # 扩展单元测试（50个）
    ├── test_database.py             # 数据库模块测试（8个）
    ├── test_cli.py                  # CLI测试（15个）
    ├── test_phase5.py               # API/AI/报告模块测试（19个）
    ├── test_e2e_gsea.py             # GSEA端到端测试（7个）
    ├── test_e2e_ssgsea.py           # ssGSEA端到端测试（8个）
    ├── test_e2e_gsva.py             # GSVA端到端测试（14个）
    ├── test_e2e_visualization.py    # 可视化集成测试（26个）
    ├── test_gmt_generation_e2e.py   # GMT生成测试（8个）
    └── fixtures/
        ├── test_genes.txt      # 测试基因列表（200个，来自 v1 example.glist）
        └── database/hsa/       # 测试数据库（截取自 v1）
            ├── hsa.kegg2gene.tab.gz   # 20个 KEGG 通路
            ├── hsa.kegg2disc.gz
            ├── hsa.GO2gene.tab.gz     # 30个 GO term
            └── hsa.GO2disc.gz
```

---

## 4. 核心模块详解

### 4.1 分析引擎 (`core/enrichment.py`)

**关键类层次**:
```
EnrichmentMethodBase (ABC)
├── FisherExactTest          # Fisher 精确检验
├── HypergeometricTest       # 超几何检验（复用 Fisher 计算统计量）
├── GSEA                     # 基因集富集分析（含置换检验）
└── SSGSEA                   # 单样本 GSEA
```

**EnrichmentAnalyzer** — 主分析器:
- `method` 属性：延迟初始化（property），首次访问时根据 `config.method` 创建
- `run_analysis(gene_set, background_set, database_data)` — 完整分析流程
- `analyze_database(gene_set, background_set, term_data, database)` — **两遍扫描逻辑**（v1.0 兼容）
- `adjust_pvalues(results, method)` — 多重检验校正（BH/BY/Bonferroni/Holm/NONE）
- `filter_results(results)` — 按 pvalue/qvalue/min_genes/max_genes 过滤
- `save_results(output_dir)` — 保存为 TSV 文件

**两遍扫描语义（与 v1.0 一致）**:
1. **第一遍**: 遍历所有条目，收集至少命中一个条目的输入基因 → `gene_total`
2. **第二遍**: 使用 `gene_total` 和 `background_total`（注释文件中所有基因数）计算 p 值

**EnrichmentResult** — 结果 dataclass:
- 字段: `term_id, term_name, database, pvalue, adjusted_pvalue, gene_count, background_count, expected_count, rich_factor, gene_list, gene_ratio, background_ratio, term_url, nes, es, fdr, leading_edge`
- `to_dict()` 返回大写下划线列名（如 `Term_ID`, `Adjusted_P_Value`, `Gene_Ratio`）

### 4.2 配置系统 (`core/config.py`)

- `Config` dataclass：所有配置项的类型安全容器
- `DEFAULT_CONFIG_YAML`：YAML 模板字符串（与 Config 默认值保持一致）
- `validate()` 方法：检查 input_file/background_file 存在性
- `SPECIES_CONFIGS`：16 个预配置物种的 KEGG 代码、taxid、名称
- `EnrichmentMethod` 枚举：`fisher, hypergeometric, gsea, ssgsea`
- `CorrectionMethod` 枚举：`BH, BY, bonferroni, holm, NONE`

### 4.3 数据库管理 (`database/manager.py`)

**DatabaseManager**:
- `load_database(name)` — 加载单个数据库（GODatabase/KEGGDatabase/ReactomeDatabase 等）
- `get_all_term_data()` — 返回 `{db_name: {term_id: {"name": str, "genes": List[str]}}}`
- `get_background_genes()` — 返回所有数据库基因的并集
- `database_dir` 参数：指向包含 `{species}.XX2gene.tab.gz` 文件的目录

**数据库文件格式**（与 v1 兼容）:
- gzip 压缩的 TSV 文件
- 第一行表头：`Gene\tTermID1\tTermID2\t...`
- 后续行：`GeneSymbol\t0\t1\t...`（0/1 表示基因是否属于该条目）

**支持的数据库**: GO, KEGG, Reactome, DO, DisGeNET, WikiPathways, MSigDB, Custom

### 4.4 物种检索 (`database/species_lookup.py`)

- `SpeciesLookup` 类：内置 16 个物种 + KEGG API 在线查询
- 支持按拉丁名、KEGG 代码、taxid 搜索
- `_load_builtin_species()` 和 `BUILTIN_TAXID_MAP` 中的 taxid 必须保持一致

### 4.5 可视化 (`visualization/` 模块)

- **Python 原生实现**：使用 matplotlib + seaborn + cutecharts，完全替代 R/ggplot2
- **模块化设计**：9 个 Python 源文件按职责拆分，各司其职
  - `plotter.py` — 统一调度入口，协调各图表生成
  - `barplot.py` / `bubble.py` — 基础富集分析图表
  - `gsea_plots.py` / `gsva_plots.py` — GSEA/GSVA/ssGSEA 专属图表
  - `common_plots.py` — 通用图表（网络图、火山图、UpSet图等）
  - `plot_theme.py` — 主题系统（10种预置主题 + 自定义主题）
  - `color_config.py` — 36组颜色风格配置 + 暗色模式支持
  - `plot_config.py` — 图表参数配置
- **支持图表**：barplot, bubble, dotplot, enrichment_map, cnet, heatmap, upset, 火山图等 12+ 种
- **输出格式**：PNG + PDF 双输出，支持 300DPI 发表级分辨率
- **`plot_all()` 方法**：为每个数据库生成所有图表类型

### 4.6 报告生成 (`report/generator.py`)

- `generate(results, output_file, gene_list, ai_interpretation)` — 生成 HTML 报告
- 导航栏根据实际数据库动态生成
- AI Interpretation 部分条件显示（仅有内容时才显示链接）

### 4.7 AI 解读 (`ai/interpreter.py`)

- 支持 OpenAI (>=1.0.0 新版 SDK) 和 Anthropic Claude
- `interpret(results, top_n)` — 对 top N 显著条目生成生物学解读

### 4.8 REST API (`api/server.py`)

- FastAPI 应用，同步后台任务（`run_analysis` 为普通函数，非 async）
- 端点: `/api/analyze`, `/api/upload`, `/api/status/{job_id}`, `/api/results/{job_id}`, `/api/results/{job_id}/plot`, `/api/results/{job_id}/report`, `/api/species`, `/api/jobs/{job_id}`
- 安全: 路径遍历防护、文件上传大小限制(10MB)、pickle 反序列化验证
- Pydantic v2: 使用 `model_dump()` 而非 `dict()`

---

## 5. 关键设计决策

| 决策 | 原因 |
|------|------|
| `method` 属性延迟初始化 | 允许构造 Analyzer 后再修改 config.method |
| 两遍扫描语义 | 与 v1.0 Perl 版本结果完全一致 |
| `subprocess.run` 替代 `os.system` | 安全（无 shell 注入）、超时控制、异常捕获 |
| 同步后台任务（非 async） | 分析操作为 CPU 密集型，async 无优势且会阻塞事件循环 |
| `calculate_enrichment` 保留为简化 API | 虽然被 `analyze_database` 绕过，但作为公共 API 仍有价值 |
| `max_genes: float = float('inf')` | YAML 中 `.inf` 可被 PyYAML 正确解析 |
| `model_dump()` 替代 `dict()` | Pydantic v2 兼容 |

---

## 6. 开发进度

### 已完成（11 个阶段）

| 阶段 | 内容 | 状态 |
|------|------|------|
| 1 | 5 个严重 Bug 修复（ssgsea 注册、CLI、urllib、OpenAI API、taxid） | ✅ |
| 2 | 10 个功能缺陷修复（配置同步、plot_all heatmap、导航栏等） | ✅ |
| 3 | 32 个单元测试补充 | ✅ |
| 4 | 代码质量优化（os.system→subprocess、清理缓存） | ✅ |
| 5 | 残余 Bug 修复（BUILTIN_TAXID_MAP、基因分隔符、API 路径等） | ✅ |
| 6 | 架构优化（延迟初始化、SpeciesInfo 重命名、死代码清理） | ✅ |
| 7 | 深度审计修复（EnrichmentResult 参数名、路径遍历、NaN pvalue 等） | ✅ |
| 8 | 最终审计修复（列名一致性、CORS 安全、导入清理、README 修正） | ✅ |
| 9 | 集成测试（v1 真实数据库，20 个测试） | ✅ |
| 10 | 推送到 GitHub | ✅ |
| 11 | **风格颜色系统重构**（R→Python 可视化迁移，Python原生 matplotlib/seaborn 实现，10种预置主题系统，36组颜色风格配置，暗色模式支持，PNG/PDF 双输出） | ✅ |

### 测试覆盖

- **总计**: 679 个测试全部通过
- **单元测试**: 64 个（算法、配置、URL 生成、物种检索、延迟初始化等）
- **集成测试**: 20 个（DatabaseManager 加载、完整分析流程、两遍扫描、结果导出、边界条件）
- **可视化端到端测试**: 26 个（柱状图、气泡图、GSEA/GSVA/ssGSEA 图表等）
- **GSEA/GSVA/ssGSEA 端到端测试**: 29 个（500 基因排序列表 + 6000×6 样本表达矩阵）

---

## 7. 待完成 / 可选改进

### 高优先级

- [x] **R 脚本→Python 迁移**: 已完成——使用 matplotlib/seaborn/cutecharts 替代 R/ggplot2，消除 R 依赖，新增主题系统和颜色配置
- [ ] **API 集成测试**: 使用 FastAPI TestClient 测试所有端点
- [ ] **CLI 集成测试**: 使用 subprocess 测试命令行入口

### 中优先级

- [ ] **CI/CD**: 添加 GitHub Actions（lint + test + 覆盖率）
- [ ] **任务过期清理**: API 服务器内存中任务自动清理机制
- [ ] **CORS 可配置化**: 通过环境变量 `ALLOWED_ORIGINS` 控制

### 低优先级

- [ ] **mypy 类型检查**: 放宽配置或补全所有类型注解
- [ ] **R 脚本注入防护**: 对嵌入 R 脚本的字符串值进行更严格的转义
- [ ] **文档**: API 文档（Swagger 已内置）、贡献指南

---

## 8. 开发指南

### 环境搭建

```bash
git clone https://github.com/zhangducsu/AllEnricher-v2.git
cd AllEnricher-v2
pip install -e ".[dev]"           # 开发模式安装（含测试依赖）
pip install -e ".[api,ai]"        # 含 API 和 AI 依赖
```

### 运行测试

```bash
pytest tests/ -v                   # 运行所有测试
pytest tests/test_integration.py -v -s  # 仅集成测试（含输出）
pytest tests/ -k "test_fisher"     # 按名称过滤
pytest --cov=allenricher           # 覆盖率报告
```

### 代码规范

- 类型注解：所有公共方法应有完整类型注解
- 文档字符串：所有公共类和方法使用中文 docstring
- 导入：仅导入实际使用的模块
- 文件末尾：PEP 8 要求有一个空行

### 添加新的分析方法

1. 在 `core/enrichment.py` 中创建 `EnrichmentMethodBase` 子类
2. 实现 `calculate_enrichment()` 和 `calculate_pvalue()` 方法
3. 在 `EnrichmentMethod` 枚举中添加新值
4. 在 `EnrichmentAnalyzer._get_method()` 中注册
5. 在 `cli.py` 的 `--method` choices 中添加
6. 编写对应的单元测试

### 添加新的数据库

1. 在 `database/manager.py` 中创建 `DatabaseBase` 子类
2. 实现 `parse()` 方法（读取 `{species}.XX2gene.tab.gz`）
3. 在 `DatabaseManager.DATABASE_CLASSES` 中注册
4. 准备测试数据库文件到 `tests/fixtures/database/hsa/`
5. 编写集成测试

### 数据库文件格式

v1 和 v2 共享相同的数据库文件格式（gzip 压缩 TSV）:
```
Gene    TermID1    TermID2    TermID3
BRCA1   1          0          1
TP53    0          1          1
```

`database_dir` 应指向包含 `{species}.XX2gene.tab.gz` 文件的目录。

---

## 9. v1 兼容性

- v1 项目位于 `/workspace/AllEnricher/`（Perl/R/Shell 实现）
- v2 可以直接加载 v1 的 `.tab.gz` 数据库文件
- v1 数据库版本: v20190612，覆盖 4 个物种（hsa/mmu/rno/ssc）
- 两遍扫描语义与 v1 完全一致
- 集成测试使用的测试数据截取自 v1 的 hsa 数据库

---

## 10. 常见陷阱

1. **`EnrichmentResult` 字段名**: `pvalue` 不是 `p_value`，`adjusted_pvalue` 不是 `q_value`，`gene_list` 不是 `genes`
2. **`to_dict()` 列名**: `Adjusted_P_Value` 不是 `Q_Value`，`Expected_Count` 不是 `Expected`
3. **`BUILTIN_TAXID_MAP` 与 `_load_builtin_species`**: 两处的 taxid 必须保持一致
4. **`Config` 与 `DEFAULT_CONFIG_YAML`**: 默认值必须保持同步
5. **`GSEA.__init__` min_size 默认值**: 应为 10（与 Config 一致），不是 15
6. **`max_genes` 类型**: 是 `float`（因为默认值是 `float('inf')`），不是 `int`
7. **API server `run_analysis`**: 是普通函数不是 async 函数
8. **Pydantic v2**: 使用 `model_dump()` 不是 `dict()`

---

## 11. 本次会话完成 — 风格颜色系统重构 (2026-05-28)

### 变更摘要

- **R→Python 可视化迁移**: 彻底移除 R/ggplot2/pheatmap/UpSetR 依赖，改用 Python 原生 matplotlib + seaborn + cutecharts 实现全部图表
- **模块化重构**: 将单一 `plotter.py` 拆分为 9 个模块文件（barplot, bubble, common_plots, color_config, gsea_plots, gsva_plots, plot_config, plot_theme, plotter），职责清晰、易于扩展
- **主题系统**: 新增 `plot_theme.py`，支持 10 种预置主题（Nature, Science, Cell, Lancet, NEJM, IEEE, ASCO, Modern, Minimal, Dark），以及自定义主题注册
- **颜色配置**: 新增 `color_config.py`，提供 36 组颜色风格（Category10~Category30 + 渐变色系 + 语义色系），支持暗色模式自动切换
- **图表配置**: 新增 `plot_config.py`，统一管理输出格式（PNG/PDF）、DPI（300 发表级）、字体大小等参数
- **输出格式增强**: 全部图表支持 PNG + PDF 双格式输出，满足不同场景需求
- **依赖消除**: 移除了对 R 运行时和 R 包（ggplot2/pheatmap/UpSetR）的依赖，简化部署
- **可视化测试**: 新增 `test_e2e_visualization.py`（26 个测试），覆盖全部图表类型的端到端验证
- **全量测试**: 679 个 pytest 测试全部通过
