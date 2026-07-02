<!--
================================================================================
AllEnricher v2 AGENTS.md -- AI Agent 协作开发上下文文件
================================================================================
版本: 基于 2026-07-03 项目状态
目的: 为 AI Agent（Codex/Claude/Trae Solo）提供完整的项目上下文，确保一致、高效的多会话协作开发
-->
# AllEnricher v2 -- Agent 协作开发指南

---

## 1. 项目概览

**AllEnricher** 是一个综合性的基因集功能富集分析工具，v2 为 Python 重构版（v1 为 Perl/R 经典版）。

- **仓库**: https://github.com/zd105/AllEnricher
- **位置**: `AllEnricher-v2/`
- **版本**: 2.3.0
- **许可证**: MIT
- **Python**: >=3.8
- **代码量**: ~12,000 行 Python（不含测试）
- **测试**: R 绘图桥接层、CLI R 分支、GSEA 核心算法、报告嵌入、API、本地全场景 E2E 和 live 外部服务 smoke 均已覆盖；最新归档见 `test_e2e_2026/99_runs/`
- **论文引用**: Zhang et al., BMC Bioinformatics, 2020

### v1 与 v2 关系

v1 位于 `AllEnricher-v1/`，是 Perl + R 经典实现，v2 必须兼容 v1 的数据库文件格式：
- v2 可直接加载 v1 构建的数据库文件（`.tab.gz` + `.2disc.gz`）
- 两遍扫描语义与 v1 `AllEnricher_v1.0.pl` 完全一致
- 统计结果（Fisher/Hypergeometric）与 v1 完全一致

---

## 2. 项目目录结构

```
AllEnricher-v2/
├── pyproject.toml              # PEP 518 项目配置：依赖、构建、工具链
├── README.md                   # 用户手册（中文）
├── DEVELOPMENT_PLAN.md         # 开发计划与 v1 代码分析
├── CLAUDE.md                   # 旧版 Agent 上下文文件
├── AGENTS.md                   # 本文件
├── config.example.yaml         # 配置文件示例
├── Dockerfile.test             # 测试用 Docker 镜像
├── allenricher/                # 【核心源码包】
│   ├── __init__.py             # 包入口，导出 5 个核心类
│   ├── cli.py                  # CLI 入口（12 个子命令）
│   ├── core/
│   │   ├── config.py           # Config dataclass + YAML 模板 + 校验
│   │   ├── enrichment.py       # 核心分析引擎：EnrichmentAnalyzer + 方法类层级
│   │   └── gsva.py             # GSVA 实现（3 种方法）
│   ├── database/
│   │   ├── manager.py          # DatabaseManager：数据库加载与调度中心
│   │   ├── builder.py          # 数据库构建调度器
│   │   ├── downloader.py       # 多线程数据下载
│   │   ├── gmt_generator.py    # GMT 文件自动生成
│   │   ├── custom_builder.py   # 自定义数据库构建器
│   │   ├── species_registry.py # 物种注册表（31,822 物种）
│   │   ├── species_lookup.py   # 物种检索
│   │   ├── ortholog_mapper.py  # 同源映射引擎
│   │   ├── version.py          # 版本管理
│   │   └── parsers/            # 各数据库解析器（14 个）
│   ├── visualization/
│   │   ├── plotter.py          # 统一调度入口
│   │   ├── barplot.py / bubble.py  # 基础图表
│   │   ├── gsea_plots.py / gsva_plots.py  # 方法专属图表
│   │   ├── common_plots.py     # 通用图表
│   │   ├── plot_theme.py       # 主题系统（10 种预设）
│   │   ├── color_config.py     # 配色系统（36 组方案）
│   │   └── plot_config.py      # 绘图参数配置
│   │   ├── r_plotter.py        # R 绘图桥接层（--use-r-plots 触发）
│   │   └── r_scripts/          # 11 个 R 脚本（发表级 GSEA 图表）
│   ├── report/
│   │   ├── generator.py        # HTML 报告生成器
│   │   └── templates/          # Jinja2 模板
│   ├── api/
│   │   ├── server.py           # FastAPI REST 服务（8 个端点）
│   │   └── static/             # Web 分析界面
│   └── ai/
│       └── interpreter.py      # AI 解读引擎（7 种后端）
├── tests/                      # 单元测试 + 集成测试（679 个）
├── database/                   # 数据库存储目录
│   ├── basic/                  # 通用基础数据
│   ├── organism/               # 物种特异性数据
│   └── versions.json           # 本地版本清单
├── docs/                       # 文档
├── examples/                   # 使用示例
├── test_data/                  # 测试数据
├── test_e2e_2026/              # 端到端测试
└── 00_input_data/              # 标准输入数据
```

---

## 3. 核心架构与数据流

### 3.1 分析管道

```
用户输入（基因列表 / 表达矩阵 / 排序基因）
    │
    ▼
CLI (cli.py)  ─── 解析参数、加载 Config
    │
    ▼
DatabaseManager ─── 加载物种数据库（.tab.gz 矩阵 + .2disc.gz 描述）
    │                    └── 自动查找 organism/v{date}/{species}/ 目录
    ▼
EnrichmentAnalyzer
    │   ├── 两遍扫描建立 gene_total / background_total
    │   ├── 执行统计检验（Fisher / Hypergeometric / GSEA / ssGSEA / GSVA）
    │   ├── 多重检验校正（BH / BY / Bonferroni / Holm / None）
    │   └── 按 p/q 值阈值过滤
    ▼
Plotter ─── 生成图表（PNG + PDF 双输出，300DPI）
    │
    ▼
ReportGenerator ─── 生成 HTML 学术报告
    │
    ▼
（可选）AIInterpreter ─── AI 驱动的结果解读
```

### 3.2 统计方法类层级

```
EnrichmentMethodBase (ABC)
├── ORA (过表达分析)
│   ├── FisherExactTest          # Fisher 精确检验
│   └── HypergeometricTest   # 超几何检验（复用 Fisher 的统计量计算）
├── GSEA                     # 基因集富集分析（含置换检验、NES、Leading Edge）
├── SSGSEA                   # 单样本 GSEA
└── GSVA                     # 基因集变异分析（Random Walk / PLAGE / Z-score）
```

### 3.3 数据库文件格式（v1 兼容）

- `{species}.{DB}2gene.tab.gz`: gzip 压缩 TSV，行为基因、列为 term，值为 0/1
- `{DB}2disc.gz`: gzip 压缩 TSV，TermID → TermName 映射
- `{species}.{DB}.gmt.gz`: GMT 格式基因集文件（build 时自动生成）

### 3.4 数据库构建流程

```
download -- 下载全物种通用原始数据 → database/basic/{type}/{date}/
    │
    ▼
build -- 提取指定物种数据 → database/organism/v{YYYYMMDD}/{species}/
    │   ├── 解析原始数据（parsers/）
    │   ├── 生成 .tab.gz + .2disc.gz
    │   └── 自动生成 .gmt.gz + build_manifest.json
    ▼
analyze -- 直接用 organism/ 目录下的数据运行分析
```

---

## 4. 支持的数据库清单

| 数据库 | 物种覆盖 | v1 兼容 | GSEA/GSVA | 完整度 |
|--------|---------|---------------|----------------------|-----------|--------|
| GO | 31,822 物种（需本地构建） | 是 | 是 | 完成 |
| KEGG | 所有 KEGG 物种 | 是 | 是 | 完成 |
| Reactome | 16 个模式物种 | 是 | 是 | 完成 |
| DO | 仅人类 | 是 | 是 | 完成 |
| DisGeNET | 仅人类 | 是 | 是 | 完成 |
| WikiPathways | 38 个模式物种 | 否 (v2 新增) | 是 | 完成 |
| TRRUST v2 | 人 + 小鼠 | 否 (v2 新增) | 是 | 完成 |
| ChEA3 | 仅人类（在线 API） | 否 (v2 新增) | 是 | 完成 |
| AnimalTFDB | 183 个动物物种 | 否 (v2 新增) | 是 | 完成 |
| hTFtarget | 仅人类 | 否 (v2 新增) | 是 | 完成 |
| Custom（用户自定义） | 不限 | 否 (v2 新增) | 是 | 完成 |

---

## 5. 可视化系统

### 5.1 主题系统（plot_theme.py）

每种主题控制 40+ 项参数：字体、字号、刻度方向（内/外）、网格线、边框、颜色等。

### 5.2 配色系统（color_config.py）

36 组配色方案：Category10~30、Paul Tol 系列、期刊风格（Nature/Science/Cell 等）、生物信息学专用、色盲友好等。
支持暗色模式自动切换。

### 5.3 图表类型总览

| 分析方法 | Python 原生图表 | R 图表（--use-r-plots） |
|---------|---------------|----------------------|
| ORA (Fisher/Hypergeometric) | barplot, bubble, dotplot | -- |
| GSEA | enrichment_curve, nes_barplot, dotplot, barplot, ridgeplot, emapplot, cnetplot, circos | dotplot, barplot, nes_plot, ridgeplot, emapplot, cnetplot, circos, enrichment, enrichment2, heatmap（10 类入口） |
| ssGSEA / GSVA | heatmap, group_comparison, dotplot, correlation | -- |
| 通用 | network, upset, volcano, method_comparison | -- |


R 图表系统通过 `--use-r-plots` CLI 参数启用，依赖 R + ggplot2/dplyr/tidyr/gridExtra/ggridges/ComplexHeatmap/circlize 等包。`enrichment`、`enrichment2` 和优化后的 `ridgeplot` 优先使用 ranked genes + GMT gene sets 生成的真实 running-ES 中间表；缺少必要数据时会明确跳过或降级为确定性摘要图，不生成模拟图。`heatmap` 现在优先展示 top GSEA 通路命中基因，无法匹配时才回退到表达矩阵中最变异基因。默认 R 图表已采用短标签、克制配色、线尾直标和更紧凑版式，适合报告与论文初稿继续精修。

---

## 6. CLI 命令总览

| 命令 | 功能 | 状态 |
|------|------|------|
| `analyze` | 运行富集分析（主命令） | 完成 |
| `download` | 下载全物种通用数据 | 完成 |
| `build` | 构建指定物种数据库（含自定义） | 完成 |
| `serve` | 启动 FastAPI Web 服务 | 完成 |
| `list` | 列出支持物种/数据库 | 完成 |
| `config` | 生成 YAML 配置文件 | 完成 |
| `check-update` | 远程版本更新检查（7 个数据源） | 完成 |
| `list-versions` | 查看已安装版本清单 | 完成 |
| `cleanup` | 清理旧版本数据 | 完成 |
| `list-species` | 列出支持物种 | 完成 |
| `query-species` | 按 taxid/KEGG 代码/拉丁名查询物种 | 完成 |
| `tf-enrich` | TF 调控网络富集分析 | 完成 |

---

## 7. 开发进度

### 7.1 已完成（按时间线）

| 日期 | 里程碑 | 关键变更 |
|------|--------|---------|
| 2026-05-15 | Bug 修复 | KEGG 层级修复、Term_Name 格式统一 |
| 2026-05-26 | GSEA/GSVA/ssGSEA | 三种方法实现 + 12+ 种图表 + GMT 生成 |
| 2026-05-27 | 物种注册表 | SpeciesRegistry（31,822 物种）、NCBI Taxonomy 集成 |
| 2026-05-28 | 风格配色重构 | R→Python 迁移、10 主题 + 36 配色 |
| 2026-05-28 | 版本管理 | check-update / list-versions / cleanup / build 血缘追溯 |
| 2026-05-30 | WikiPathways | 38 物种 + GPML 解析器 + NCBI Gene ID 转换 |
| 2026-05-30 | TF 调控网络 | TRRUST v2 + ChEA3 + 多库整合（TFMetaAnalyzer） |
| 2026-05-31 | 多物种 TF | AnimalTFDB（183 物种）+ hTFtarget + 同源映射 |
| 2026-06-13 | R-based GSEA 绘图修复 | 修复 R 参数路径、`%||%`、heatmap 输出格式、CLI heatmap/enrichment2 接入；enrichment 图改为真实 running-ES 中间表 |
| 2026-06-28 | R-based GSEA 绘图发表级优化 | 统一发表主题、标签压缩和色盲友好配色；ridgeplot 使用真实 hit-rank 分布；emapplot 改为通路重叠气泡矩阵；cnetplot 改为确定性通路-基因网络 |
| 2026-06-29 | R/GSEA 绘图真实性与报告回归 | 修复 signed ES/NES 方向、GSEA matrix 解包、running-ES 文件覆盖、报告只嵌入少量 R 图、Plotter 假成功路径；补齐 heatmap 通路基因筛选与 Matplotlib Agg 后端 |
| 2026-07-03 | v2.3.0 发布收尾 | 版本统一、package-data 打包 R/HTML/static 资源、config CLI 覆盖与 custom DB E2E 修复、全场景 E2E 归档 |

### 7.2 最新验证快照

- 已通过：`python -m py_compile allenricher/cli.py allenricher/visualization/r_plotter.py allenricher/visualization/plotter.py allenricher/visualization/plot_theme.py allenricher/visualization/barplot.py allenricher/visualization/bubble.py allenricher/visualization/common_plots.py allenricher/visualization/gsea_plots.py allenricher/visualization/gsva_plots.py allenricher/report/generator.py allenricher/core/enrichment.py`
- 已通过：`python -m pytest tests/test_enrichment.py tests/test_enrichment_extended.py tests/test_gsea_extended.py tests/test_visualization.py tests/test_report_integration.py tests/test_r_plotter.py tests/test_cli_r_plots.py -q`
- 最新验证：WSL2 已安装 R/R 包并完成 R 脚本冒烟、CLI `--use-r-plots` E2E 和视觉审计；local `20260702_235737` 为 `68 PASS / 5 EXPECTED_FAIL / 2 SKIP / 0 FAIL`，live `20260703_003638` 为 `2 PASS / 1 SKIP / 1 SKIPPED_MISSING_SECRET / 0 FAIL`，两者视觉审计均为 0 issue。
- 当前限制：download 用例按用户要求继续跳过；OpenAI live 用例在未提供 key 时标记为 `SKIPPED_MISSING_SECRET`。

### 7.3 待完成

| 优先级 | 任务 | 说明 |
|--------|------|------|
| 中 | CI/CD (GitHub Actions) | lint + test + coverage |
| 中 | 任务过期清理 | API 服务端内存中任务自动清理 |
| 中 | CORS 可配置化 | 通过环境变量 `ALLOWED_ORIGINS` 控制 |
| 低 | mypy 类型检查 | 放宽配置或补全全部类型注解 |
| 低 | API 文档补充 | Swagger 已内置，需补充贡献指南 |

---

## 8. 关键技术细节与约定

### 8.1 两遍扫描语义（与 v1 兼容）

1. **第一遍**: 遍历所有 term，收集至少命中一个输入基因的 term → 计算 `gene_total`
2. **第二遍**: 使用 `gene_total` 和 `background_total`（注释文件中所有基因数）计算 p 值

这是与 v1 结果完全一致的根本保证，不得修改此逻辑。

### 8.2 EnrichmentResult 字段名约定（易错点）

- Python dataclass 字段: `pvalue`, `adjusted_pvalue`, `gene_list`, `gene_count`
- `to_dict()` 输出列名（大写+下划线）: `P_Value`, `Adjusted_P_Value`, `Gene_Count`, `Genes`
- **注意**: 不是 `p_value`、不是 `q_value`、不是 `Q_Value`、不是 `Expected`

### 8.3 Config 与 YAML 一致性

- `Config` dataclass 的默认值必须与 `DEFAULT_CONFIG_YAML` 字符串模板保持一致
- `max_genes` 类型是 `float`（因为默认值是 `float('inf')`）
- `GSEA.__init__` 的 `min_size` 默认值是 10（与 Config 一致）

### 8.4 Pydantic v2

- 全部使用 `model_dump()` 替代 `dict()`
- API 服务使用 FastAPI + Pydantic v2

### 8.5 子进程调用

- 使用 `subprocess.run` 替代 `os.system`（安全性、超时控制、异常捕获）
- API 服务中 `run_analysis` 是普通函数（非 async），因为分析操作是 CPU 密集型

### 8.6 懒初始化模式

- `EnrichmentAnalyzer.method` 属性为懒初始化（property），首次访问时才根据 `config.method` 创建方法实例
- 这允许构造 Analyzer 后再修改 `config.method`

### 8.7 代码规范

- 所有公共方法应有完整类型注解
- 所有公共类和方法使用中文 docstring
- 文件末尾保留一个空行（PEP 8）
- 仅导入实际使用的模块
- 代码风格遵循 Black（line-length=100），isort 使用 Black 兼容配置

---

## 9. 测试概览

### 9.1 测试命令

```bash
pytest tests/ -v                   # 运行全部测试
pytest tests/ -k "test_fisher"     # 按名称过滤
pytest --cov=allenricher           # 覆盖率报告
docker build -t allenricher:v2 -f Dockerfile.test .   # Docker 环境测试
```

### 9.2 测试文件与数量

| 测试文件 | 数量 | 覆盖范围 |
|---------|------|---------|
| `test_enrichment.py` | 14 | 基础算法 |
| `test_enrichment_extended.py` | 50 | 扩展功能 |
| `test_database.py` | 8 | 数据库模块 |
| `test_cli.py` | 15 | CLI 接口 |
| `test_phase5.py` | 19 | API/AI/报告 |
| `test_e2e_gsea.py` | 7 | GSEA 端到端 |
| `test_e2e_ssgsea.py` | 8 | ssGSEA 端到端 |
| `test_e2e_gsva.py` | 14 | GSVA 端到端 |
| `test_e2e_visualization.py` | 26 | 可视化集成 |
| `test_gmt_generation_e2e.py` | 7 | GMT 生成 |
| `test_r_plotter.py` | 3 | R 绘图 Python 桥接层 |
| 其他专项测试 | ~500+ | TF/自定义数据库/WikiPathways 等 |
| **合计** | **~682+** | 以本地实际 pytest 输出为准 |

---

## 10. 安装与运行

```bash
# 开发模式安装（核心依赖）
pip install -e .

# 含 API + AI + 可视化
pip install -e ".[api,ai,visualization]"

# 含全部开发和测试工具
pip install -e ".[all]"

# 快速使用
allenricher analyze -i genes.txt -s hsa -d GO,KEGG -o results/
```

---

## 11. 多会话协作注意事项

1. **先读再写**: 修改任何文件前，先阅读该文件及其 `__init__.py` 的导入关系
2. **测试先行**: 新增功能前编写测试，修改逻辑后立即运行相关测试确认
3. **不分叉规范**: 遵循项目现有代码风格（Black、isort、中文 docstring）
4. **v1 兼容**: 修改分析引擎逻辑时必须保持两遍扫描语义不变
5. **不随意格式化**: 只修改必要代码，不格式化无关行
6. **检查点机制**: 多步骤任务每完成一个阶段，复盘已完成/已验证/待做事项
7. **显性暴露未完成**: 所有未完成、占位、伪代码、半成品必须在输出中明确标识
8. **优先使用 WSL2**: 本地任务优先在 WSL2 中运行

---

*最后更新: 2026-07-03（v2.3.0 发布收尾、R GSEA 绘图与全场景 E2E 验证）*
