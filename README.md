# AllEnricher v2.0

**基因集功能富集分析工具 - Python 重构版**

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

---

## 目录

- [简介](#简介)
- [功能特性](#功能特性)
- [系统要求](#系统要求)
- [安装方式](#安装方式)
- [快速开始](#快速开始)
- [详细使用说明](#详细使用说明)
  - [CLI 命令详解](#31-cli-命令详解)
  - [本地 Web 服务](#serve---启动本地-web-分析服务)
  - [自定义数据库构建](#自定义数据库构建)
  - [Python API 使用](#32-python-api-使用)
  - [GSEA/GSVA/ssGSEA 使用示例](#33-gseagsvasssgsea-使用示例)
- [数据库版本管理](#数据库版本管理)
- [测试指南](#测试指南)
- [输出文件说明](#输出文件说明)
- [配置文件](#配置文件)
- [常见问题](#常见问题)

---

## 简介

AllEnricher v2.0 是 AllEnricher v1.0 的 Python 重构版本，保留了 v1 的核心功能并增加了多项改进：

- **兼容 v1 数据库**：可直接使用 v1 构建的 GO、KEGG、Reactome、DO、DisGeNET 数据库
- **多种统计方法**：Fisher 精确检验、超几何检验、GSEA、ssGSEA、GSVA
- **多重检验校正**：BH、BY、Bonferroni、Holm
- **丰富可视化**：柱状图、气泡图、富集曲线图、热图、网络图等12+种发表级图表，Python原生实现，支持PNG/PDF双输出、6种学术风格主题、19种专业配色方案
- **本地 Web 服务**：基于 FastAPI 的交互式分析界面，浏览器打开即可使用
- **AI 智能解读**：支持 OpenAI、Claude、DeepSeek、GLM、MiniMax、Ollama
- **REST API**：完整的 RESTful API，支持程序化调用
- **交互式报告**：学术风格的 HTML 报告

---

## 功能特性

### 支持的数据库

| 数据库 | 描述 | 物种支持 |
|---------|------|----------|
| GO | Gene Ontology | 31822 种（需构建本地数据库） |
| KEGG | KEGG Pathway | 所有 KEGG 物种 |
| Reactome | Reactome Pathway | 16 种模式生物 |
| DO | Disease Ontology | 人类 |
| DisGeNET | 疾病-基因关联 | 人类 |

### 统计分析方法

| 方法 | 描述 | 适用场景 |
|------|------|----------|
| Fisher | Fisher 精确检验 | 标准富集分析 |
| Hypergeometric | 超几何检验 | Fisher 的替代方案 |
| **GSEA** | 基因集富集分析 | 排序基因列表（case vs control） |
| **ssGSEA** | 单样本 GSEA | 单样本通路活性评分 |
| **GSVA** | 基因集变异分析 | 多样本通路活性矩阵（3种方法变体） |

**GSEA/GSVA/ssGSEA 特性：**
- 🧬 **GSEA**: 基于置换检验的富集分析，支持NES标准化和Leading Edge识别
- 📊 **ssGSEA**: 单样本水平计算通路富集得分，适用于免疫浸润分析
- 🔄 **GSVA**: 三种计算方法（Random Walk/PLAGE/Z-score），适用于样本聚类和异质性分析
- 🎨 **发表级可视化**: 富集曲线图、热图、气泡图、网络图等12+种图表类型，Python原生matplotlib/seaborn实现，PNG/PDF双输出，6种学术风格主题，19种专业配色方案

### 可视化风格与配色系统

支持 **6 种学术风格主题** 和 **19 种专业配色方案**，覆盖主流期刊要求和展示场景。

**6 种学术风格：**

| 风格 | 适用场景 | 字体 | 字号 | 刻度方向 | 边框 |
|------|----------|------|------|----------|------|
| nature | Nature系列期刊投稿 | Helvetica | 8pt | 朝内 | 仅左、下 |
| science | Science期刊投稿 | Times New Roman | 9pt | 朝外 | 四边 |
| cell | Cell系列期刊投稿 | Arial | 9pt | 朝内 | 仅左、下 |
| colorblind | 色盲友好受众 | sans-serif | 10pt | 朝内 | 仅左、下 |
| presentation | 学术报告/PPT | sans-serif | 14pt | 朝外 | 四边+网格 |
| omicshare | 中文期刊投稿 | 微软雅黑 | 10pt | 朝内 | 四边 |

**19 种配色方案：**

- **Paul Tol 系列(8种)**: bright, high_contrast, vibrant, muted, medium_contrast, light, sunset, burga — 专为科学可视化优化
- **色盲友好(1种)**: okabe_ito — 8色可被常见色盲类型区分
- **科研期刊色板(6种)**: nature, science, cell, lancet, nejm, jama — 匹配期刊官方风格
- **生物信息学工具(2种)**: gsea, omicshare
- **中国风格(1种)**: china_style
- **默认色板(1种)**: default

**使用方式：**

```bash
# 指定风格和配色
allenricher analyze -i genes.txt --style nature --palette nature

# 学术报告风格
allenricher analyze -i genes.txt --style presentation --palette tol_vibrant

# 色盲友好
allenricher analyze -i genes.txt --style colorblind --palette okabe_ito
```

```python
from allenricher.visualization import PlotTheme

# 应用Nature风格
PlotTheme.apply('nature')

# 临时切换（with上下文）
with PlotTheme.context('presentation'):
    fig, ax = plt.subplots()
    # 此区域使用presentation风格
```

### AI 后端支持

- **OpenAI**: GPT-4、GPT-3.5
- **Anthropic**: Claude 3
- **DeepSeek**: deepseek-chat（国产，性价比高）
- **GLM**: glm-4（智谱 AI，中文能力强）
- **MiniMax**: abab6.5s-chat
- **Ollama**: 本地部署的开源模型
- **Mock**: 测试用模拟后端

---

## 系统要求

### 必需依赖

- Python 3.8+

### 可选依赖

- FastAPI + Uvicorn（API 服务）
- OpenAI/Anthropic Python 包（AI 解读）
- pytest（运行测试）

### Python 可视化依赖

可视化模块使用 Python 原生实现，自动安装以下依赖（随 `pip install -e .` 一起安装）：
- matplotlib>=3.4.0 — 基础绘图引擎
- seaborn>=0.11.0 — 统计图形美化
- cutecharts>=0.2.0 — 轻量图表补充

---

## 安装方式

### 方式一：Docker（推荐，最简单）

```bash
# 克隆项目
git clone https://github.com/zd105/AllEnricher.git
cd AllEnricher/AllEnricher-v2

# 构建 Docker 镜像（包含 Python + 所有依赖）
docker build -t allenricher:v2 -f Dockerfile.test .

# 运行测试
docker run --rm -v "$(pwd):/workspace" -w /workspace allenricher:v2 bash -c "python3 -m pytest tests/ -v"
```

> 注意：Dockerfile.test 已包含 Python、pytest 等所有依赖。

### 方式二：pip 安装

```bash
# 克隆项目
git clone https://github.com/zd105/AllEnricher.git
cd AllEnricher/AllEnricher-v2

# 安装基础依赖
pip install -e .

# 安装可选依赖（API 服务）
pip install -e ".[api]"

# 安装 AI 解读依赖
pip install -e ".[ai]"

# 安装所有依赖
pip install -e ".[all]"
```

### 方式三：本地开发环境

```bash
# 1. 安装 Python 依赖
pip install numpy pandas scipy statsmodels pyyaml requests tqdm openpyxl

# 2. 安装 pytest（测试用）
pip install pytest pytest-cov

# 3. 安装项目
pip install -e .
```

---

## 快速开始

### 1. 准备输入文件

一个简单的基因列表文件（每行一个基因符号）：

```text
BRCA1
TP53
EGFR
MYC
KRAS
BRAF
PTEN
RB1
CDK4
MDM2
```

保存为 `genes.txt`。

### 2. 使用 CLI 进行分析

```bash
# 基本用法（使用 v1 构建的数据库）
allenricher analyze \
    -i genes.txt \
    -s hsa \
    -d GO,KEGG \
    -o results/

# 指定数据库目录（使用 v1 已有的数据库）
allenricher analyze \
    -i genes.txt \
    -s hsa \
    -d GO,KEGG \
    -o results/ \
    --database-dir /path/to/v1/database/organism/v20190612/hsa

# 指定 p 值阈值
allenricher analyze \
    -i genes.txt \
    -s hsa \
    -d GO \
    -p 0.01 \
    -q 0.01 \
    -o results/
```

### 3. 查看输出结果

```bash
# 输出目录结构
results/
├── GO_enrichment.tsv      # GO 富集结果
├── KEGG_enrichment.tsv    # KEGG 富集结果
├── plots/                 # 可视化图表（PNG + PDF 双输出）
│   ├── GO_barplot.png
│   ├── GO_barplot.pdf
│   ├── GO_bubble.png
│   ├── GO_bubble.pdf
│   ├── KEGG_barplot.png
│   ├── KEGG_barplot.pdf
│   ├── KEGG_bubble.png
│   └── KEGG_bubble.pdf
└── report.html            # HTML 交互报告
```

---

## 详细使用说明

### 3.1 CLI 命令详解

#### analyze - 富集分析

```bash
allenricher analyze [选项]

必选参数：
  -i, --input FILE          输入基因列表文件（每行一个基因符号）

可选参数：
  -s, --species TEXT        物种代码（默认：hsa）
                            支持：hsa, mmu, rno, ssc, dre, cel, dme, ath, ...
  -d, --databases TEXT      数据库列表，逗号分隔（默认：GO,KEGG）
                            支持：GO, KEGG, Reactome, DO, DisGeNET
  -o, --output DIR          输出目录（默认：./results）
  -m, --method TEXT         统计方法（默认：fisher）
                            fisher, hypergeometric, gsea, ssgsea, gsva
  -e, --expression-matrix   表达矩阵文件（ssGSEA/GSVA必需，TSV格式）
  -r, --ranked-genes        排序基因列表文件（GSEA必需，含权重）
  --gmt FILE                GMT格式基因集文件（GSEA/ssGSEA/GSVA可选）
  --plot-types LIST         可视化图表类型，逗号分隔
                            GSEA: enrichment,nes_barplot,dotplot
                            ssGSEA/GSVA: heatmap,group_comparison,dotplot,correlation
  --plot-format FORMAT      图表格式（png/pdf/svg，默认png；支持同时输出PNG+PDF，用逗号分隔如"png,pdf"）
  --plot-dpi INT           图表分辨率（默认300）
  -c, --correction TEXT     多重检验校正方法（默认：BH）
                            BH, BY, bonferroni, holm, none
  -p, --pvalue FLOAT        P 值阈值（默认：0.05）
  -q, --qvalue FLOAT        Q 值阈值（默认：0.05）
  -n, --min-genes INT       每条目最小基因数（默认：2）
  -j, --jobs INT            并行任务数（默认：1）
  --database-dir DIR        数据库目录（必选，如使用 v1 数据库）
  --config FILE             YAML 配置文件
  --no-plot                跳过绘图
  --no-report              跳过 HTML 报告
  --verbose                详细日志输出
```

#### list - 列出资源

```bash
# 列出支持的物种
allenricher list species

# 列出支持的数据库
allenricher list databases
```

#### config - 生成配置

```bash
# 生成默认配置文件
allenricher config -o my_config.yaml

# 编辑配置文件后使用
allenricher analyze --config my_config.yaml -i genes.txt
```

#### download - 下载全体物种通用数据

```bash
# 下载 GO 全体物种基础数据（gene2go.gz + gene_info.gz + go-basic.obo）
allenricher download -d go

# 下载 Reactome 全体物种基础数据（NCBI2Reactome + gene_info.gz）
allenricher download -d reactome

# 同时下载多个
allenricher download -d go,reactome

# 可选参数
allenricher download -d go,reactome,do,disgenet --database-dir ./database
```

**下载加速选项：**

```bash
# 多线程加速（默认 4 线程，大文件 >100MB 自动启用）
allenricher download -d go --workers 8

# 禁用多线程（网络不稳定时使用）
allenricher download -d go --no-multi-thread

# 跳过完整性校验（加快下载，但不推荐）
allenricher download -d go --no-verify

# 强制重新下载（覆盖已存在文件）
allenricher download -d go --force
```

> **注意**：download 下载的是全体物种的通用原始数据，存入 `database/basic/{type}/{date}/`。
> 后续需要 build 命令从中提取指定物种的数据。
> 
> **特性**：大文件自动多线程分片下载 | 镜像源自动切换 | gzip 完整性校验 | 断点续传

#### build - 构建指定物种数据库

```bash
# 为人类 (hsa, taxid=9606) 构建 GO 和 Reactome 数据库
allenricher build -s hsa -t 9606 -d GO,Reactome

# 输出到 database/organism/v{date}/hsa/
# 然后可用于 analyze
allenricher analyze -i genes.txt -s hsa --database-dir database/organism/v{date}/hsa/
```

#### 自定义数据库构建

`allenricher build` 支持通过用户提供的注释文件构建自定义数据库，并自动生成 GMT 文件，可直接用于后续的富集分析（包括 GSEA/ssGSEA/GSVA）。

**支持的注释文件格式（TSV，Tab 分隔）：**

| 格式 | 列定义 | 示例 |
|------|--------|------|
| 四列（带层级） | `gene<TAB>term_id<TAB>term_name<TAB>hierarchy` | `BRCA1<TAB>TERM001<TAB>Cell Cycle<TAB>Biology\|Cell Biology\|Cell Cycle` |
| 三列 | `gene<TAB>term_id<TAB>term_name` | `BRCA1<TAB>TERM001<TAB>Cell Cycle` |
| 两列 | `gene<TAB>term` | `BRCA1<TAB>Cell Cycle`（term 同时作为 term_id 和 term_name） |

> 默认自动检测文件格式（根据首行列数），也可通过 `--annot-format` 手动指定。

**CLI 用法示例：**

```bash
# 基本用法（自动检测格式）
allenricher build \
    --species hsa \
    --taxid 9606 \
    --custom-annot annotation.txt \
    --custom-db-name MyDB

# 指定四列格式和层级分隔符
allenricher build \
    --species hsa \
    --taxid 9606 \
    --custom-annot annotation.txt \
    --custom-db-name MyDB \
    --annot-format four_column \
    --hierarchy-sep "|"
```

**Python API 用法：**

```python
from allenricher.database.custom_builder import CustomDatabaseBuilder

# 创建构建器
builder = CustomDatabaseBuilder(root_dir="./database")

# 从注释文件构建自定义数据库
output_dir = builder.build_from_annotation(
    annotation_file="annotation.txt",
    species="hsa",
    taxid=9606,
    db_name="MyDB"
)

# 指定格式和层级分隔符
output_dir = builder.build_from_annotation(
    annotation_file="annotation.txt",
    species="hsa",
    taxid=9606,
    db_name="MyDB",
    format_type="four_column",
    hierarchy_separator="|"
)

# 构建完成后即可用于分析
# allenricher analyze -i genes.txt -s hsa -d MyDB --database-dir <output_dir>
```

**输出文件说明：**

构建完成后，在 `database/organism/v{YYYYMMDD}/{species}/` 目录下生成以下 3 个文件：

| 文件 | 说明 |
|------|------|
| `{species}.{db_name}2gene.tab.gz` | 基因-条目 0/1 矩阵（TSV 格式，gzip 压缩） |
| `{db_name}2disc.gz` | 条目描述文件（含层级信息，gzip 压缩） |
| `{species}.{db_name}.gmt.gz` | GMT 基因集文件（**程序自动生成，用户无需提供**） |

> **关键说明**：GMT 文件由程序在构建时自动生成，用户只需提供注释文件即可。生成的自定义数据库可直接用于 `allenricher analyze` 的 `-d` 参数，也支持 GSEA/ssGSEA/GSVA 分析。

#### serve - 启动本地 Web 分析服务

```bash
# 安装 API 依赖
pip install -e ".[api]"

# 启动服务（默认端口 8000）
allenricher serve --port 8000 --host 0.0.0.0

# 浏览器打开
# http://localhost:8000
```

启动后，用户可以在浏览器中打开 `http://localhost:8000`，通过 **交互式 Web 界面** 完成富集分析，无需编写任何代码。

**Web 界面功能：**
- 📝 基因列表输入（文本粘贴或文件上传）
- 🧬 物种选择（支持所有已构建数据库的物种）
- 🗃️ 数据库勾选（GO、KEGG、Reactome、DO、DisGeNET）
- ⚙️ 参数配置（P值阈值、校正方法、最小基因数等）
- 📊 结果展示（交互式表格，支持排序和搜索）
- 📥 结果下载（TSV结果、PDF图表、HTML报告）
- 🔄 异步任务（提交后实时显示进度）

**API 端点一览：**

| 端点 | 方法 | 描述 |
|------|------|------|
| `/` | GET | Web 分析界面 |
| `/api/species` | GET | 获取支持的物种列表 |
| `/api/databases` | GET | 获取可用数据库信息 |
| `/api/analyze` | POST | 提交富集分析任务 |
| `/api/upload` | POST | 上传基因列表文件并分析 |
| `/api/status/{job_id}` | GET | 查询任务状态和进度 |
| `/api/results/{job_id}` | GET | 获取分析结果（JSON/TSV） |
| `/api/results/{job_id}/plot` | GET | 下载可视化图表 |
| `/api/results/{job_id}/report` | GET | 下载 HTML 报告 |
| `/api/jobs/{job_id}` | DELETE | 删除任务 |
| `/docs` | GET | Swagger API 文档 |
| `/redoc` | GET | ReDoc API 文档 |

**API 调用示例：**

```bash
# 提交分析任务
curl -X POST "http://localhost:8000/api/analyze" \
  -H "Content-Type: application/json" \
  -d '{
    "genes": ["BRCA1", "TP53", "EGFR", "MYC", "KRAS"],
    "species": "hsa",
    "databases": ["GO", "KEGG"],
    "method": "fisher",
    "pvalue_cutoff": 0.05,
    "qvalue_cutoff": 0.05
  }'
# 返回: {"job_id": "xxx-xxx-xxx", "status": "pending"}

# 查询任务状态
curl "http://localhost:8000/api/status/xxx-xxx-xxx"
# 返回: {"job_id": "...", "status": "completed", "progress": 1.0, ...}

# 获取结果
curl "http://localhost:8000/api/results/xxx-xxx-xxx"
# 返回: {"GO": [...], "KEGG": [...], ...}

# 下载 HTML 报告
curl -o report.html "http://localhost:8000/api/results/xxx-xxx-xxx/report"
```

### 3.2 Python API 使用

```python
from allenricher import EnrichmentAnalyzer, Config, DatabaseManager

# 1. 配置
config = Config(
    species="hsa",
    databases=["GO", "KEGG"],
    method="fisher",
    pvalue_cutoff=0.05,
    qvalue_cutoff=0.05,
    min_genes=2
)

# 2. 加载数据库
db_manager = DatabaseManager(
    db_dir="./database/organism/v20190612/hsa",
    species="hsa"
)
db_manager.load_databases(["GO", "KEGG"])

# 3. 获取数据
gene_set = {"BRCA1", "TP53", "EGFR", "MYC", "KRAS"}
background = db_manager.get_background_genes()
database_data = db_manager.get_all_term_data()

# 4. 运行分析
analyzer = EnrichmentAnalyzer(config)
results = analyzer.run_analysis(
    gene_set=gene_set,
    background_set=background,
    database_data=database_data,
    parallel=False
)

# 5. 保存结果
analyzer.save_results("./output")

# 6. 生成报告
from allenricher.report.generator import ReportGenerator
report_gen = ReportGenerator("./output")
report_gen.generate(results, "./output/report.html", gene_list=list(gene_set))
```

### 3.3 GSEA/GSVA/ssGSEA 使用示例

#### GSEA 分析（排序基因列表）

```bash
# 准备排序基因列表文件（两列：gene, weight）
# ranked_genes.tsv:
# gene<TAB>weight
# BRCA1<TAB>2.5
# TP53<TAB>-1.8
# ...

# 运行GSEA分析
allenricher analyze \
    --method gsea \
    -r ranked_genes.tsv \
    -s hsa \
    --gmt hsa.KEGG.gmt.gz \
    --plot-types enrichment,nes_barplot,dotplot \
    --plot-format png,pdf \
    -o gsea_results/
```

#### ssGSEA 分析（单样本通路活性）

```bash
# 准备表达矩阵（基因×样本，TSV格式）
# expression.tsv:
# gene<TAB>Sample1<TAB>Sample2<TAB>...
# BRCA1<TAB>8.5<TAB>9.2<TAB>...
# ...

# 运行ssGSEA分析
allenricher analyze \
    --method ssgsea \
    -e expression.tsv \
    -s hsa \
    --gmt hsa.GO.gmt.gz \
    --plot-types heatmap,group_comparison \
    --groups "Tumor:Sample1,Sample2,Sample3;Normal:Sample4,Sample5,Sample6" \
    -o ssgsea_results/
```

#### GSVA 分析（三种方法变体）

```bash
# GSVA (Random Walk - 默认)
allenricher analyze \
    --method gsva \
    -e expression.tsv \
    -s hsa \
    --gmt hsa.Reactome.gmt.gz \
    --plot-types heatmap,correlation \
    -o gsva_results/

# GSVA (PLAGE方法)
allenricher analyze \
    --method gsva \
    -e expression.tsv \
    --gmt hsa.Reactome.gmt.gz \
    --gsva-method plage \
    -o gsva_plage_results/

# GSVA (Z-score方法)
allenricher analyze \
    --method gsva \
    -e expression.tsv \
    --gmt hsa.Reactome.gmt.gz \
    --gsva-method zscore \
    -o gsva_zscore_results/
```

---

## 数据库版本管理

v2.1.0 新增完整的数据库版本管理体系，支持远程更新检测、本地版本追溯、构建血缘追踪、版本锁定和冗余清理。

### check-update - 检查远程数据源更新

```bash
# 检查所有数据源是否有更新
allenricher check-update

# 指定数据库目录
allenricher check-update --database-dir ./database

# JSON 格式输出
allenricher check-update --json
```

支持检测 7 个数据源：NCBI gene2go、NCBI gene_info、GO Ontology、EBI GOA、KEGG、Reactome、NCBI Taxonomy。

### list-versions - 查看本地已安装版本

```bash
# 查看所有已安装的基础数据和物种数据库版本
allenricher list-versions

# 查看构建血缘追踪（依赖链 + 源数据版本）
allenricher list-versions --lineage

# JSON 格式输出
allenricher list-versions --json
```

### cleanup - 清理旧版本

```bash
# 预览将删除的旧版本（保留最新 2 个版本）
allenricher cleanup --dry-run

# 保留最新 1 个版本
allenricher cleanup --keep 1

# 实际执行清理
allenricher cleanup --keep 2
```

### download --force - 强制重新下载

```bash
# 下载前自动检查更新，无更新时跳过
allenricher download -d go

# 强制重新下载，跳过更新检查
allenricher download -d go --force
```

### analyze --use-version - 使用指定版本分析

```bash
# 使用指定版本的数据库进行分析（确保结果可复现）
allenricher analyze -i genes.txt -s hsa -d GO --use-version v20260515 -o results/
```

分析结果（TSV/HTML）自动嵌入数据库版本信息，确保分析结果可追溯和可复现。

---

## 测试指南

### 4.1 运行所有测试

```bash
# 在项目根目录
cd AllEnricher-v2

# 运行所有测试
pytest tests/ -v

# 带覆盖率
pytest tests/ --cov=allenricher --cov-report=html

# 运行特定测试文件
pytest tests/test_enrichment.py -v
pytest tests/test_database.py -v
pytest tests/test_cli.py -v
```

### 4.2 Docker 中运行测试

```bash
# 构建镜像
docker build -t allenricher:v2 -f Dockerfile.test .

# 运行测试
docker run --rm \
    -v "$(pwd):/workspace" \
    -w /workspace \
    allenricher:v2 \
    python3 -m pytest tests/ -v
```

### 4.3 使用 v1 数据进行端对端测试

```bash
# v1 数据库路径
V1_DB="/path/to/v1/database/organism/v20190612/hsa"

# v1 示例基因列表
V1_EXAMPLE="/path/to/v1/example/example.glist"

# 运行分析
allenricher analyze \
    -i "$V1_EXAMPLE" \
    -s hsa \
    -d GO,KEGG \
    -o ./test_results \
    --database-dir "$V1_DB"

# 对比 v1 结果
diff ./test_results/GO_enrichment.tsv "$V1_EXAMPLE/../fisher/Q0.05/example.glist.GO.xls"
```

### 4.4 测试文件说明

| 测试文件 | 测试内容 |
|----------|----------|
| `test_enrichment.py` | 核心富集分析（14 个测试） |
| `test_enrichment_extended.py` | 扩展功能测试（50 个测试） |
| `test_database.py` | 数据库模块（8 个测试） |
| `test_cli.py` | 命令行接口（15 个测试） |
| `test_phase5.py` | API/AI/报告模块（19 个测试） |
| `test_e2e_gsea.py` | GSEA端到端测试（7 个测试） |
| `test_e2e_ssgsea.py` | ssGSEA端到端测试（8 个测试） |
| `test_e2e_gsva.py` | GSVA端到端测试（14 个测试） |
| `test_e2e_visualization.py` | 可视化集成测试（26 个测试，覆盖全部图表类型） |
| `test_gmt_generation_e2e.py` | GMT文件生成测试（8 个测试） |

---

## 输出文件说明

### 5.1 TSV 结果文件

```
Term_ID     Term_Name                   Gene_Count  Rich_Factor  P_Value        Adjusted_P_Value  Genes           Term_URL
GO:0005576  extracellular region        10          0.0523        1.23e-05      0.00234           GENE1;GENE2;...  https://amigo.geneont...
GO:0051301  cell division                8           0.0418        4.56e-04      0.01823           GENE3;GENE4;...  https://amigo.geneont...
```

**列说明：**

| 列名 | 说明 |
|------|------|
| Term_ID | GO/Pathway 标识符 |
| Term_Name | GO/Pathway 名称 |
| Gene_Count | 条目中的基因数 |
| Rich_Factor | 富集比率（Obs/Exp） |
| P_Value | 原始 P 值 |
| Adjusted_P_Value | 校正后的 P 值（Q 值） |
| Genes | 富集基因列表（分号分隔） |
| Term_URL | GO/Pathway 网页链接 |

### 5.2 PNG + PDF 可视化

可视化模块使用 Python 原生 matplotlib/seaborn 实现，支持 PNG 和 PDF 双格式输出（300DPI 发表级分辨率）。

- **柱状图**：横轴为 -log10(Q-value)，纵轴为条目名称，按类别着色，支持 10 种学术主题风格
- **气泡图**：横轴为 Rich Factor，纵轴为条目名称，点大小为基因数，颜色为 -log10(Q-value)

### 5.3 HTML 报告

学术风格的交互式报告，包含：

- 分析摘要统计
- 交互式数据表格（支持排序、搜索、下载）
- 可视化图表链接
- 可选的 AI 解读内容

### 5.4 GSEA/GSVA/ssGSEA 可视化图表

**GSEA 专属图表：**

| 图表类型 | 描述 | 文件（PNG + PDF 双输出） |
|----------|------|------|
| 富集曲线图 | 三面板图：Running ES曲线、基因位置标记、基因排序度量 | `*_enrichment.{png,pdf}` |
| NES条形图 | 水平条形图展示各通路NES值，按显著性排序 | `*_nes_barplot.{png,pdf}` |
| 气泡图 | 气泡大小表示基因数，颜色表示显著性 | `*_dotplot.{png,pdf}` |

**ssGSEA/GSVA 专属图表：**

| 图表类型 | 描述 | 文件（PNG + PDF 双输出） |
|----------|------|------|
| 通路活性热图 | 样本×通路的活性得分热图，支持聚类 | `*_heatmap.{png,pdf}` |
| 组间比较图 | 箱线图/小提琴图比较不同组别通路活性 | `*_group_comparison.{png,pdf}` |
| 样本相关性热图 | 样本间通路活性的相关性矩阵 | `*_correlation.{png,pdf}` |

**通用图表：**

| 图表类型 | 描述 | 文件（PNG + PDF 双输出） |
|----------|------|------|
| 通路网络图 | 基于基因重叠的通路关系网络 | `enrichment_network.{png,pdf}` |
| 火山图 | NES vs -log10(pvalue)的散点图 | `*_volcano.{png,pdf}` |
| 方法比较图 | 不同方法结果的相关性散点图 | `method_comparison.{png,pdf}` |

**图表配置：**

```bash
# 指定图表类型
--plot-types enrichment,heatmap,dotplot

# 指定输出格式（支持PNG/PDF/SVG，可同时输出多个）
--plot-format png,pdf

# 指定分辨率（发表级300DPI）
--plot-dpi 300

# 选择主题风格（默认modern）
--plot-theme nature

# 选择颜色风格（默认Set2）
--color-style Nature

# 暗色模式
--dark-mode

# 生成所有可用图表
--plot-types all

# 支持的 --plot-theme 选项：
# nature, science, cell, lancet, nejm, ieee, asco, modern, minimal, dark

# 支持的 --color-style 选项：
# Category10~Category30, Nature, Science, Cell, Lancet, NEJM, IEEE, ASCO, Set1~Set3, Pastel1~Pastel2, Dark2, Accent
```

---

## 配置文件

### 6.1 YAML 配置示例

```yaml
# 基本配置
species: "hsa"
databases:
  - "GO"
  - "KEGG"
  - "Reactome"

# 分析参数
method: "fisher"
correction: "BH"
pvalue_cutoff: 0.05
qvalue_cutoff: 0.05
min_genes: 2

# GSEA 参数
gsea_min_size: 10
gsea_max_size: 500

# 可视化
top_terms: 20
plot_formats:
  - "png"
  - "pdf"
plot_theme: "modern"          # 主题：nature/science/cell/lancet/nejm/ieee/asco/modern/minimal/dark
color_style: "Set2"           # 颜色风格：Category10~30, Nature, Science, Cell 等
dark_mode: false              # 暗色模式

# 性能
n_jobs: 4

# AI 解读（可选）
ai_interpretation: false
ai_backend: "mock"  # mock, openai, claude, deepseek, glm, minimax, ollama
ai_model: "gpt-4"

# 数据库路径
database_dir: "/path/to/v1/database/organism/v20190612/hsa"
```

### 6.2 环境变量

| 变量 | 说明 | 用于 |
|------|------|------|
| `OPENAI_API_KEY` | OpenAI API 密钥 | AI 解读 |
| `ANTHROPIC_API_KEY` | Anthropic API 密钥 | AI 解读 |
| `DEEPSEEK_API_KEY` | DeepSeek API 密钥 | AI 解读 |
| `GLM_API_KEY` | 智谱 AI API 密钥 | AI 解读 |
| `MINIMAX_API_KEY` | MiniMax API 密钥 | AI 解读 |
| `MINIMAX_GROUP_ID` | MiniMax Group ID | AI 解读 |

---

## 常见问题

### Q1: 如何使用 v1 构建的数据库？

```bash
# 方式一：命令行指定
allenricher analyze -i genes.txt --database-dir /path/to/v1/database/organism/v20190612/hsa

# 方式二：配置文件
allenricher analyze --config config.yaml -i genes.txt
# config.yaml 中设置 database_dir
```

### Q2: 提示 "Database directory not found"

确保：
1. 数据库目录路径正确
2. 目录中存在 `.tab.gz` 和 `.2disc.gz` 文件

```bash
# 检查目录内容
ls /path/to/v1/database/organism/v20190612/hsa/
# 应该看到：hsa.GO2gene.tab.gz, GO2disc.gz, hsa.kegg2gene.tab.gz, ...
```

### Q3: 可视化图表生成失败

确保已安装 Python 可视化依赖：

```bash
pip install matplotlib seaborn cutecharts
```

### Q4: 如何构建自己的数据库？

```bash
# 需要准备 NCBI gene_info.gz 文件
# 从 https://ftp.ncbi.nlm.nih.gov/gene/DATA/ 下载

allenricher build \
    -s hsa \
    -t 9606 \
    -d GO \
    --database-dir ./my_database \
    --gene-info /path/to/gene_info.gz
```

### Q5: AI 解读如何使用？

```bash
# 使用 Mock（无需配置）
allenricher analyze -i genes.txt --ai mock

# 使用 OpenAI
export OPENAI_API_KEY=sk-xxx
allenricher analyze -i genes.txt --ai openai --ai-key $OPENAI_API_KEY

# 使用 DeepSeek（国产）
export DEEPSEEK_API_KEY=sk-xxx
allenricher analyze -i genes.txt --ai deepseek --ai-key $DEEPSEEK_API_KEY
```

### Q6: 测试失败怎么办？

1. 确保所有依赖已安装
2. 使用 Docker 环境排除环境问题：

```bash
docker build -t allenricher:v2 -f Dockerfile.test .
docker run --rm -v "$(pwd):/workspace" -w /workspace allenricher:v2 \
    python3 -m pytest tests/test_enrichment.py -v
```

---

## 许可证

MIT License

## 联系方式

- **GitHub Issues**: https://github.com/zd105/AllEnricher/issues
- **邮箱**: allenricher@example.com

## 引用

如果 AllEnricher 对您的研究有帮助，请引用：

```bibtex
@article{zhang2020allenricher,
  title={AllEnricher: a comprehensive gene set function enrichment tool for both model and non-model species},
  author={Zhang, Du and Hu, Qian and Liu, Xiang and others},
  journal={BMC Bioinformatics},
  volume={21},
  pages={106},
  year={2020}
}
```

---

## 最近更新 (v2.2.0)

### 2026-05-28: 风格颜色系统重构

- **R→Python 可视化迁移**: 彻底移除 R/ggplot2/pheatmap/UpSetR 依赖，改用 Python 原生 matplotlib + seaborn 实现全部图表
- **模块化重构**: 将单一 `plotter.py` 拆分为 9 个模块（barplot, bubble, common_plots, color_config, gsea_plots, gsva_plots, plot_config, plot_theme, plotter），职责清晰、易于扩展
- **6 种学术风格主题**: nature（Helvetica, 8pt, 刻度朝内, 仅左下边框）、science（Times New Roman, 9pt, 刻度朝外, 四边边框）、cell（Arial, 9pt, 刻度朝内）、colorblind（sans-serif, 10pt）、presentation（14pt, 网格线, 适合投影）、omicshare（微软雅黑, 中文友好），每种风格独立控制字体、字号、粗细、边框、内外刻度、网格线等 40 项参数
- **19 种专业配色方案**: Paul Tol 系列 8 种（bright/vibrant/muted/light/sunset 等）、色盲友好 okabe_ito 1 种、科研期刊风格 6 种（nature/science/cell/lancet/nejm/jama）、生物信息学 2 种（gsea/omicshare）、中国风格 1 种（china_style）、默认 1 种
- **GO/KEGG 分类颜色动态生成**: 不再硬编码分类颜色，根据当前配色方案自动从色板提取
- **300DPI + PNG/PDF 双输出**: 全部图表支持 300 DPI 发表级分辨率，同时输出 PNG（位图）和 PDF（矢量图）
- **气泡图三色渐变**: palette 首尾两个不连续颜色 + 白色构成三色渐变，灰色尺寸图例与 Q 值色条并列展示
- **条形图 GO 分类前缀清理**: 纵坐标自动移除 BP/CC/MF 分类前缀，分类信息由颜色图例提供
- **CLI 新增 --style/--palette 参数**: 支持命令行直接指定图表风格和配色方案
- **依赖简化**: 移除了对 R 运行时的依赖，部署更简单
- **全量测试通过**: 679 个 pytest 测试全部通过

### 2026-05-28: 数据库版本管理系统

- **🔄 远程更新检测**: 新增 `check-update` 命令，支持 7 个数据源的远程版本检测（NCBI gene2go/gene_info、GO Ontology、KEGG、Reactome、EBI GOA、NCBI Taxonomy）
- **📋 本地版本清单**: 新增 `list-versions` 命令，查看已安装的基础数据和物种数据库版本，支持 `--lineage` 血缘追踪和 `--json` 输出
- **🧹 旧版本清理**: 新增 `cleanup` 命令，支持 `--dry-run` 预览和 `--keep N` 保留策略
- **📥 智能下载**: `download` 命令自动检查远程更新，无更新时跳过；新增 `--force` 强制重新下载
- **🔒 版本锁定**: `analyze` 命令新增 `--use-version` 参数，支持指定任意已安装版本进行分析，确保结果可复现
- **📝 构建血缘追踪**: 每次 `build` 自动生成 `build_manifest.json`，记录完整的源数据依赖链和版本号
- **📊 分析结果版本记录**: TSV 输出文件头部自动嵌入版本注释（AllEnricher 版本、数据库版本、源数据版本），HTML 报告动态显示版本信息
- **🧪 测试覆盖**: 版本管理模块 25 个单元测试全部通过，全量 651 个测试 0 失败

### 2026-05-27: 统一物种注册表与数据库管理优化

- **🧬 统一物种注册表**: 新增 `SpeciesRegistry` 模块，支持 31,822 个物种的统一查询和管理
- **🔍 CLI 物种查询**: 新增 `list-species` 和 `query-species` 命令，支持按 taxid/KEGG 代码/拉丁名查询
- **🌐 NCBI Taxonomy 集成**: 自动下载 NCBI Taxonomy 数据库，确保物种拉丁名准确可靠
- **📦 数据库自动查找**: `DatabaseManager` 支持自动查找 `organism/v{date}/{species}/` 目录结构，无需手动指定 `--database-dir`
- **🧬 GOA 文件名修复**: 修正 EBI GOA 文件名格式（`{taxid}.goa`），确保下载和构建流程正常
- **🔄 gene2go 全文件扫描**: 修复 `_check_gene2go_has_taxid()` 仅检查前 10000 行的 bug，支持人类等大数据集物种
- **🧬 Genome 背景基因集**: 新增 `get_genome_genes()` 方法，支持 `gene_info.gz` 全基因组基因作为背景集
- **🔗 数据合并去重优化**: gene2go 与 GOA 合并以 taxid 为唯一标识去重，优先使用 gene2go 数据
- **✅ 全量测试通过**: 628 个 pytest 测试全部通过，E2E 验证三种 background-mode 均正常

### 2026-05-26: GSEA/GSVA/ssGSEA 功能扩展

- **🧬 GSEA 富集分析**: 基于置换检验的标准GSEA实现，支持NES标准化和Leading Edge识别
- **📊 ssGSEA 单样本分析**: 单样本水平计算通路富集得分，适用于免疫浸润和多样本比较
- **🔄 GSVA 基因集变异分析**: 三种计算方法（Random Walk/PLAGE/Z-score），适用于样本聚类和异质性分析
- **🎨 发表级可视化**: 新增12+种可视化图表类型
  - GSEA: 富集曲线图、NES条形图、气泡图
  - ssGSEA/GSVA: 通路活性热图、组间比较图、样本相关性热图
  - 通用: 通路网络图、火山图、方法比较图、UpSet图
- **📁 GMT文件生成**: 数据库构建时自动生成GMT格式基因集文件，支持GSEA/ssGSEA/GSVA分析
- **✅ 全量端到端测试**: 新增63个E2E测试，覆盖500基因排序列表和6000×6样本表达矩阵

### 2026-05-15: Bug 修复

- **KEGG 层级修复**: 修复批量 API 查询 CLASS 字段解析逻辑，KEGG Term_Name 正确显示三层分类（如 `Cellular Processes|Cell Growth And Death|Cell Cycle`）
- **Term_Name 格式统一**: GO/KEGG/Reactome/DO 各数据库 Term_Name 统一为首字母大写层级格式
- **AI 解读优化**: 英文简洁分点格式（~250 词），HTML 正确渲染 Markdown（加粗/换行）
- **免责声明**: 每个 HTML 报告底部统一添加 AI 解读免责声明

---

*最后更新：2026-05-28（v2.2.0）*
