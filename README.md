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
- [测试指南](#测试指南)
- [输出文件说明](#输出文件说明)
- [配置文件](#配置文件)
- [常见问题](#常见问题)

---

## 简介

AllEnricher v2.0 是 AllEnricher v1.0 的 Python 重构版本，保留了 v1 的核心功能并增加了多项改进：

- **兼容 v1 数据库**：可直接使用 v1 构建的 GO、KEGG、Reactome、DO、DisGeNET 数据库
- **多种统计方法**：Fisher 精确检验、超几何检验、GSEA、ssGSEA
- **多重检验校正**：BH、BY、Bonferroni、Holm
- **丰富可视化**：柱状图、气泡图（复用 v1 R 脚本）
- **AI 智能解读**：支持 OpenAI、Claude、DeepSeek、GLM、MiniMax、Ollama
- **REST API**：FastAPI 驱动的 Web 服务
- **交互式报告**：学术风格的 HTML 报告

---

## 功能特性

### 支持的数据库

| 数据库 | 描述 | 物种支持 |
|---------|------|----------|
| GO | Gene Ontology | 15464 种（需构建本地数据库） |
| KEGG | KEGG Pathway | 所有 KEGG 物种 |
| Reactome | Reactome Pathway | 16 种模式生物 |
| DO | Disease Ontology | 人类 |
| DisGeNET | 疾病-基因关联 | 人类 |

### 统计分析方法

| 方法 | 描述 | 适用场景 |
|------|------|----------|
| Fisher | Fisher 精确检验 | 标准富集分析 |
| Hypergeometric | 超几何检验 | Fisher 的替代方案 |
| GSEA | 基因集富集分析 | 排序基因列表 |
| ssGSEA | 单样本 GSEA | 单样本分析 |

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
- R（用于可视化，条形图和气泡图）
- Rscript（R 命令行工具）

### 可选依赖

- FastAPI + Uvicorn（API 服务）
- OpenAI/Anthropic Python 包（AI 解读）
- pytest（运行测试）

### R 包要求

基础 R 安装即可，条形图和气泡图使用 R base 原生函数和 ggplot2，不依赖其他第三方 R 包。

---

## 安装方式

### 方式一：Docker（推荐，最简单）

```bash
# 克隆项目
git clone https://github.com/zd105/AllEnricher.git
cd AllEnricher/AllEnricher-v2

# 构建 Docker 镜像（包含 R + Python + 所有依赖）
docker build -t allenricher:v2 -f Dockerfile.test .

# 运行测试
docker run --rm -v "$(pwd):/workspace" -w /workspace allenricher:v2 bash -c "python3 -m pytest tests/ -v"
```

> 注意：Dockerfile.test 已包含 R、Python、pytest 等所有依赖。

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
# 1. 安装 R
# Ubuntu/Debian:
sudo apt-get install r-base r-base-dev

# macOS:
brew install r

# Windows: 从 https://cran.r-project.org/bin/windows/base/ 下载安装

# 2. 安装 Python 依赖
pip install numpy pandas scipy statsmodels pyyaml requests tqdm openpyxl

# 3. 安装 pytest（测试用）
pip install pytest pytest-cov

# 4. 安装项目
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
├── plots/                 # 可视化图表
│   ├── GO_barplot.pdf
│   ├── GO_bubble.pdf
│   ├── KEGG_barplot.pdf
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
                            fisher, hypergeometric, gsea, ssgsea
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

#### serve - 启动 API 服务

```bash
allenricher serve --port 8000 --host 0.0.0.0
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

### 3.3 REST API 使用

```bash
# 启动服务
allenricher serve --port 8000

# 提交分析任务
curl -X POST "http://localhost:8000/api/analyze" \
  -H "Content-Type: application/json" \
  -d '{
    "genes": ["BRCA1", "TP53", "EGFR", "MYC", "KRAS"],
    "species": "hsa",
    "databases": ["GO", "KEGG"],
    "method": "fisher",
    "pvalue_cutoff": 0.05,
    "qvalue_cutoff": 0.05,
    "database_dir": "/path/to/v1/database/organism/v20190612/hsa"
  }'

# 查询任务状态
curl "http://localhost:8000/api/status/{job_id}"

# 获取结果
curl "http://localhost:8000/api/results/{job_id}"
```

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
    bash -c "python3 -m pytest tests/ -v"
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

### 5.2 PDF 可视化

- **柱状图**：横轴为 -log10(Q-value)，纵轴为条目名称，按类别着色
- **气泡图**：横轴为 Rich Factor，纵轴为条目名称，点大小为基因数，颜色为 -log10(Q-value)

### 5.3 HTML 报告

学术风格的交互式报告，包含：

- 分析摘要统计
- 交互式数据表格（支持排序、搜索、下载）
- 可视化图表链接
- 可选的 AI 解读内容

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
  - "pdf"

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

### Q3: R 脚本执行失败

确保 R 和 Rscript 在 PATH 中：

```bash
# 检查 R 是否可用
R --version
Rscript --version

# 如果不可用，添加到 PATH 或安装 R
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
2. 检查 R 是否正确安装
3. 使用 Docker 环境排除环境问题：

```bash
docker build -t allenricher:v2 -f Dockerfile.test .
docker run --rm -v "$(pwd):/workspace" -w /workspace allenricher:v2 \
    bash -c "python3 -m pytest tests/test_enrichment.py -v"
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

## 最近更新 (v2.0.1)

### 2026-05-15: Bug 修复

- **KEGG 层级修复**: 修复批量 API 查询 CLASS 字段解析逻辑，KEGG Term_Name 正确显示三层分类（如 `Cellular Processes|Cell Growth And Death|Cell Cycle`）
- **Term_Name 格式统一**: GO/KEGG/Reactome/DO 各数据库 Term_Name 统一为首字母大写层级格式
- **AI 解读优化**: 英文简洁分点格式（~250 词），HTML 正确渲染 Markdown（加粗/换行）
- **免责声明**: 每个 HTML 报告底部统一添加 AI 解读免责声明

---

*最后更新：2026-05-15*
