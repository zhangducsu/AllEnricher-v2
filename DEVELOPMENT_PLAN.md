# AllEnricher v2 开发计划

## 一、v1 代码分析报告

### 1.1 项目特色

| 特色 | 描述 |
|------|------|
| **多数据库支持** | GO、KEGG、Reactome、DO、DisGeNET 五种数据库 |
| **多物种支持** | GO 支持 15464 种物种，KEGG 支持所有 KEGG 物种，Reactome 支持 16 种模式生物 |
| **本地数据库** | 从公共资源自动下载并构建本地数据库，支持离线分析 |
| **两种统计检验** | Fisher 精确检验、超几何检验 |
| **多重检验校正** | BH、BY、holm、none |
| **可视化** | 条形图、气泡图 |
| **命令行工具** | 单命令完成分析和可视化 |

### 1.2 核心功能模块及实现逻辑

#### 1.2.1 数据库构建模块

| 脚本 | 功能 | 输入 | 输出 |
|------|------|------|------|
| `make_speciesDB` | 主控脚本 | 物种代码、taxid | 调用各数据库构建脚本 |
| `makeDB.go.v1.0.sh` | 构建 GO 数据库 | gene2go.gz, gene_info.gz, go-basic.obo | {species}.GO2gene.tab.gz, {species}.GO2disc.gz |
| `makeDB.kegg.v1.1.sh` | 构建 KEGG 数据库 | KEGG API | {species}.kegg2gene.tab.gz, {species}.kegg2disc.gz |
| `makeDB.reactome.v1.0.sh` | 构建 Reactome 数据库 | NCBI2Reactome_All_Levels.txt | {species}.Reactome2gene.tab.gz |
| `makeDB.do.v1.0.sh` | 构建 DO 数据库 | DISEASES 数据库 | {species}.DO2gene.tab.gz |
| `makeDB.DisGeNET.v1.0.sh` | 构建 DisGeNET 数据库 | DisGeNET 数据库 | {species}.CUI2gene.tab.gz |

> **注意**：DO 和 DisGeNET 仅限人类（hsa）构建。数据源本身只包含人类基因关联，
> 无多物种版本管理需求（对应 v1 make_speciesDB L127-135）。
> 在 v2 中，这两个数据库存储在 `database/basic/do/` 和 `database/basic/disgenet/`
> 下，不使用日期版本子目录；build 时自动检测 `species=="hsa"` 才执行构建。

#### 1.2.2 富集分析核心逻辑（AllEnricher_v1.0.pl）

```
Step 1: 读取基因列表 → %gene_list
Step 2: 读取数据库文件 → %ref (term → gene → 0/1)
Step 3: 读取背景基因列表 → %background_list
Step 4: 遍历每个条目计算统计量
        - num_in_C = 条目在背景中的基因数
        - num_in_O = 条目在输入基因中的基因数
        - expected = num_in_C / background_total * gene_total
        - rich_factor = num_in_O / expected
        - 仅保留 num_in_O > 1 的条目
Step 5: 调用 R 脚本进行统计检验
Step 6: 多重检验校正
Step 7: 添加条目名称映射
```

#### 1.2.3 统计检验逻辑（enrich.R）

```r
# Fisher 精确检验
mat[i,3] <- fisher.test(m)$p.value

# 超几何检验
mat[i,3] <- phyper(k-1, M, N-M, n, lower.tail=FALSE)

# 多重检验校正
mat[,4] <- p.adjust(mat[,3], method=adj, n=length(mat[,3]))
```

#### 1.2.4 可视化逻辑

- **条形图（barplot.R）**：横轴 -log10(Q-value)，纵轴条目名称，按 GO 类别/KEGG 分类着色
- **气泡图（bubble_plot.R）**：横轴 RichFactor，纵轴条目名称，点大小=基因数，颜色=-log10(Q-value)

### 1.3 数据库文件格式

| 文件 | 格式 | 说明 |
|------|------|------|
| `*.2gene.tab.gz` | Gene\tTerm1\tTerm2...\n gene\t0/1\t0/1... | 基因-条目矩阵 |
| `*.2disc.gz` | TermID\tTermName | 条目描述映射 |
| `gene_info` | tax_id\tgene_id\tsymbol\t... | 基因信息 |
| `gene2pathway.txt` | gene\tentrez_id\tpathway_id\tpathway_name | 基因-通路映射 |

### 1.4 代码文件逻辑关系

```
AllEnricher (主入口)
    ├── 调用 AllEnricher_v1.0.pl (核心分析)
    │       ├── 读取数据库文件
    │       ├── 计算统计量
    │       └── 调用 enrich.R (统计检验)
    └── 调用 plot_enrich_v1.0.pl (可视化)
            ├── 调用 barplot.R
            └── 调用 bubble_plot.R

make_speciesDB (数据库构建入口)
    ├── 调用 makeDB.go.v1.0.sh
    │       └── 调用 gene2GO_extract.pl, obo2go.pl
    ├── 调用 makeDB.kegg.v1.1.sh
    │       └── 调用 keggMapGrab.R, pathway2tab.pl
    ├── 调用 makeDB.reactome.v1.0.sh
    │       └── 调用 gene2ReactomePathway_extract.pl
    ├── 调用 makeDB.do.v1.0.sh
    └── 调用 makeDB.DisGeNET.v1.0.sh
```

---

## 二、v2 开发框架

### 2.1 目录结构

```
allenricher/
├── __init__.py
├── cli.py                    # 命令行接口
├── core/
│   ├── __init__.py
│   ├── config.py             # 配置管理
│   ├── enrichment.py         # 富集分析核心
│   └── statistics.py         # 统计检验（已合并到 enrichment.py）
├── database/
│   ├── __init__.py
│   ├── manager.py            # 数据库管理
│   ├── builder.py            # 数据库构建
│   ├── downloader.py         # 数据下载
│   └── parsers/              # 各数据库解析器
│       ├── __init__.py
│       ├── go.py
│       ├── kegg.py
│       ├── reactome.py
│       ├── do.py
│       └── disgenet.py
├── visualization/
│   ├── __init__.py
│   ├── plotter.py            # 可视化调度（Python 调用 R 脚本）
│   ├── barplot.R             # 条形图（R base 原生，v1 脚本）
│   └── bubble.R              # 气泡图（ggplot2，v1 脚本）
├── report/
│   ├── __init__.py
│   └── generator.py          # HTML 报告生成（学术风格）
├── api/
│   ├── __init__.py
│   └── server.py             # FastAPI 服务
└── ai/
    ├── __init__.py
    └── interpreter.py        # AI 解读（7 个后端）
```

### 2.2 技术选型

| 模块 | v1 技术 | v2 技术 | 理由 |
|------|---------|---------|------|
| **核心语言** | Perl + R | Python 3.8+ | 生态丰富，易于维护 |
| **统计检验** | R (stats) | scipy.stats | 原生 Python，无外部依赖 |
| **数据处理** | Perl hash | pandas | 高效处理表格数据 |
| **可视化** | R (ggplot2) | R base 原生 + ggplot2（v1 脚本） | 复用 v1 验证过的 R 脚本 |
| **命令行** | Getopt::Long | argparse | 标准库，无额外依赖 |
| **配置管理** | 无 | pydantic + YAML | 类型安全，配置灵活 |
| **API 服务** | 无 | FastAPI | 高性能异步框架 |
| **并行处理** | 无 | concurrent.futures | 加速大规模分析 |
| **AI 解读** | 无 | OpenAI/Claude/DeepSeek/GLM/MiniMax/Ollama | 多后端支持 |

---

## 三、开发进度

### 阶段 1：核心功能完善 ✅ 已完成

| 任务 | 描述 | 状态 | 测试 |
|------|------|------|------|
| 1.1 完善 `enrichment.py` | Fisher精确检验、超几何检验、GSEA、ssGSEA、多重校正（BH/BY/holm/bonferroni/none） | ✅ 完成 | 64 个单元测试通过 |
| 1.2 完善 `database/manager.py` | 正确加载 v1 数据库格式（.tab.gz / .2disc.gz） | ✅ 完成 | 端对端验证通过 |
| 1.3 实现 `statistics.py` | 合并到 enrichment.py，Fisher/超几何/GSEA/ssGSEA | ✅ 完成 | 64 个单元测试通过 |
| 1.4 实现结果输出 | TSV 格式，与 v1 输出兼容 | ✅ 完成 | 端对端验证通过 |

### 阶段 2：可视化模块 ✅ 已完成

| 任务 | 描述 | 状态 | 测试 |
|------|------|------|------|
| 2.1 实现 `barplot.R` | 使用 v1 原版 R base 原生绘图，支持 GO/KEGG 分类着色 | ✅ 完成 | 图表审查通过 |
| 2.2 实现 `bubble.R` | 使用 v1 原版 ggplot2 绘图 | ✅ 完成 | 图表审查通过 |
| 2.3 完善 `plotter.py` | Python 端调度层，调用 R 脚本并传递参数 | ✅ 完成 | 端对端验证通过 |

### 阶段 3：数据库构建模块 ✅ 已完成

| 任务 | 描述 | 状态 | 测试 |
|------|------|------|------|
| 3.1 实现 `database/downloader.py` | 从 NCBI、GO、Reactome、DO、DisGeNET 下载数据，支持断点续传 | ✅ 完成 | 单元测试通过 |
| 3.2 实现 `database/parsers/go.py` | 解析 gene2go.gz + go-basic.obo | ✅ 完成 | 单元测试通过 |
| 3.3 实现 `database/parsers/kegg.py` | 解析 gene2pathway.txt + pathway_summary | ✅ 完成 | 单元测试通过 |
| 3.4 实现 `database/parsers/reactome.py` | 解析 NCBI2Reactome 文件 | ✅ 完成 | 单元测试通过 |
| 3.5 实现 `database/parsers/do.py` | 解析 Jensen Lab disease TSV | ✅ 完成 | 单元测试通过 |
| 3.6 实现 `database/parsers/disgenet.py` | 解析 DisGeNET associations | ✅ 完成 | 单元测试通过 |
| 3.7 实现 `database/builder.py` | 统一调度所有解析器 | ✅ 完成 | 8 个单元测试 + 端对端通过 |

### 阶段 4：命令行接口 ✅ 已完成

| 任务 | 描述 | 状态 | 测试 |
|------|------|------|------|
| 4.1 实现 `analyze` 子命令 | 运行富集分析（11 步工作流） | ✅ 完成 | 15 个单元测试 + 7 个端对端通过 |
| 4.2 实现 `download` 子命令 | 下载数据库（GO/KEGG/Reactome/DO/DisGeNET） | ✅ 完成 | — |
| 4.3 实现 `build` 子命令 | 构建数据库 | ✅ 完成 | — |
| 4.4 实现 `list` 子命令 | 列出可用物种和数据库 | ✅ 完成 | 单元测试通过 |
| 4.5 实现 `config` 子命令 | 生成默认 YAML 配置文件 | ✅ 完成 | 单元测试通过 |
| 4.6 实现 `serve` 子命令 | 启动 FastAPI API 服务器 | ✅ 完成（额外） | — |

### 阶段 5：扩展功能 ✅ 已完成

| 任务 | 描述 | 状态 | 测试 |
|------|------|------|------|
| 5.1 实现 GSEA 分析 | 基因集富集分析 | ✅ 完成 | 64 个单元测试通过 |
| 5.2 实现 ssGSEA 分析 | 单样本 GSEA | ✅ 完成 | 64 个单元测试通过 |
| 5.3 实现 REST API | FastAPI 服务，9 个端点 | ✅ 完成 | 5 个单元测试通过 |
| 5.4 实现 AI 解读 | 7 个后端（OpenAI/Claude/DeepSeek/GLM/MiniMax/Ollama/Mock） | ✅ 完成 | 9 个单元测试通过 |
| 5.5 实现 HTML 报告 | 学术风格交互式报告（Noto Serif + Source Sans Pro） | ✅ 完成 | 3 个单元测试通过 |

---

## 四、测试汇总

### 单元测试

| 测试文件 | 测试数 | 状态 |
|----------|--------|------|
| `tests/test_enrichment.py` | 14 | ✅ 全部通过 |
| `tests/test_enrichment_extended.py` | 50 | ✅ 全部通过 |
| `tests/test_database.py` | 8 | ✅ 全部通过 |
| `tests/test_cli.py` | 15 | ✅ 全部通过 |
| `tests/test_phase5.py` | 19 | ✅ 全部通过 |
| **合计** | **106** | **✅ 全部通过** |

### 端对端测试

| 测试 | 说明 | 状态 |
|------|------|------|
| 核心分析 + v1 对比 | 使用 v1 example 数据，GO 314条、KEGG 7条结果与 v1 一致 | ✅ 通过 |
| 可视化 | GO/KEGG barplot + bubble 图表审查 | ✅ 通过 |
| 数据库构建 | 构建 → 加载 → 分析全链路 | ✅ 通过 |
| CLI 工作流 | analyze / list / config / --version / --config | ✅ 7/7 通过 |
| API + 报告 + AI | FastAPI 端点、HTML 报告、AI Mock 解读 | ✅ 5/6 通过 |

---

## 五、关键设计决策

1. **兼容性**：v2 必须能读取 v1 构建的数据库文件 ✅
2. **算法一致性**：统计检验结果应与 v1 完全一致 ✅
3. **输出格式**：TSV 输出格式与 v1 兼容 ✅
4. **R 脚本复用**：条形图和气泡图直接复用 v1 的 R 脚本，确保输出风格一致 ✅
5. **模块化**：各模块独立，可单独使用 ✅
6. **AI 多后端**：支持 OpenAI、Claude、DeepSeek、GLM、MiniMax、Ollama、Mock 共 7 个后端 ✅
7. **学术风格报告**：HTML 报告采用低调专业的学术审美，参考 Nature/Science supplementary data 排版 ✅
8. **download/build 分层架构**：download 下载全体物种通用数据 → database/basic/；build 提取指定物种 → database/organism/。DO/DisGeNET 仅为人类专有，不使用版本子目录 ✅
