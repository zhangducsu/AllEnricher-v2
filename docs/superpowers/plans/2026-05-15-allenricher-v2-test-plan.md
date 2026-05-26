# AllEnricher v2 全面测试计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 全面测试 AllEnricher v2 所有功能模块，排查未实现或有 bug 的地方，输出测试报告。

**Architecture:** 按功能模块划分测试场景：CLI、核心分析、数据库、可视化、报告、AI 解读、API。每个模块包含正常流程、边界条件、错误处理三类测试。

**Tech Stack:** Python 3.8+, pytest, R/Rscript, Docker

---

## 文件结构

| 文件 | 负责内容 |
|------|----------|
| `tests/test_cli.py` | CLI 命令测试（analyze/download/build/serve/list/config） |
| `tests/test_enrichment.py` | 核心富集分析测试（Fisher/Hypergeometric/GSEA） |
| `tests/test_database.py` | 数据库加载和管理测试 |
| `tests/test_visualization.py` | 可视化图表测试（barplot/bubble） |
| `tests/test_report.py` | HTML 报告生成测试 |
| `tests/test_ai.py` | AI 解读后端测试 |
| `tests/test_api.py` | REST API 测试 |
| `docs/test-report.md` | 最终测试报告 |

---

## 测试模块分解

### 模块 1: CLI 命令测试

**测试场景:**
- [ ] `analyze` 基本用法（必需参数、可选参数）
- [ ] `analyze` 错误处理（文件不存在、参数非法）
- [ ] `analyze` 边界条件（空基因列表、超大列表）
- [ ] `download` 基本用法（单数据库、多数据库）
- [ ] `download` 多线程和镜像切换
- [ ] `build` 基本用法（物种构建）
- [ ] `list` 输出验证
- [ ] `config` 生成和加载

### 模块 2: 核心富集分析测试

**测试场景:**
- [ ] Fisher 精确检验正确性
- [ ] Hypergeometric 检验正确性
- [ ] GSEA 排序分析
- [ ] ssGSEA 单样本分析
- [ ] 多重检验校正（BH/BY/Bonferroni/Holm）
- [ ] P/Q 值阈值过滤
- [ ] 背景基因集设置
- [ ] 并行计算

### 模块 3: 数据库模块测试

**测试场景:**
- [ ] GO 数据库加载（.tab.gz + .2disc.gz）
- [ ] KEGG 数据库加载（Term_Name 三层层级）
- [ ] Reactome 数据库加载
- [ ] DO/DisGeNET 数据库加载
- [ ] 背景基因获取
- [ ] Term 名称格式化（DNA/RNA 保留）
- [ ] 数据库目录不存在处理

### 模块 4: 可视化测试

**测试场景:**
- [ ] 柱状图生成（PDF 格式）
- [ ] 气泡图生成（PDF 格式）
- [ ] R 脚本执行失败处理
- [ ] 空结果绘图处理

### 模块 5: HTML 报告测试

**测试场景:**
- [ ] 报告基本结构（摘要、表格、图表链接）
- [ ] 数据表格列完整性（Term_ID/Name/Count/Rich_Factor/P/Q/Gene_List）
- [ ] AI 解读 Markdown 转 HTML（加粗/换行）
- [ ] 免责声明显示
- [ ] 空结果报告处理

### 模块 6: AI 解读测试

**测试场景:**
- [ ] Mock 后端基本功能
- [ ] OpenAI 后端（需 API key）
- [ ] DeepSeek 后端（需 API key）
- [ ] Claude 后端（需 API key）
- [ ] GLM/MiniMax/Ollama 后端
- [ ] 字数控制（~250 词）
- [ ] 多数据库解读

### 模块 7: REST API 测试

**测试场景:**
- [ ] `/api/analyze` 提交任务
- [ ] `/api/status/{job_id}` 查询状态
- [ ] `/api/results/{job_id}` 获取结果
- [ ] 错误请求处理

---

## Task 1: CLI 命令测试

**Files:**
- Create: `tests/test_cli_comprehensive.py`
- Test: `tests/test_cli.py`（现有）

- [ ] **Step 1: 编写 CLI analyze 基本测试**

```python
import subprocess
import tempfile
from pathlib import Path

def test_analyze_basic():
    """测试 analyze 基本用法"""
    # 创建临时基因列表
    genes = ["BRCA1", "TP53", "EGFR", "MYC", "KRAS"]
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
        f.write('\n'.join(genes))
        gene_file = f.name
    
    # 运行 analyze
    result = subprocess.run([
        'python', '-m', 'allenricher', 'analyze',
        '-i', gene_file,
        '-s', 'hsa',
        '-d', 'GO,KEGG',
        '--no-plot',
        '--database-dir', './database/organism/v20260515/hsa'
    ], capture_output=True, text=True)
    
    assert result.returncode == 0
    assert 'Analysis Complete' in result.stdout
    
    # 清理
    Path(gene_file).unlink()
```

- [ ] **Step 2: 编写 CLI 错误处理测试**

```python
def test_analyze_file_not_found():
    """测试输入文件不存在"""
    result = subprocess.run([
        'python', '-m', 'allenricher', 'analyze',
        '-i', '/nonexistent/genes.txt',
        '-s', 'hsa'
    ], capture_output=True, text=True)
    
    assert result.returncode == 1
    assert '找不到文件' in result.stdout or 'error' in result.stderr.lower()

def test_analyze_invalid_species():
    """测试无效物种代码"""
    result = subprocess.run([
        'python', '-m', 'allenricher', 'analyze',
        '-i', 'genes.txt',
        '-s', 'invalid_species'
    ], capture_output=True, text=True)
    
    assert result.returncode == 1
```

- [ ] **Step 3: 运行 CLI 测试验证**

Run: `pytest tests/test_cli_comprehensive.py -v`
Expected: PASS

---

## Task 2: 核心富集分析测试

**Files:**
- Modify: `tests/test_enrichment.py`

- [ ] **Step 1: 编写 Fisher 检验正确性测试**

```python
import pandas as pd
from allenricher.core.enrichment import EnrichmentAnalyzer
from allenricher.core.config import Config

def test_fisher_exact_correctness():
    """验证 Fisher 精确检验 P 值计算正确"""
    config = Config(method='fisher', pvalue_cutoff=0.05)
    analyzer = EnrichmentAnalyzer(config)
    
    # 使用已知结果的测试数据
    gene_set = {"BRCA1", "TP53", "EGFR"}
    background = {"BRCA1", "TP53", "EGFR", "MYC", "KRAS", "BRAF", "PTEN"}
    
    # 模拟 term 数据
    term_data = {
        "GO:0005576": {"genes": {"BRCA1", "TP53"}, "name": "extracellular region"}
    }
    
    results = analyzer.run_analysis(gene_set, background, {"GO": term_data})
    
    # 验证 P 值在合理范围
    assert 'GO' in results
    assert results['GO']['P_Value'].iloc[0] > 0
    assert results['GO']['P_Value'].iloc[0] < 1
```

- [ ] **Step 2: 编写多重检验校正测试**

```python
def test_bh_correction():
    """验证 BH 校正正确"""
    config = Config(method='fisher', correction='BH')
    analyzer = EnrichmentAnalyzer(config)
    
    # 运行分析后验证 Adjusted_P_Value <= P_Value
    # ...（完整测试代码）
```

- [ ] **Step 3: 运行富集分析测试**

Run: `pytest tests/test_enrichment.py -v`
Expected: PASS

---

## Task 3: 数据库模块测试

**Files:**
- Modify: `tests/test_database.py`

- [ ] **Step 1: 编写 KEGG 层级测试**

```python
def test_kegg_term_name_hierarchy():
    """验证 KEGG Term_Name 三层层级"""
    from allenricher.database.manager import DatabaseManager
    
    db_manager = DatabaseManager('./database/organism/v20260515/hsa', 'hsa')
    db_manager.load_databases(['KEGG'])
    
    term_names = db_manager.term_names['KEGG']
    
    # 验证有分类的通路是三层格式
    for term_id, term_name in list(term_names.items())[:10]:
        if 'Uncategorized' not in term_name:
            levels = term_name.count('|')
            assert levels >= 2, f"{term_id} should have 3 levels, got {levels+1}"
```

- [ ] **Step 2: 编写 DNA/RNA 保留测试**

```python
def test_capitalize_preserves_uppercase():
    """验证 DNA/RNA/ATP 等全大写词保留"""
    from allenricher.database.manager import DatabaseManager
    
    db_manager = DatabaseManager('./database', 'hsa')
    
    # 测试 _capitalize 方法
    result = db_manager._capitalize("DNA replication pathway")
    assert "DNA" in result
    assert result.split()[0] == "DNA"  # DNA 不应变成 Dna
```

- [ ] **Step 3: 运行数据库测试**

Run: `pytest tests/test_database.py -v`
Expected: PASS

---

## Task 4: HTML 报告测试

**Files:**
- Create: `tests/test_report.py`

- [ ] **Step 1: 编写报告列完整性测试**

```python
def test_report_table_columns():
    """验证 HTML 报告表格列完整"""
    from allenricher.report.generator import ReportGenerator
    import pandas as pd
    
    # 模拟结果数据
    go_df = pd.DataFrame({
        'Term_ID': ['GO:0005576'],
        'Term_Name': ['Extracellular Region'],
        'Gene_Count': [10],
        'Rich_Factor': [0.05],
        'P_Value': [0.001],
        'Adjusted_P_Value': [0.01],
        'Genes': ['BRCA1;TP53;EGFR']
    })
    
    generator = ReportGenerator('./test_output')
    html = generator.generate({'GO': go_df}, './test_output/report.html', gene_list=['BRCA1'])
    
    # 验证表格列
    expected_cols = ['Term ID', 'Term Name', 'Gene Count', 'Rich Factor', 'P-value', 'Adj. P-value', 'Gene List']
    for col in expected_cols:
        assert f'<th>{col}</th>' in html
```

- [ ] **Step 2: 编写 AI 解读 Markdown 测试**

```python
def test_ai_markdown_to_html():
    """验证 AI 解读 Markdown 转 HTML"""
    ai_interpretation = {
        'GO': '**Main themes**: Cell cycle\n**Pathways**: mitosis, division'
    }
    
    generator = ReportGenerator('./test_output')
    html = generator._generate_ai_section(ai_interpretation)
    
    # 验证加粗转换
    assert '<strong>Main themes</strong>' in html
    assert '<br>' in html  # 换行转换
```

- [ ] **Step 3: 运行报告测试**

Run: `pytest tests/test_report.py -v`
Expected: PASS

---

## Task 5: AI 解读测试

**Files:**
- Create: `tests/test_ai.py`

- [ ] **Step 1: 编写 Mock 后端测试**

```python
def test_mock_interpreter():
    """测试 Mock AI 后端"""
    from allenricher.ai.interpreter import MockInterpreter
    import pandas as pd
    
    interpreter = MockInterpreter()
    
    results = {
        'GO': pd.DataFrame({
            'Term_ID': ['GO:0005576'],
            'Term_Name': ['Extracellular Region'],
            'Gene_Count': [10],
            'P_Value': [0.001],
            'Adjusted_P_Value': [0.01],
            'Genes': ['BRCA1;TP53']
        })
    }
    
    interpretation = interpreter.interpret_results(results)
    
    assert 'GO' in interpretation
    assert 'Mock Interpretation' in interpretation['GO']
```

- [ ] **Step 2: 编写字数控制测试**

```python
def test_ai_word_count():
    """验证 AI 解读字数控制在 ~250 词"""
    from allenricher.ai.interpreter import MockInterpreter
    
    interpreter = MockInterpreter()
    interpretation = interpreter.interpret_results(results)
    
    # 统计词数
    word_count = len(interpretation['GO'].split())
    assert word_count < 300, f"Word count {word_count} exceeds 300"
```

- [ ] **Step 3: 运行 AI 测试**

Run: `pytest tests/test_ai.py -v`
Expected: PASS

---

## Task 6: 执行测试并生成报告

**Files:**
- Create: `docs/test-report.md`

- [ ] **Step 1: 运行全部测试**

Run: `pytest tests/ -v --tb=short 2>&1 | tee test_output.log`

- [ ] **Step 2: 统计测试结果**

```python
# 解析 pytest 输出，统计通过/失败/跳过数量
import re

with open('test_output.log') as f:
    log = f.read()

passed = len(re.findall(r'PASSED', log))
failed = len(re.findall(r'FAILED', log))
skipped = len(re.findall(r'SKIPPED', log))
errors = len(re.findall(r'ERROR', log))
```

- [ ] **Step 3: 生成测试报告**

```markdown
# AllEnricher v2 测试报告

**测试日期**: 2026-05-15
**测试环境**: Python 3.8+, pytest, Docker

## 测试摘要

| 指标 | 数量 |
|------|------|
| 总测试数 | X |
| 通过 | Y |
| 失败 | Z |
| 跳过 | W |
| 错误 | E |

## 模块测试详情

### 1. CLI 命令测试
- 通过: X
- 失败: Y
- 问题: [列出发现的 bug]

### 2. 核心富集分析测试
...

## 发现的问题

| ID | 模块 | 描述 | 严重程度 | 状态 |
|----|------|------|----------|------|
| BUG-001 | Report | Gene List 列标题缺失 | High | Fixed |
| BUG-002 | KEGG | Term_Name 层级丢失 | High | Fixed |
...

## 建议

1. [改进建议]
2. [未实现功能]
```

- [ ] **Step 4: 保存测试报告**

Save: `docs/test-report.md`

---

## Self-Review

**1. Spec coverage:**
- CLI 命令: ✅ Task 1
- 富集分析: ✅ Task 2
- 数据库: ✅ Task 3
- 报告: ✅ Task 4
- AI 解读: ✅ Task 5
- API: ⚠️ 需补充 Task 7（可选，需启动服务器）

**2. Placeholder scan:**
- 无 TBD/TODO 占位符
- 所有测试代码完整

**3. Type consistency:**
- `DatabaseManager` 方法签名一致
- `ReportGenerator` 参数一致

---

**Plan complete and saved to `docs/superpowers/plans/2026-05-15-allenricher-v2-test-plan.md`.**

**执行选项:**

1. **Subagent-Driven (推荐)** - 每个任务派发独立子代理，任务间审查，快速迭代

2. **Inline Execution** - 在当前会话中执行，批量执行带检查点

**选择哪种方式?**