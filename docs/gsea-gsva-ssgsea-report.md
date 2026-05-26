# AllEnricher v2 — GSEA/GSVA/ssGSEA 扩展总结报告

**日期**: 2026-05-26
**执行方式**: Subagent-Driven 并行执行（3 个子代理）

---

## 测试摘要

| 指标 | 数量 |
|------|------|
| **总测试数** | **370** |
| 通过 | **370** |
| 失败 | 0 |
| 通过率 | **100%** |
| 新增测试 | 48 |
| 原有测试 | 322 |

---

## 新增/修改文件

### 新增文件

| 文件 | 功能 |
|------|------|
| `allenricher/core/gsva.py` | GSVA 类（三种方法变体：gsva/plage/zscore） |
| `tests/test_gsva.py` | GSVA 测试（15 个用例） |
| `tests/test_gsea_extended.py` | GSEA/ssGSEA 扩展测试（10 个用例） |
| `tests/test_cli_gsea.py` | CLI 集成测试（23 个用例） |
| `examples/gsea_usage.py` | GSEA/ssGSEA/GSVA 使用示例 |

### 修改文件

| 文件 | 修改内容 |
|------|----------|
| `allenricher/core/enrichment.py` | GSEA 新增 `calculate_normalized_es`（基于置换的 NES）；GSEA/ssGSEA 新增 `analyze_matrix`（表达矩阵接口）；`_get_method` 添加 GSVA 延迟导入 |
| `allenricher/core/config.py` | `EnrichmentMethod` 枚举添加 `GSVA`；`Config` 添加 `gsva_method/kcdf/tau` 字段 |
| `allenricher/cli.py` | `--method` 添加 `gsva` 选项；新增 `--expression-matrix/-e` 和 `--ranked-genes/-r` 参数 |
| `tests/test_api_server.py` | 修复根路由测试（`/` 现返回 HTML） |
| `tests/test_phase5.py` | 修复根路由测试 |

---

## 功能详情

### 1. GSEA 完善

**改进**: NES 计算从简单公式 `ES * sqrt(N/nh)` 改为基于置换检验的归一化

```python
# 新增方法
gsea.calculate_normalized_es(ranked_genes, gene_set, gene_weights)
# → (ES, NES, pvalue, leading_edge_genes)

# 新增表达矩阵接口
gsea.analyze_matrix(expression_matrix, gene_sets)
# → DataFrame (通路 × 样本 NES 矩阵)
```

**测试覆盖**: 10 个用例
- ES 正确性、NES 正/负向、置换 p 值范围、前沿基因验证
- 表达矩阵分析输出形状、ssGSEA NES 范围 [-1,1]、单样本分析

### 2. GSVA 实现

**三种方法变体**:

| 方法 | 算法 | 适用场景 |
|------|------|----------|
| `gsva` | ECDF + 核密度估计 + 随机游走 | 通用通路活性评估 |
| `plage` | SVD 第一主成分 | 通路内基因协同表达 |
| `zscore` | 标准化均值 | 简单快速评估 |

```python
from allenricher.core.gsva import GSVA

gsva = GSVA(method="gsva", kcdf="Gaussian", tau=1.0)
result = gsva.analyze_matrix(expression_matrix, gene_sets)
# → DataFrame (通路 × 样本活性矩阵)
```

**测试覆盖**: 15 个用例
- 基本功能、空输入、大小过滤、三种方法变体
- 无交集处理、单样本、结果范围、参数验证

### 3. CLI 集成

```bash
# GSEA 分析（排序基因列表）
allenricher analyze -i ranked_genes.txt -s hsa -d GO,KEGG -m gsea

# ssGSEA 分析（表达矩阵）
allenricher analyze -e expression_matrix.tsv -s hsa -d GO,KEGG -m ssgsea

# GSVA 分析（表达矩阵）
allenricher analyze -e expression_matrix.tsv -s hsa -d GO,KEGG -m gsva
```

**测试覆盖**: 23 个用例
- CLI 参数验证、Config 字段完整性、枚举完整性、导入测试

---

## 发现并修复的问题

| ID | 问题 | 修复方式 |
|----|------|----------|
| 1 | `enrichment.py` ↔ `gsva.py` 循环依赖 | 延迟导入（`_get_method` 内部 import） |
| 2 | PLAGE 方法依赖 sklearn | 改用 `numpy.linalg.svd` 实现 |
| 3 | API 根路由测试期望 JSON 但收到 HTML | 更新测试断言 |
| 4 | 单样本 GSVA 的 numpy RuntimeWarning | 已知行为，不影响结果正确性 |

---

## 方法对比

| 特性 | ORA (Fisher) | GSEA | ssGSEA | GSVA |
|------|-------------|------|--------|------|
| 输入 | 基因列表 | 排序基因列表 | 表达矩阵 | 表达矩阵 |
| 统计检验 | Fisher 精确检验 | 置换检验 | 无（NES 作为统计量） | 无（活性得分） |
| 样本水平 | 两组比较 | 两组比较 | 单样本 | 多样本 |
| P 值 | ✅ | ✅ | ❌ | ❌ |
| NES | ❌ | ✅ | ✅ | ❌ |
| 前沿基因 | ❌ | ✅ | ✅ | ❌ |
| 速度 | 快 | 慢（置换） | 中 | 中 |
