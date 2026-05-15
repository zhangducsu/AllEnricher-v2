# KEGG Term_Name 层级异常问题排查

## 问题描述
KEGG 富集分析结果的 TSV 文件中，Term_Name 层级显示异常。需要从根本上排查数据流，不使用硬编码修复。

## 数据流分析

```
KEGG REST API (get/pathway/hsa{id})
    ↓ 解析 CLASS 字段
_fetch_pathway_categories() → pathway_categories dict
    ↓ 保存
pathway_summary.txt (格式: hsa00010\tCategory|Subcategory|PathwayName)
    ↓ KEGGParser 读取
build_database() → kegg2disc.gz
    ↓ manager.py 加载
_load_term_names() → term_names["KEGG"]
    ↓ 格式化
_format_term_name() → Term_Name 列
```

## 排查步骤

### 步骤1: 检查 API 返回的真实数据
- 测试 KEGG API `get/hsa04110` 返回的 CLASS 字段
- 确认是否包含三层分类信息

### 步骤2: 检查 pathway_summary.txt 内容
- 查看实际生成的 pathway_summary.txt 格式
- 确认分类信息是否正确保存

### 步骤3: 检查 kegg2disc.gz 生成
- 查看 KEGGParser 如何读取 pathway_summary.txt
- 确认写入 kegg2disc.gz 的格式

### 步骤4: 检查 manager.py 加载
- 确认 _load_term_names() 是否正确读取
- 确认 _format_term_name() 是否正确处理

### 步骤5: 检查 enrichment 结果
- 查看最终 TSV 中的 Term_Name 列

## 当前状态
- 之前使用了硬编码映射 `_get_hardcoded_categories()` 作为备用
- 需要找到为什么 API 获取失败或数据丢失的真正原因
