# KEGG 层级 & AI 解读问题排查与修复

**创建时间**: 2026-05-15
**状态**: ✅ 已修复并验证

## 问题列表

### 问题 1: KEGG 层级丢失 - 只有一层

**现象**: KEGG 应显示三层结构 `Cellular Processes|Cell growth and death|Cell Cycle`，但 160 个通路只显示一层（如 `Virion - Ebolavirus...`）

**根因**: `kegg_fetcher.py` 中批量 API 查询的 CLASS 字段解析逻辑错误

```python
# 错误代码：假设 ENTRY 下一行就是 CLASS
if line.startswith(f"ENTRY") and pw_id in line:
    if j + 1 < len(lines) and lines[j + 1].startswith("CLASS"):
```

**实际 API 返回结构**:
```
ENTRY       hsa03250                    Pathway
NAME        Viral life cycle - HIV-1
CLASS       Genetic Information Processing; Information processing in viruses
```

ENTRY → NAME → CLASS，不是 ENTRY → CLASS！

**修复**: 使用状态机正确解析，检测 ENTRY 开始、/// 结束，在中间查找 CLASS 字段

**修复效果**:
| 指标 | 修复前 | 修复后 |
|------|--------|--------|
| Uncategorized | 160 (43%) | 7 (2%) |
| 有分类 | 211 (57%) | 362 (98%) |
| 3层结构 | 部分 | 290 个 |

### 问题 2: AI 解读字数超限

**现象**: 约 300 词期望，实际输出 350-550 词

**修复**: Prompt 改为 `IMPORTANT: Keep total response under 250 words`，使用精确的分点结构

### 问题 3: HTML Markdown 格式未生效

**现象**: AI 解读在 HTML 中显示 Markdown 源码

**修复**:
- `\n` → `<br>` 换行
- `**text**` → `<strong>text</strong>` 加粗（正则替换）

### 问题 4: Disclaimer 重复

**修复**: 各模块末尾移除，只保留 HTML 底部一个全局 Disclaimer

### 问题 5: Term_Name 层级格式

**修复**: `manager.py` 新增 `_format_term_name()` + `_capitalize()` 方法
- GO: `Biological Process|Chromosome Segregation`
- KEGG: `Cellular Processes|Cell Growth And Death|Cell Cycle`
- 保留 DNA/RNA/ATP 等全大写词

## 修改文件

| 文件 | 修改内容 |
|------|----------|
| `allenricher/database/kegg_fetcher.py` | 修复 CLASS 字段解析逻辑（状态机）、新增 hardcoded 备用映射 |
| `allenricher/database/manager.py` | 新增 `_format_term_name()` + `_capitalize()` 方法 |
| `allenricher/ai/interpreter.py` | 统一所有后端 prompt 为英文简洁格式 |
| `allenricher/report/generator.py` | HTML Markdown 转 HTML、全局 Disclaimer |
| `README.md` | 文档更新 |
