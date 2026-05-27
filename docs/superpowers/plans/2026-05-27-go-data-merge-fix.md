# GO 数据合并修复计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 修复 GO gene2go 和 UniProt GOA 数据合并逻辑，确保使用 taxid 作为唯一标识去重，拉丁名统一从 NCBI Taxonomy 获取，优先使用 gene2go 数据

**Architecture:** 修改 downloader.py 中的注册表构建和合并方法，确保数据一致性和准确性

**Tech Stack:** Python, pytest, AllEnricher v2

---

## 问题背景

根据讨论结论，当前代码存在以下问题：

1. **拉丁名来源不一致**: gene2go 和 GOA 使用各自的拉丁名，未统一从 NCBI Taxonomy 获取
2. **去重逻辑不完善**: 合并时可能重复计数物种
3. **优先级不明确**: 未明确优先使用 gene2go 数据

## 修复要求

1. ✅ 统一使用 NCBI 十进制 taxid 作为物种唯一识别编码
2. ✅ 所有物种拉丁名统一从 NCBI Taxonomy 获取
3. ✅ 合并两个 GO 数据集时以 taxid 为唯一标识去重
4. ✅ 重复物种优先保留 gene2go 数据，仅当 gene2go 无数据时补充 GOA 数据

---

## Task 1: 验证当前代码状态

**Files:**
- Read: `allenricher/database/downloader.py`
- Read: `allenricher/database/species_registry.py`

- [ ] **Step 1: 检查当前 latin_name 来源**

```bash
cd f:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2
grep -n "latin_name" allenricher/database/downloader.py | head -20
```

- [ ] **Step 2: 检查当前合并逻辑**

```bash
grep -n "_merge_go_registries\|_build_go_registry\|_merge_all_registries" allenricher/database/downloader.py
```

- [ ] **Step 3: 验证当前注册表数据**

```bash
python -c "
from allenricher.database.species_registry import SpeciesRegistry
registry = SpeciesRegistry.load_default()
print(f'Total species: {len(registry.entries)}')

# 检查几个关键物种
for taxid in [9606, 10090, 10116]:
    entry = registry.query_by_taxid(taxid)
    if entry:
        print(f'{taxid}: {entry.latin_name}, GO={entry.has_go}, source={entry.go_source}')
"
```

---

## Task 2: 修复 gene2go 注册表构建逻辑

**Files:**
- Modify: `allenricher/database/downloader.py:565-657`

- [ ] **Step 1: 验证 `_build_go_registry` 使用 NCBI Taxonomy**

当前代码应该已经使用 NCBI Taxonomy，验证如下：

```python
# 在 _build_go_registry 方法中，检查第 603-614 行
taxonomy_dir = self.basic_dir / "taxonomy"
taxid_to_name = self._load_taxid_to_name_map(taxonomy_dir)

if not taxid_to_name:
    logger.info("NCBI Taxonomy 映射为空，尝试下载...")
    taxonomy_tsv = self._download_taxonomy_names(taxonomy_dir)
    if taxonomy_tsv:
        taxid_to_name = self._load_taxid_to_name_map(taxonomy_dir)
```

- [ ] **Step 2: 确保 gene2go 不使用 gene_info 的 Symbol 作为 latin_name**

检查第 615-634 行的回退逻辑，确保仅在 NCBI Taxonomy 完全不可用时才使用 gene_info：

```python
# 如果仍然为空，回退到 gene_info（仅用于无法获取 taxonomy 的情况）
if not taxid_to_name:
    logger.warning("无法获取 NCBI Taxonomy，使用 gene_info 作为备选")
    # ... 回退逻辑
```

- [ ] **Step 3: 运行测试验证 gene2go 注册表**

```bash
python -c "
from allenricher.database.downloader import DataDownloader
from pathlib import Path
import csv

d = DataDownloader()
go_dir = Path('./database/basic/go/GO20260527')

# 检查 gene2go 注册表中的拉丁名
with open(go_dir / 'go_species_registry.tsv', 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        if row['taxid'] in ['9606', '10090', '10116']:
            print(f\"gene2go: taxid={row['taxid']}, name={row['latin_name']}\")
"
```

Expected:
```
gene2go: taxid=9606, name=Homo sapiens
gene2go: taxid=10090, name=Mus musculus
gene2go: taxid=10116, name=Rattus norvegicus
```

---

## Task 3: 修复 GOA 索引构建逻辑

**Files:**
- Modify: `allenricher/database/downloader.py:659-727`

- [ ] **Step 1: 检查 GOA 索引的 latin_name 来源**

当前 GOA 索引从文件名解析拉丁名：

```python
# 第 704 行
latin_name = species_part.replace("_", " ")
```

这会导致问题，例如 `30675.S.goa` 解析为 `S`。

- [ ] **Step 2: 修改 GOA 索引使用 NCBI Taxonomy**

在 `_download_goa_index` 方法中，添加使用 NCBI Taxonomy 获取 latin_name 的逻辑：

```python
def _download_goa_index(self, output_dir: Path) -> Path:
    """从 UniProt GOA FTP 获取物种索引"""
    # ... 现有代码 ...
    
    # 解析 HTML 提取 .goa 文件链接
    goa_entries: List[Dict[str, str]] = []
    for match in re.finditer(r'href="([^"]+\.goa(?:\.gz)?)"', resp.text):
        filename = match.group(1)
        # ... 解析 taxid 和 species_part ...
        
        # 临时存储，稍后统一查询 NCBI Taxonomy
        goa_entries.append({
            "taxid": int(taxid_str),
            "species_part": species_part,  # 临时保存，用于回退
            "filename": filename,
        })
    
    # 从 NCBI Taxonomy 获取 latin_name
    taxonomy_dir = self.basic_dir / "taxonomy"
    taxid_to_name = self._load_taxid_to_name_map(taxonomy_dir)
    
    # 为每个 entry 添加 latin_name
    final_entries = []
    for entry in goa_entries:
        taxid = entry["taxid"]
        # 优先使用 NCBI Taxonomy
        latin_name = taxid_to_name.get(taxid, "")
        # 如果 NCBI Taxonomy 中没有，使用文件名中的名称
        if not latin_name:
            latin_name = entry["species_part"].replace("_", " ")
            logger.warning(f"Taxid {taxid} not in NCBI Taxonomy, using filename: {latin_name}")
        
        final_entries.append({
            "taxid": taxid,
            "latin_name": latin_name,
            "filename": entry["filename"],
        })
    
    # ... 写入 goa_species_index.tsv ...
```

- [ ] **Step 3: 运行测试验证 GOA 索引**

```bash
python -c "
import csv

# 检查 GOA 索引中的拉丁名
with open('database/basic/go/GO20260527/goa_species_index.tsv', 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        if row['taxid'] in ['9606', '10090', '10116']:
            print(f\"GOA: taxid={row['taxid']}, name={row['latin_name']}\")
"
```

Expected:
```
GOA: taxid=9606, name=Homo sapiens
GOA: taxid=10090, name=Mus musculus
GOA: taxid=10116, name=Rattus norvegicus
```

---

## Task 4: 修复合并逻辑 - 优先使用 gene2go

**Files:**
- Modify: `allenricher/database/downloader.py:779-815`

- [ ] **Step 1: 检查当前 `_merge_go_registries` 逻辑**

当前逻辑（已修复）：

```python
if in_g2g and in_goa:
    source = "both"
    # 优先使用 gene2go 的 latin_name（来自 NCBI Taxonomy）
    g2g_name = g2g_data[taxid]["latin_name"]
    goa_name = goa_data[taxid]["latin_name"]
    if g2g_name:
        latin_name = g2g_name
    else:
        latin_name = goa_name
    gene_count = g2g_data[taxid]["gene_count"]
    term_count = g2g_data[taxid]["term_count"]
```

- [ ] **Step 2: 验证合并逻辑正确性**

确保以下优先级：
1. 如果 taxid 只在 gene2go 中 → 使用 gene2go 数据
2. 如果 taxid 只在 GOA 中 → 使用 GOA 数据
3. 如果 taxid 在两者中 → 使用 gene2go 的 latin_name 和统计数据

- [ ] **Step 3: 运行测试验证合并结果**

```bash
python -c "
import csv

# 检查合并后的注册表
with open('database/basic/go/GO20260527/go_species_registry.tsv', 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')
    both_count = 0
    g2g_only_count = 0
    goa_only_count = 0
    
    for row in reader:
        source = row['source']
        if source == 'both':
            both_count += 1
        elif source == 'ncbi_gene2go':
            g2g_only_count += 1
        elif source == 'uniprot_goa':
            goa_only_count += 1
    
    print(f'Both: {both_count}')
    print(f'Gene2GO only: {g2g_only_count}')
    print(f'GOA only: {goa_only_count}')
    print(f'Total: {both_count + g2g_only_count + goa_only_count}')
"
```

---

## Task 5: 验证最终 supported_species.tsv

**Files:**
- Read: `database/supported_species.tsv`

- [ ] **Step 1: 检查最终注册表数据质量**

```bash
python -c "
import csv

with open('database/supported_species.tsv', 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')
    rows = list(reader)

print(f'Total species: {len(rows)}')

# 检查拉丁名质量
empty_names = [r for r in rows if not r.get('latin_name')]
short_names = [r for r in rows if r.get('latin_name') and len(r['latin_name']) <= 2]

print(f'Empty latin_name: {len(empty_names)}')
print(f'Short latin_name (<=2 chars): {len(short_names)}')

# 检查关键物种
for taxid in ['9606', '10090', '10116']:
    for r in rows:
        if r['taxid'] == taxid:
            print(f\"{taxid}: {r['latin_name']}, GO={r['has_go']}, source={r['go_source']}\")
            break
"
```

- [ ] **Step 2: 验证物种数量合理**

```bash
python -c "
# 验证物种数量在合理范围内
# Gene2GO: ~2342 物种
# GOA: ~29501 物种
# 合并后应该接近 31822 物种（交集约 29501，gene2go 独有约 2321）

from allenricher.database.species_registry import SpeciesRegistry
registry = SpeciesRegistry.load_default()
print(f'Total species in registry: {len(registry.entries)}')

# 统计 GO 支持情况
go_count = sum(1 for e in registry.entries.values() if e.has_go)
print(f'Species with GO: {go_count}')
"
```

---

## Task 6: 运行完整测试套件

**Files:**
- Test: `tests/`

- [ ] **Step 1: 运行所有测试**

```bash
cd f:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2
python -m pytest tests/ -v --tb=short 2>&1 | tail -20
```

- [ ] **Step 2: 验证关键测试通过**

重点关注以下测试：
- `tests/test_cli.py::TestCmdList::test_list_species`
- `tests/test_species_registry.py`
- `tests/test_enrichment.py`

- [ ] **Step 3: 运行 E2E 测试**

```bash
python -m allenricher analyze -i '../AllEnricher-v1/example/example.glist' -s hsa -d GO --background-mode annotated -o test_output/final_test
```

Expected: 成功运行，生成富集结果

---

## 执行选项

**Plan complete and saved to `docs/superpowers/plans/2026-05-27-go-data-merge-fix.md`. Two execution options:**

**1. Subagent-Driven (recommended)** - I dispatch a fresh subagent per task, review between tasks, fast iteration

**2. Inline Execution** - Execute tasks in this session using executing-plans, batch execution with checkpoints

**Which approach?**
