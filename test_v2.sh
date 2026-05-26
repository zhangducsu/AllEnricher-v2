#!/bin/bash
# AllEnricher v2 Docker 测试脚本

set -e

# 检查 Docker 是否运行
if ! docker info > /dev/null 2>&1; then
    echo "错误: Docker 未运行，请先启动 Docker Desktop"
    exit 1
fi

# 工作目录
WORKDIR="$(cd "$(dirname "$0")" && pwd)"

echo "=========================================="
echo "AllEnricher v2 Docker 测试"
echo "使用 v1 的示例数据进行测试"
echo "=========================================="

# 检查文件
if [ ! -f "$WORKDIR/../AllEnricher-v1/example/example.glist" ]; then
    echo "错误: 基因列表文件不存在"
    exit 1
fi

if [ ! -d "$WORKDIR/../AllEnricher-v1/database/organism/v20190612/hsa" ]; then
    echo "错误: 数据库目录不存在"
    exit 1
fi

# 创建结果目录
mkdir -p "$WORKDIR/results"

# 运行测试容器
docker run --rm \
    -v "$WORKDIR:/workspace/allenricher" \
    -v "$WORKDIR/../AllEnricher-v1/example:/workspace/example" \
    -v "$WORKDIR/../AllEnricher-v1/database:/workspace/database" \
    -v "$WORKDIR/results:/workspace/results" \
    -w /workspace \
    my-bio-env:latest \
    python3 << 'PYTHON_SCRIPT'
import sys
import os

# 添加 v2 到路径
sys.path.insert(0, '/workspace/allenricher')

from pathlib import Path
from allenricher import EnrichmentAnalyzer, Config
from allenricher.database.manager import DatabaseManager

print("=" * 60)
print("AllEnricher v2 测试")
print("=" * 60)

# 配置
gene_list = "/workspace/example/example.glist"
database_dir = "/workspace/database/organism/v20190612/hsa"
output_dir = "/workspace/results"
species = "hsa"

# 检查文件
if not os.path.exists(gene_list):
    print(f"错误: 基因列表文件不存在: {gene_list}")
    sys.exit(1)

if not os.path.exists(database_dir):
    print(f"错误: 数据库目录不存在: {database_dir}")
    sys.exit(1)

# 创建配置
config = Config(
    species=species,
    databases=["GO", "KEGG"],
    method="fisher",
    pvalue_cutoff=0.05,
    qvalue_cutoff=0.05,
    min_genes=2
)

# 创建输出目录
Path(output_dir).mkdir(parents=True, exist_ok=True)

# 加载基因列表
print(f"\n[1] 加载基因列表: {gene_list}")
analyzer = EnrichmentAnalyzer(config)
gene_set = analyzer.load_gene_list(gene_list)
print(f"    基因数量: {len(gene_set)}")

# 加载数据库
print(f"\n[2] 加载数据库: {database_dir}")
db_manager = DatabaseManager(database_dir, species)
db_manager.load_databases(config.databases)
background_set = db_manager.get_background_genes()
print(f"    背景基因数量: {len(background_set)}")

# 获取数据库数据
database_data = db_manager.get_all_term_data()
print(f"    已加载数据库: {list(database_data.keys())}")

# 运行分析
print("\n[3] 运行富集分析...")
results = analyzer.run_analysis(gene_set, background_set, database_data)

# 检查结果
if not results or len(results) == 0:
    print("    未找到显著富集结果")
    sys.exit(0)

# 保存结果
print(f"\n[4] 保存结果到: {output_dir}/")
analyzer.save_results(output_dir)

# 打印摘要
print("\n[5] 分析结果摘要:")
print("-" * 60)
for db_name, df in results.items():
    print(f"    {db_name}: {len(df)} 个显著富集条目")
    if len(df) > 0:
        print(f"        Top 5:")
        for idx, row in df.head(5).iterrows():
            term_id = row.get('Term_ID', row.get('Term_ID', 'N/A'))
            term_name = row.get('Term_Name', 'N/A')
            if len(str(term_name)) > 50:
                term_name = str(term_name)[:50] + "..."
            print(f"            - {term_id}: {term_name}")

print("\n" + "=" * 60)
print("测试完成!")
print("=" * 60)
PYTHON_SCRIPT

echo ""
echo "结果已保存到: $WORKDIR/results/"
ls -la "$WORKDIR/results/"
