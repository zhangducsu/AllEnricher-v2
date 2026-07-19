#!/bin/bash
# AllEnricher v2 Docker Test Script

set -e

# Check Docker for running
if ! docker info > /dev/null 2>&1; then
    echo "Error: Docker is not running, please start Docker Desktop first"
    exit 1
fi

# Working Directory
WORKDIR="$(cd "$(dirname "$0")" && pwd)"

echo "=========================================="
echo "AllEnricher v2 Docker Test"
echo "Test with v1 example data"
echo "=========================================="

# Check Files
if [ ! -f "$WORKDIR/../AllEnricher-v1/example/example.glist" ]; then
    echo "Error: Gene List File does not exist"
    exit 1
fi

if [ ! -d "$WORKDIR/../AllEnricher-v1/database/organism/v20190612/hsa" ]; then
    echo "Error: Database directory does not exist"
    exit 1
fi

# Create Results Directory
mkdir -p "$WORKDIR/results"

# Run test container
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

# Add v2 to Path
sys.path.insert(0, '/workspace/allenricher')

from pathlib import Path
from allenricher import EnrichmentAnalyzer, Config
from allenricher.database.manager import DatabaseManager

print("=" * 60)
print("AllEnricher v2 test")
print("=" * 60)

# Configure
gene_list = "/workspace/example/example.glist"
database_dir = "/workspace/database/organism/v20190612/hsa"
output_dir = "/workspace/results"
species = "hsa"

# Check Files
if not os.path.exists(gene_list):
    print(f"Error: Genetic list file does not exist: {gene_list}")
    sys.exit(1)

if not os.path.exists(database_dir):
    print(f"Error: Database directory does not exist: {database_dir}")
    sys.exit(1)

# Create Configuration
config = Config(
    species=species,
    databases=["GO", "KEGG"],
    method="fisher",
    pvalue_cutoff=0.05,
    qvalue_cutoff=0.05,
    min_genes=2
)

# Create Output Directory
Path(output_dir).mkdir(parents=True, exist_ok=True)

# Loading list of genes
print(f"\n[1] Loading of gene lists: {gene_list}")
analyzer = EnrichmentAnalyzer(config)
gene_set = analyzer.load_gene_list(gene_list)
print(f"Gene count: {len(gene_set)}")

# Loading Database
print(f"\n[2] Loading of database: {database_dir}")
db_manager = DatabaseManager(database_dir, species)
db_manager.load_databases(config.databases)
background_set = db_manager.get_background_genes()
print(f"Number of background genes: {len(background_set)}")

# Access to database data
database_data = db_manager.get_all_term_data()
print(f"Loaded database: {list(database_data.keys())}")

# Run Analysis
print("\n[3] Run the enrichment analysis...")
results = analyzer.run_analysis(gene_set, background_set, database_data)

# Inspection results
if not results or len(results) == 0:
    print("No significant enrichment results found")
    sys.exit(0)

# Save Results
print(f"\n[4] Save results to: {output_dir}/")
analyzer.save_results(output_dir)

# Print Summary
print("\n[5] Summary of the results of the analysis:")
print("-" * 60)
for db_name, df in results.items():
    print(f"{db_name}: {len(df)}A significant rich entry")
    if len(df) > 0:
        print(f"        Top 5:")
        for idx, row in df.head(5).iterrows():
            term_id = row.get('Term_ID', row.get('Term_ID', 'N/A'))
            term_name = row.get('Term_Name', 'N/A')
            if len(str(term_name)) > 50:
                term_name = str(term_name)[:50] + "..."
            print(f"            - {term_id}: {term_name}")

print("\n" + "=" * 60)
print("Test complete!")
print("=" * 60)
PYTHON_SCRIPT

echo ""
echo "Result saved to: $WORKDIR/results/"
ls -la "$WORKDIR/results/"
