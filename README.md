# AllEnricher v2.0

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Version](https://img.shields.io/badge/Version-2.0.0-orange.svg)](https://github.com/zd105/AllEnricher)

**A comprehensive gene set function enrichment tool for multiple species**

AllEnricher is a modern, high-performance tool for gene set enrichment analysis, supporting multiple databases, algorithms, and species with AI-powered result interpretation.

## ✨ Features

### 🧬 Comprehensive Analysis Support
- **Multiple Databases**: GO, KEGG, Reactome, WikiPathways, MSigDB, DO, DisGeNET
- **Multiple Algorithms**: Fisher's exact test, Hypergeometric test, GSEA, ssGSEA (single-sample GSEA)
- **Multiple Species**: 16 pre-configured model organisms with species lookup support (KEGG API for additional species)
- **Custom Annotation**: Build custom gene set libraries from 3-column TSV files (gene, term_id, term_name)
- **Species Lookup**: Search species by Latin name, KEGG abbreviation, or NCBI taxid interchangeably
- **Database Versioning**: Select specific database versions for reproducible analysis

### 🚀 Performance & Usability
- **Parallel Processing**: Multi-core support for faster analysis (configurable via `n_jobs`)
- **Progress Tracking**: Real-time progress bars and status updates
- **No Gene Limit**: Unlimited gene input size (no hard cap on gene count)

### 📊 Rich Visualization
- **Multiple Plot Types**: Bar plots, bubble plots, dot plots, enrichment maps, network plots, heatmaps, UpSet plots
- **Interactive Reports**: HTML reports with sortable tables and embedded plots
- **Publication-Ready**: High-quality PDF/PNG output
- **Pathway Hyperlinks**: GO/KEGG pathway URLs included in results for quick reference

### 🤖 AI-Powered Interpretation
- **Multiple AI Backends**: OpenAI GPT-4, Anthropic Claude, Local Ollama
- **Biological Insights**: Automated interpretation of enrichment results
- **Term Summaries**: AI-generated descriptions of biological terms

### 🔌 Flexible Integration
- **REST API**: FastAPI-based web service
- **Command Line**: Intuitive CLI with extensive options
- **Python API**: Direct integration into Python workflows

## 📦 Installation

### From PyPI (Recommended)
```bash
pip install allenricher
```

### From Source
```bash
git clone https://github.com/zd105/AllEnricher-v2.git
cd AllEnricher-v2
pip install -e .
```

### With Optional Dependencies
```bash
# For API support
pip install allenricher[api]

# For AI interpretation
pip install allenricher[ai]

# For all features
pip install allenricher[all]
```

## 🚀 Quick Start

### Command Line

```bash
# Basic analysis
allenricher analyze -i genes.txt -s hsa -d GO,KEGG -o results/

# With AI interpretation
allenricher analyze -i genes.txt -s hsa --ai openai --ai-key YOUR_KEY

# Download databases
allenricher download -d GO,KEGG -s hsa

# Start API server
allenricher serve --port 8000
```

### Python API

```python
from allenricher import EnrichmentAnalyzer, Config, DatabaseManager

# Configure analysis
config = Config(
    species="hsa",
    databases=["GO", "KEGG", "Reactome"],
    method="fisher",
    qvalue_cutoff=0.05
)

# Load databases
db_manager = DatabaseManager("./database", "hsa")
db_manager.load_databases(["GO", "KEGG"])

# Run analysis
analyzer = EnrichmentAnalyzer(config)
gene_set = {"BRCA1", "TP53", "EGFR", "MYC", "KRAS"}
background = db_manager.get_background_genes()

results = analyzer.run_analysis(
    gene_set=gene_set,
    background_set=background,
    database_data=db_manager.get_all_term_data()
)

# Save results
analyzer.save_results("./output")
```

### REST API

```bash
# Start server
allenricher serve --port 8000

# Submit analysis
curl -X POST "http://localhost:8000/api/analyze" \
  -H "Content-Type: application/json" \
  -d '{
    "genes": ["BRCA1", "TP53", "EGFR"],
    "species": "hsa",
    "databases": ["GO", "KEGG"]
  }'

# Check status
curl "http://localhost:8000/api/status/{job_id}"
```

## 📖 Documentation

### Supported Databases

| Database | Description | Species Support |
|----------|-------------|-----------------|
| GO | Gene Ontology | 16 model organisms (expandable via build) |
| KEGG | KEGG Pathways | 16 model organisms (expandable via build) |
| Reactome | Reactome Pathways | 16 model organisms |
| WikiPathways | WikiPathways | Multiple species |
| MSigDB | Molecular Signatures | Human |
| DO | Disease Ontology | Human |
| DisGeNET | Disease-Gene Associations | Human |

### Supported Methods

| Method | Description | Use Case |
|--------|-------------|----------|
| Fisher | Fisher's exact test | Standard enrichment |
| Hypergeometric | Hypergeometric test | Alternative to Fisher |
| GSEA | Gene Set Enrichment Analysis | Ranked gene lists |
| SSGSEA | Single-sample GSEA | Single sample analysis |

### CLI Commands

```bash
# Analyze gene list
allenricher analyze [OPTIONS]

Options:
  -i, --input FILE       Input gene list file [required]
  -s, --species TEXT     Species code (default: hsa)
  -d, --databases TEXT   Comma-separated databases (default: GO,KEGG)
  -o, --output DIR       Output directory (default: ./results)
  -m, --method TEXT      Enrichment method: fisher/hypergeometric/gsea/ssgsea (default: fisher)
  -c, --correction TEXT  Multiple testing correction (default: BH)
  -p, --pvalue FLOAT     P-value cutoff (default: 0.05)
  -q, --qvalue FLOAT     Q-value cutoff (default: 0.05)
  --ai TEXT              AI backend (openai/claude/ollama/mock)
  --ai-key TEXT          AI API key

# Download databases
allenricher download -d GO,KEGG -s hsa

# Build species databases (download and build local GO/KEGG/Reactome databases)
allenricher build -s hsa -t 9606 -d GO,KEGG

Options:
  -s, --species TEXT     Species KEGG code (e.g., hsa, mmu) [required]
  -t, --taxonomy INT     NCBI Taxonomy ID (e.g., 9606 for human) [required]
  -d, --databases TEXT   Comma-separated databases to build (default: GO,KEGG,Reactome)
  --database-dir DIR     Database output directory (default: ./database)
  --gene-info FILE       Path to NCBI gene_info.gz file (for GO database building)

# List resources
allenricher list species               # List pre-configured model organisms
allenricher list databases             # List available databases

# Generate config
allenricher config -o my_config.yaml
```

## 📊 Output

### Result Files
- `*_enrichment.tsv` - Tab-separated enrichment results
- `plots/*.pdf` - Publication-ready plots
- `report.html` - Interactive HTML report
- `ai_interpretation.json` - AI-generated interpretations

### Result Columns

| Column | Description |
|--------|-------------|
| Term_ID | GO/Pathway identifier |
| Term_Name | Term/Pathway name |
| Gene_Count | Number of genes in term |
| Rich_Factor | Enrichment ratio |
| P_Value | Raw p-value |
| Adjusted_P_Value | Multiple testing corrected p-value |
| Genes | List of enriched genes |
| Term_URL | Direct hyperlink to the GO/KEGG pathway page for quick reference |

## 🔧 Configuration

Create a YAML configuration file:

```yaml
# allenricher.yaml
species: "hsa"
databases:
  - "GO"
  - "KEGG"
  - "Reactome"
method: "fisher"
correction: "BH"
pvalue_cutoff: 0.05
qvalue_cutoff: 0.05
min_genes: 2

# GSEA/ssGSEA parameters
gsea_min_size: 10       # Minimum gene set size for GSEA/ssGSEA (default: 10)
gsea_max_size: 500      # Maximum gene set size for GSEA/ssGSEA (default: 500)

# Gene set filtering
min_genes: 2           # Minimum genes per term (default: 2)
max_genes: .inf        # Maximum genes per term, .inf means unlimited (default: unlimited)

# Visualization
plot_formats:
  - "pdf"
  - "png"
top_terms: 20

# Performance
n_jobs: 1               # Number of parallel workers (1 = sequential, -1 for all cores)

# AI
ai_interpretation: true
ai_model: "gpt-4"
```

Use with:
```bash
allenricher analyze --config allenricher.yaml -i genes.txt
```

## 🧪 Testing

```bash
# Run tests
pytest tests/

# With coverage
pytest --cov=allenricher tests/

# Type checking
mypy allenricher
```

## 📝 Citation

If you use AllEnricher in your research, please cite:

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

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🤝 Contributing

Contributions are welcome! Please read our [Contributing Guide](CONTRIBUTING.md) for details.

## 📧 Contact

- **Issues**: [GitHub Issues](https://github.com/zd105/AllEnricher/issues)
- **Email**: allenricher@example.com

## 🙏 Acknowledgments

- Original AllEnricher v1.0 by Du Zhang et al.
- Gene Ontology Consortium
- KEGG Database
- Reactome Project
- All contributors and users
