import pandas as pd
import json
from pathlib import Path

# 加载结果
go_df = pd.read_csv('results/GO_enrichment.tsv', sep='\t')
kegg_df = pd.read_csv('results/KEGG_enrichment.tsv', sep='\t')

# 读取 AI 解读
with open('results/ai_interpretation.json', 'r', encoding='utf-8') as f:
    ai_interpretation = json.load(f)

results = {'GO': go_df, 'KEGG': kegg_df}

with open('example_genes.txt', 'r') as f:
    gene_list = [line.strip() for line in f if line.strip()]

from allenricher.report.generator import ReportGenerator
generator = ReportGenerator(output_dir='results')
html = generator.generate(
    results=results,
    output_file='results/report.html',
    gene_list=gene_list,
    ai_interpretation=ai_interpretation
)
print('HTML 已重新生成')
