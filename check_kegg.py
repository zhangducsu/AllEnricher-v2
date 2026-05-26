import pandas as pd

kegg = pd.read_csv('results/KEGG_enrichment.tsv', sep='\t')

print('=== KEGG Term_Name 层级分布 ===')
for _, row in kegg.head(20).iterrows():
    name = row['Term_Name']
    levels = name.count('|') + 1
    print(f"  {row['Term_ID']}: {levels}层 - {name[:70]}")

print()
print('=== 统计 ===')
one_level = sum(1 for name in kegg['Term_Name'] if name.count('|') == 0)
two_level = sum(1 for name in kegg['Term_Name'] if name.count('|') == 1)
three_level = sum(1 for name in kegg['Term_Name'] if name.count('|') == 2)
print(f"  1层: {one_level}")
print(f"  2层: {two_level}")
print(f"  3层: {three_level}")
