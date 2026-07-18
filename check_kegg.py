import pandas as pd

kegg = pd.read_csv('results/KEGG_enrichment.tsv', sep='\t')

print('= = KEG Term_Name Tier Distribution = =')
for _, row in kegg.head(20).iterrows():
    name = row['Term_Name']
    levels = name.count('|') + 1
    print(f"{row['Term_ID']}: {levels}Layer-{name[: 70]}")

print()
print('== sync, corrected by elderman ==')
one_level = sum(1 for name in kegg['Term_Name'] if name.count('|') == 0)
two_level = sum(1 for name in kegg['Term_Name'] if name.count('|') == 1)
three_level = sum(1 for name in kegg['Term_Name'] if name.count('|') == 2)
print(f"Floor 1: {one_level}")
print(f"Floor 2: {two_level}")
print(f"3rd Floor: {three_level}")
