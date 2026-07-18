#!/usr/bin/env python3
"""
Compares the drawing output of v1 and v2
"""

from pathlib import Path

v1_dir = Path(r"F:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v1\example\allenricher\fisher\Q0.05")
v2_dir = Path(r"F:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2\test_output\v1v2_compare\plots")

print("="*70)
print("Draw File Comparison")
print("="*70)

# Drawing file for v1
v1_plots = list(v1_dir.glob("*.pdf"))
print(f"\nDrawing file for v1 ({len(v1_plots)}(a) The number of persons:")
for f in sorted(v1_plots):
    size = f.stat().st_size
    print(f"  {f.name:<50} {size:>10,} bytes")

# Drawing file for v2
v2_plots = list(v2_dir.glob("*.pdf"))
print(f"\nv2 Drawing Files ({len(v2_plots)}(a) The number of persons:")
for f in sorted(v2_plots):
    size = f.stat().st_size
    print(f"  {f.name:<50} {size:>10,} bytes")

# Comparison
print("\n" + "="*70)
print("A comparative summary of the mapping output:")
print("="*70)

# Grouped by Database
databases = ['GO', 'KEGG', 'Reactome', 'DO', 'DisGeNET']
for db in databases:
    v1_bar = list(v1_dir.glob(f"*.{db}_barplot*.pdf"))
    v1_bubble = list(v1_dir.glob(f"*.{db}_bubbleplot*.pdf"))
    v2_bar = list(v2_dir.glob(f"{db}_barplot*.pdf"))
    v2_bubble = list(v2_dir.glob(f"{db}_bubble*.pdf"))
    
    print(f"\n{db}:")
    print(f"v1: column charts ({len(v1_bar)}), bubble charts ({len(v1_bubble)})")
    print(f"v2: column ({len(v2_bar)}), bubble ({len(v2_bubble)})")
    
    if len(v1_bar) == len(v2_bar) and len(v1_bubble) == len(v2_bubble):
        print(f"* Same type of drawing")
    else:
        print(f"The number of drawings is inconsistent")

print("\n" + "="*70)
