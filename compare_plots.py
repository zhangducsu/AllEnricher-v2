#!/usr/bin/env python3
"""
对比v1和v2的绘图输出
"""

from pathlib import Path

v1_dir = Path(r"F:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v1\example\allenricher\fisher\Q0.05")
v2_dir = Path(r"F:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2\test_output\v1v2_compare\plots")

print("="*70)
print("绘图文件对比")
print("="*70)

# v1的绘图文件
v1_plots = list(v1_dir.glob("*.pdf"))
print(f"\nv1绘图文件 ({len(v1_plots)} 个):")
for f in sorted(v1_plots):
    size = f.stat().st_size
    print(f"  {f.name:<50} {size:>10,} bytes")

# v2的绘图文件
v2_plots = list(v2_dir.glob("*.pdf"))
print(f"\nv2绘图文件 ({len(v2_plots)} 个):")
for f in sorted(v2_plots):
    size = f.stat().st_size
    print(f"  {f.name:<50} {size:>10,} bytes")

# 对比
print("\n" + "="*70)
print("绘图产出对比总结:")
print("="*70)

# 按数据库分组
databases = ['GO', 'KEGG', 'Reactome', 'DO', 'DisGeNET']
for db in databases:
    v1_bar = list(v1_dir.glob(f"*.{db}_barplot*.pdf"))
    v1_bubble = list(v1_dir.glob(f"*.{db}_bubbleplot*.pdf"))
    v2_bar = list(v2_dir.glob(f"{db}_barplot*.pdf"))
    v2_bubble = list(v2_dir.glob(f"{db}_bubble*.pdf"))
    
    print(f"\n{db}:")
    print(f"  v1: 柱状图({len(v1_bar)}个), 气泡图({len(v1_bubble)}个)")
    print(f"  v2: 柱状图({len(v2_bar)}个), 气泡图({len(v2_bubble)}个)")
    
    if len(v1_bar) == len(v2_bar) and len(v1_bubble) == len(v2_bubble):
        print(f"  ✓ 绘图类型一致")
    else:
        print(f"  ⚠ 绘图数量不一致")

print("\n" + "="*70)
