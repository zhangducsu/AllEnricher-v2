#!/usr/bin/env python3
"""检查数据库文件的内容"""

import gzip
from pathlib import Path


def check_file():
    v1_db_dir = Path(__file__).parent.parent / "AllEnricher-v1" / "database" / "organism" / "v20190612" / "hsa"
    
    print("=" * 80)
    print("检查 hsa.GO2gene.tab.gz")
    print("=" * 80)
    filepath = v1_db_dir / "hsa.GO2gene.tab.gz"
    
    with gzip.open(filepath, 'rt', encoding='utf-8') as f:
        lines = [next(f) for _ in range(10)]
        print("\n".join(lines))
    
    print("\n" + "=" * 80)
    print("检查 hsa.KEGG2gene.tab.gz")
    print("=" * 80)
    filepath = v1_db_dir / "hsa.KEGG2gene.tab.gz"
    
    with gzip.open(filepath, 'rt', encoding='utf-8') as f:
        lines = [next(f) for _ in range(10)]
        print("\n".join(lines))


if __name__ == "__main__":
    check_file()
