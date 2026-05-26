#!/usr/bin/env python3
"""AllEnricher v2.0 - 命令行入口模块

支持通过 `python -m allenricher` 直接运行
"""

import sys
from allenricher.cli import main

if __name__ == '__main__':
    sys.exit(main())
