#!/usr/bin/env python3
"""AllEnricher v2.0 - Command Line Entry Module

Support for adoption `python -m allenricher` Run Directly
"""

import sys
from allenricher.cli import main

if __name__ == '__main__':
    sys.exit(main())
