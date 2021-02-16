#!/usr/bin/env bash
set -euo pipefail

sw zap |
sw serial 2.0 '10 µM' to 0.1 7 -0 -m f85 -d 'nuclease-free water' |
sw ivtt f85 -p nebexpress -v 5 -n 8 -c 500 -C 10000  |
sw gel bolt 8

