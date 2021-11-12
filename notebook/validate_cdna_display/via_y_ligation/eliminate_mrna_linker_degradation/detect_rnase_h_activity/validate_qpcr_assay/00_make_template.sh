#!/usr/bin/env bash
set -euo pipefail

sw make f108 |
sw spin_cleanup |
sw ivt f108 -C |
sw zymo_clean -d pre -s 20 |
sw aliquot '5 µL' '10 µM'
