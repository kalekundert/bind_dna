#!/usr/bin/env zsh
set -euo pipefail
alias sw=stepwise

sw zap |
sw purex f88 f89 -w 30m -v 6 -r -T 333 |
sw cdna/couple -n 2 -v 2.5 |
sw gel sds/o194 4 -S |
sw laser blue red
