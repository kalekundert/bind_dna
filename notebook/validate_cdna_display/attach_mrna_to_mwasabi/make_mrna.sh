#!/usr/bin/env zsh
set -euo pipefail
alias sw=stepwise

sw zap |
sw cdna/anneal 1 -v 16 |
sw cdna/ligate 1 -v 160
