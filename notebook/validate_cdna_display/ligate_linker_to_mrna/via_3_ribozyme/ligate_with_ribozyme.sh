#!/usr/bin/env zsh
set -euo pipefail
alias sw=stepwise

sw reactions ligate_with_ribozyme.xlsx |
sw zap |
sw cdna/anneal 2 f92,f94 o129 -m link | 
sw cdna/ligate 2 -v 10 -MQ |
sw gel urea/o194 7
