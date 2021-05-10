#!/usr/bin/env zsh
set -euo pipefail
alias sw=stepwise

sw zap |\
  sw cdna/anneal 1 f11 o129 -v 2.2 |\
  sw samples.txt |\
  sw cdna/liga 2 -v 10 -P -m lig,rna |\
  sw gel anneal 3
