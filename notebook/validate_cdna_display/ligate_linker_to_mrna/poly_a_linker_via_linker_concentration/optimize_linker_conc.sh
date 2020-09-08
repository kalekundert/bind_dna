#!/usr/bin/env zsh
set -euo pipefail
alias sw=stepwise

sw zap |\
  sw samples.txt |\
  sw serial 1 80µM 10 4 -m o129 |\
  sw cdna/anneal 4 f11 o129 -v 4 -m mrna -L '10–80 µM' |\
  sw cdna/liga 4 -v 10 |\
  sw gel urea/o194 3
