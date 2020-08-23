#!/usr/bin/env zsh
set -euo pipefail
alias sw=stepwise

sw zap |
  sw cdna/anneal 1 f11 o129 -v 4 |
  sw cdna/ligate 7 -v 4.4 -i '0 min, 10 min, 2h, 16h' -I '16°C, 25°C' -m lig,rna -Q |
  sw custom "When each timepoint is reached, quench the corresponding reactions as follows:~Add 0.2 µL 500 mM EDTA~Freeze at -20°C" |
  sw gel anneal 7
