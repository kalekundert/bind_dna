#!/usr/bin/env zsh
set -euo pipefail
alias sw=stepwise

sw zap |
sw cdna/anneal 2 f85 o128,o129 -m mrna -v 8 -x 0.6 -R 5 |
sw note "Using 0.6x linker reduces the amount of unligated linker, see expt #1." "Setup \\d+ annealing reactions?" |
sw cdna/ligate 2 -v 80 |
sw cdna/wash |
sw custom "Label the products f88,f89" |
sw aliquot '4 µL' '333 nM' |
sw note "This concentration refers to the ligated species (f88 and f89).  There is also an approximately equal amount of unligated mRNA (f85) in these reactions."
