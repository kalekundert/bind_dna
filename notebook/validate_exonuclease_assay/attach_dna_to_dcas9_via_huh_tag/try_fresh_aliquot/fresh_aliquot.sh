#!/usr/bin/env zsh
set -euo pipefail
alias sw=stepwise

sw reactions fresh_aliquot.csv |
sw dilute_dcas9.txt |
sw huh f12 -n 9 -c 1 -m buffer,sgrna,dna |
sw gel sdsmax 10 -S |
sw gelgreen_coom



