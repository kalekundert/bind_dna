#!/usr/bin/env zsh
set -euo pipefail
alias sw=stepwise

sw reactions fresh_aliquot.csv -k "-80: Stored at -80°C and freeze-thawed ≈5-10 times over ≈2 months." -k "-20: Thawed once and stored at -20°C for a few days." -k "fresh: Freshly thawed for the first time for this experiment." |
sw dilute_dcas9.txt |
sw huh f12 -n 9 -c 1 -m buffer,sgrna,dna |
sw gel sdsmax 10 -r 60 -S |
sw gelgreen_coom



