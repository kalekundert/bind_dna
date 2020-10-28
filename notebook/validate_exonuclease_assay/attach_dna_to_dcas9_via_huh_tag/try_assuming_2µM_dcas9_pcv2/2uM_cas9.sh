#!/usr/bin/env zsh
set -euo pipefail
alias sw=stepwise

sw make f12 f16 |
sw pcr_cleanup |
sw reactions 2uM_cas9_reactions.xlsx |
sw huh f12 f16 -n 4 -S -c 1.84 -b "Mn²⁺ buffer" -B 2.5 |
sw gel sdsmax 7 -r60 -S |
sw coomassie_gelgreen.txt |
sw laser blue red
