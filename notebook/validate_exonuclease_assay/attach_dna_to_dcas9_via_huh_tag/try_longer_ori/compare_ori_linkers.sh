#!/usr/bin/env zsh
set -euo pipefail
alias sw=stepwise

sw reactions compare_ori_linkers.csv -k 'PCV: Rep ori with native context from the PCV2 genome (f106)' -k 'A15: Rep ori with a 15-adenosine spacer (f107)' |
sw vegarocha |
sw huh f105 f106 f107 -n 7 -c 1.84 -S -b 'Vega-Rocha Rep buffer' -B 2.5 |
sw gel sdsmax -S 10 |
sw gelgreen -r |
sw simplyblue -fr
