#!/usr/bin/env zsh
set -euo pipefail
alias sw=stepwise

sw reactions repeat_a15_linker.csv -k 'A15: Rep ori with a 15-adenosine spacer (f107)' |
sw vegarocha |
sw huh f107 -n 3 -c 1.84 -S -b 'Vega-Rocha Rep buffer' -B 2.5 |
sw gel sdsmax -S 4 |
sw gelgreen_coom
