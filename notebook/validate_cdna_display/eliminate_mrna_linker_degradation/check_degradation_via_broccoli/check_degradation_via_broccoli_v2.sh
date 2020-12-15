#!/usr/bin/env zsh
set -euo pipefail
alias sw=stepwise

sw zap |
sw cdna/anneal f85,f97 o129 -R5 -x 0.6 |
sw cdna/ligate -n 2 |
sw step "Label the products: f89,f98" |
sw reactions check_degradation_via_broccoli_v2.csv |
sw purex f89,f98 -n 5 -r -T 125 -v 3.6 |
sw gel sds/o194 10 -S |
sw laser blue red |
sw filonov2015 |
sw gelgreen

