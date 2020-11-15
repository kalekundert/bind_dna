#!/usr/bin/env zsh
set -euo pipefail
alias sw=stepwise

sw zap |
sw reactions 
sw cdna/make f97 o129 |
sw step "Label the product: f98" |
sw purex f98 -r -T333 -v 3.6 |
