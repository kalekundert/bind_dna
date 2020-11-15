#!/usr/bin/env zsh
set -euo pipefail
alias sw=stepwise

# - Clone plasmid
#   
#   sw make p175
#
# - Make IVT
#
#   sw make f94 f95
#
# - Test broccoli
#
#   - PAGE
#   - Stain with DFHBI-1T
#
# - Test PURExpress
#   
#   - Attach linker (f98)
#     - sw cdna/make f97 o129
#   - PURExpress
#     - sw purex
#     - timecourse?
#   - Image
#     - sw stain_dfhbi
#
sw zap |
sw cdna/make f97 o129 |
sw step "Label the product: f98" |
sw reactions check_degradation_via_broccoli.csv |
sw purex f98 -n 3 -r -v 3.6 |
sw gel sds/o194 6 -S |
sw laser blue red |
sw filonov2015

