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
#     - sw cdna/make f95 o129
#   - PURExpress
#     - sw purex
#     - timecourse?
#   - Image
#     - sw stain_dfhbi
#
sw zap |
sw 

