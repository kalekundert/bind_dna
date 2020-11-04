#!/usr/bin/env zsh
set -euo pipefail
alias sw=stepwise

sw make f105 f106 f107 |
sw pcr_cleanup |
sw e_gel |
sw step "Enter concentrations in the database."
