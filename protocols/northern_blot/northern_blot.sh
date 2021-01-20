#!/usr/bin/env zsh
set -euo pipefail
alias sw=stepwise

sw northern_blot/transfer_wet |
sw northern_blot/crosslink |
sw northern_blot/hybridize_ultrahyb o135
