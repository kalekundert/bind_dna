#!/usr/bin/env zsh
set -euo pipefail
alias sw=stepwise

sw gel urea f15 --mix-volume 60 -S |
sw northern/transfer_semidry |
sw northern/optimize_uv |
sw northern/hybridize o135 |
sw laser nir

