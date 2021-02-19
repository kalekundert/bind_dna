#!/usr/bin/env zsh
set -euo pipefail
alias sw=stepwise

# Don't do a whole gel just to optimize cross-linking; way faster to do a dot 
# blot.

#sw gel urea f15 --mix-volume 60 -S |
#sw northern/transfer_semidry |

# f15
sw serial 12 -m f15 -d nuclease-free water
sw step | 

sw northern/optimize_uv |
sw northern/hybridize o135 |
sw laser nir

# https://www.abcam.com/protocols/rna-dot-blot-protocol

