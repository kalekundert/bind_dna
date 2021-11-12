#!/usr/bin/env zsh
set -euo pipefail
alias sw=stepwise

# Don't do a whole gel just to optimize cross-linking; way faster to do a dot 
# blot.

#sw gel urea f15 --mix-volume 60 -S |
#sw northern/transfer_semidry |

# RNA quantity
# - Typical gel:
#   - 2.5 µL 160 nM mRNA
# - Dot blot: 2 µL [1] → 200 nM
#
# - 4 steps 10x serial dilution
#   - 10x good for assaying dynamic range.
sw serial 12µL 200nM / 10 4 -m f15 -d 'nuclease-free water' |

sw step "Incubate at 95°C for 3 min." |

sw step "Soak a sheet of HyBond N+ membrane in 1x TBE for ≈1 min" |
sw step "Pipet 2 µL of each dilution onto the membrane." |
sw note "https://tinyurl.com/h1qqi3db" |

sw northern/optimize_uv |
sw northern/hybridize o135 |
sw laser nir

# [1] https://www.abcam.com/protocols/rna-dot-blot-protocol

