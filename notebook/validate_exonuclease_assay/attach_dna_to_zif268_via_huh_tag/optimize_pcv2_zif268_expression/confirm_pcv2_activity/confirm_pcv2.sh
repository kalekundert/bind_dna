#!/usr/bin/env bash
set -euo pipefail

sw cond confirm_pcv2_cond.xlsx |

# Buffer:
# - Don't use Vega-Rocha buffer.
# - Use Zif268 storage buffer with MnCl₂.
sw zif_mncl2_buffer.txt |

# Controls:
# - Include the fusion protein in the master mix, because it's very 
#   concentrated.
sw step "Prepare the +DNA control:~7.26 µL water~2 µL 5x Zif268 storage buffer + MnCl₂~0.74 µL 1360 nM f12" |
sw huh PCV2-Zif268 f12 -n4 -P 1.81 -b "Zif268 storage buffer + MnCl₂" -C 5 |

sw gel bolt/pcv2 4 -S | sw stain gelgreen | sw stain sypro-ruby/microwave
