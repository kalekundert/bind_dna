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
sw step "Prepare the +DNA control:~7.24 µL water~2 µL 5x Zif268 storage buffer + MnCl₂~0.76 µL 1314 nM f107" |
sw huh PCV2-Zif268 f107 -n7 -P 16 -b "Zif268 storage buffer + MnCl₂" -C 5 |

# Stain:
# - I've tried using SYPRO Ruby to more sensitively detect the protein, but 
#   didn't get good results.  Not sure why.
# - Once I make f134, I can try staining with SYPRO Orange and directly imaging 
#   both channels.  I don't know if this will work, though.
sw gel bolt/pcv2 7
