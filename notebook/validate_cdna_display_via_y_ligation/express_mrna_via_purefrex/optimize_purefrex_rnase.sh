#!/usr/bin/env bash
set -euo pipefail

sw zap |

# Volume:
# - f15 aliquots are 5 µL (I think), so I tried to increase the volume as much 
#   as possible (to make pipetting more accurate) without using more than one 
#   aliquot.
sw serial 2.6 '5000 nM' to 50 7 -0 -m f15 -d 'nuclease-free water' |

# Reaction volume:
# - I could reduce to 2.5 µL and still keep all pipetting steps above 0.5 µL, 
#   but I think I'm just a bit more comfortable with 5 µL.
# mRNA concentration:
# - GeneFrontier recommends 10-1000 nM.
# - I'll test exactly that; all the numbers work out well.

sw ivtt f15 -p frex/rnase -v 5 -n 8 -c 1000 -C 5000  |

sw gel bolt 8 -S |
sw laser blue red

