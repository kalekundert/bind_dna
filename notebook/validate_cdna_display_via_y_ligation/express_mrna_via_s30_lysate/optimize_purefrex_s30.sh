#!/usr/bin/env bash
set -euo pipefail

sw zap |

# Volume:
# - f15 aliquots are 5 µL (I think), so I tried to increase the volume as much 
#   as possible (to make pipetting more accurate) without using more than one 
#   aliquot.
# - I need 1 µL for each reaction.  2.6 µL total is cutting it close, but 
#   should be enough.
sw serial 2.6 '5000 nM' to 50 7 -0 -m f15 -d 'nuclease-free water' |

# Reaction volume:
# - I could reduce to 2.5 µL and still keep all pipetting steps above 0.5 µL, 
#   but I think I'm just a bit more comfortable with 5 µL.
# mRNA concentration:
# - GeneFrontier recommends 10-1000 nM.
# - I'll test exactly that; all the numbers work out well.

sw ivtt f15 -p purefrex -v 5 -n 8 -c 1000 -C 5000  |

# mRNA concentration:
# - Promega recommends 10-100 µg for a 50 µL reaction.
# - That corresponds to 800-8000 nM f15 (MW=249047.71 Da).
# - I'd need a stock of 40 µM to reach a final concentration of 8 µM.
# - My stock is only 5 µM.  I'd have make new RNA specifically for this 
#   experiment, and even then it'd be hard to make it so concentrated.  I'd 
#   probably need to do an ethanol precipitation.
# - And after that I'd also have to get f89 really concnetrated, as well.
# - I'm just going to use what I have (e.g. the dilutions from above), and see 
#   how things work.

sw ivtt f15 -p promega-s30-linear -v 5 -n 8 -c 1000 -C 5000  |

sw gel bolt 16 -S |
sw laser blue red

