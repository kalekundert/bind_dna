#!/usr/bin/env bash
set -euo pipefail

# I went back and forth on which conditions to test, and ultimately I decided 
# to simply compare my established protocol to the [Reyes2021] protocol.  
# Depending on what I see, I can make more fine-grained comparisons later.
sw cond ligate_undiluted_salt_cond.xlsx \
  -k "F: FLAG (f115 = f111 + o237)" \
  -k "Z: Zif268 (f119 = f11 + o129)" \
  |

sw zap |

sw reaction \
  -s "+linker control/s" \
  'nuclease-free water  ;         ; to 17.6 µL ; +' \
  'linker               ;   10 µM ;     0.8 µL ; -' \
  |

sw reaction \
  -s "+mRNA control/s" \
  'nuclease-free water  ;         ; to 17.6 µL ; +' \
  'mRNA                 ;   10 µM ;     0.8 µL ; -' \
  |

sw step "Keep these controls on ice until ready to run the gel." |

# mRNA/linker concentrations:
# - [Reyes2021] has final concentrations:
#   - mRNA: 4 µM
#   - linker: 6 µM
#
# - My stocks aren't concentrated enough to reach this, so I'll do:
#   - mRNA: 2 µM
#   - linker: 3 µM
#
# mRNA/linker volumes:
# - I chose the volumes to use a full aliquot of each mRNA.
# - This is wasteful in that I'll make much more annealed mRNA than I need to 
#   run a gel, but:
#   - Bigger volumes make pipetting easier/more accurate.
#   - I can save the left-overs and attempt to purify them later.
#
# mRNA/linker ratio:
# - Normally I use 0.6x linker, to minimize the amount of unreacted linker.
# - But for these optimization experiments, I should use 1x to be sure that the 
#   linker concentration is not what's limiting yield.
# - I could even argue for using 1.5x linker, like [Reyes2021].  But I think 1x 
#   is fine.
sw cdna/anneal f11,f111 o129,o237 -n 2 -r 3 |

# Extra master mix:
# - Scaled to match the concentrated reaction below.
sw reaction \
  -s "−ligase control/s" \
  'nuclease-free water  ;         ;  to 10 µL ; +' \
  'T4 RNA ligase buffer ;     10x ;      1 µL ; +' \
  'ATP                  ;   10 mM ;      1 µL ; +' \
  'annealed mRNA/linker ; 4.54 µM ;    4.4 µL ; -' \
  -v 4 \
  -X 2.5 \
  -n 2 |

# Extra master mix:
# - Scale master mix to make pipetting the ligase more accurate.
# - All of the master mix components are relatively cheap, so I prefer to make 
#   a lot of extra master mix if it means more accurate pipetting.
sw reaction \
  -s "concentrated ligation reaction/s" \
  'nuclease-free water  ;         ;  to 10 µL ; +' \
  'T4 RNA ligase buffer ;     10x ;      1 µL ; +' \
  'ATP                  ;   10 mM ;      1 µL ; +' \
  'T4 RNA ligase        ; 10 U/µL ;    0.2 µL ; +' \
  'annealed mRNA/linker ; 4.54 µM ;    4.4 µL ; -' \
  -v 4 \
  -X 0.5 \
  -n 2 |

# Ligase volume:
# - Very small volume of ligase; I won't be able to pipet this accurately.
# - But I don't want to scale the reaction; that would use too much material.
# - Adding more ligase doesn't seem to make the reaction go further, so I might 
#   just try to be sure to err on the side of too much.
sw cdna/ligate 1.76 4.54 -n 2 |

sw step "Label the products: f119, f115" |

sw step "Dilute the concentrated reactions to 17.60 µL~4 µL reactions: Add 13.6 µL nuclease-free water~Final mRNA/linker concentration: 454 nM" |

sw gel urea/o194 10 -c 454 -s sybr_cy5

