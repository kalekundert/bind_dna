#!/usr/bin/env bash
set -euo pipefail

cat <<EOF
- Add mRNA/linker only controls
- Specify that all samples should be diluted to roughly the same concentration, 
  otherwise it's hard to get a good image for densiometry
EOF

# I went back and forth on which conditions to test, and ultimately I decided 
# to simply compare my established protocol to the [Reyes2021] protocol.  
# Depending on what I see, I can make more fine-grained comparisons later.
sw cond ligate_reyes2021_cond.xlsx |

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
#   - Bigger volumes make pipetting easier.
#   - I can save the left-overs and attempt to purify them later.
sw reaction \
  -t "ligation reaction/s" \
  'nuclease-free water  ;         ; to 9.8 µL ; +' \
  'T4 RNA ligase buffer ;     10x ;      1 µL ; +' \
  'ATP                  ;   10 mM ;      1 µL ; +' \
  'f11,f111             ;   10 µM ;      2 µL ; -' \
  'o129,o237            ;   10 µM ;      3 µL ; -' \
  -n 2 |

sw step "Incubate as follows:~90°C for 30s~Cool to 25°C at 1°C/s" |

# I'd rather include the enzyme in a master mix (e.g. anneal with just the 
# oligos and 1x T4 RNA ligase buffer; add water, ATP, more T4 RNA ligase 
# buffer, and enzyme after incubation).  This would avoid having to pipet such 
# small volumes, but for now I'm just going to follow [Reyes2021] exactly.
sw step "Add to each reaction:~0.2 µL 10 U/µL T4 RNA ligase" |

sw step "Incubate at 25°C for 30 min." |

sw step "Label the products: f119, f115" |

sw cdna/make f11,f111 o129,o239 -n 2 -v 2 -c 0 -W -A |

sw step "Label the products: f119, f113" |

sw gel urea/o194 4


