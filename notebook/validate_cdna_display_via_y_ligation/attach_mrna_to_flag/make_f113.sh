#!/usr/bin/env bash
set -euo pipefail

# My standard protocol for ligating the linker to the mRNA doesn't work for 
# f113, because it's small enough to pass through to 100 kDa MWCO spin filter.
# This is an alternative protocol that uses ethanol precipitation to 
# concentrate the mRNA after the ligation.  This protocol is modeled after 
# [Reyes2021], except that they also include a gel extraction step which I 
# don't have the equipment for at the moment.

sw zap |

# mRNA quantity:
# - [Reyes2021] used 2 µg mRNA/linker in a 100 µL translation reaction.
# - f113 MW: ≈57257.87
# - final conc: 350 nM
# - this is in-line with my results from expt #99 (i.e. that 200-2000 nM 
#   gives >70% expression), albeit on the low end.
# - the stock concentration will have to be >1000 nM to reach this.
# - probably want at least 2000 nM to be comfortable, maybe 5000 nM to try 
#   higher concentrations.
sw cdna/anneal f111 o237 -l 0.6 -r 10 |

# mRNA conc:
# - Calculated by hand for '-r 10' in previous step.
# - 5.62 µM ≈ 10 µM * 10 µL / 17.78 µL
sw cdna/ligate 17.78 5.624296962879639 |

sw step "Label the product: f113" |

# Preset:
# - I was choosing between:
#   - 'pcr': right length, wrong backbone/strandedness
#   - 'microrna': right backbone/strandedness, wrong length
#
# - It turns out that both protocols are very similar, differing only in how 
#   much ethanol is added (microrna: more) and the minimum centrifugation speed 
#   (microrna: faster).  And I'll do the centrifugation steps at max speed 
#   either way.
#
# - I decided to use 'pcr' because I think the differences were mainly geared 
#   towards the shorter length of micrornas, not the different backbone.
sw ethanol -p pcr -v 10 |

sw aliquot "5 µL" "2 µM"
