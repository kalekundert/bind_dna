#!/usr/bin/env bash
set -euo pipefail

# Outline
# - Make mRNA
# - Add toehold
# - incubate
# - treat with RNase H
# - PAGE

sw zap |
sw cond \
  protect_mrna_with_o224_cond.xlsx \
  -k "dna/dna oligo values are molar ratios relative to mRNA+linker" |

# If I don't have any f89; I might have some.
# KBK: I have some, but it's only 333 nM, not 1 µM.
# KBK: I made some, but only got 500 nM.
#      I modified the rest of the protocol to account for that.
#sw cdna/make_mrna f85 o129 -A |
#sw step "Label the product: f89" |

sw serial 9 '100µM' / 10 3 -0 -m o224 |

# 7.4 µL total volume:
# - 3.2 µL for each RNase H reaction (2)
# - 1 µL extra
#
# 300 nM final mRNA concentration
# - Would rather this be higher; as high as possible.
sw reaction \
  "PBS  ;      10x ; 0.30 µL ; +" \
  "f89  ;   500 nM ; 1.80 µL ; +" \
  "o224 ; 0-100 µM ; 0.90 µL ; -" \
  -n 4 -v 7.4 \
  -t "annealing" |

sw step "Incubate at 95°C for 2 min, then let cool at room temperature." |

# 2 master mixes, +/- RNase H
# - Same mRNA conc as expression reactions [#18]
#
#   - Optimal expression [#18]:
#     160 nM = 0.8 µL 1 µM mRNA in a 5 µL rxn
#
#   - This reaction:
#     3.2 µL = 6 µL 160 nM mRNA with 300 nM stock.
#
# - RNase H
#   - 1 U RNase H: digest 20 pmol in 20 min at 37°C [#77]
#   - 0.96 pmol mRNA = 3.2 µL * 0.3 µM
#   - 0.048 U = 0.96 pmol / 20 pmol/U
#   - 10x excess: 0.048 U/rxn
#   - round up: 0.5 U/rxn = 0.1 µL
sw reaction \
  "nuclease-free water ;           ; to 6.00 µL; +" \
  "RNase H buffer      ;       10x ;    0.60 µL; +" \
  "f89 + o224          ;    300 nM ;    3.20 µL; -" \
  "RNase H             ; 5000 U/mL ;    0.10 µL; +" \
  -t "+RNase H" -n 4 |

sw note "mRNA volume and concentration derived from expt #18.  Expts #61,62,65 have a different ratio, but I'm not sure why." |
sw note "100 U/mL RNase H final; need at least 8 U/mL to digest all template, see expt #77." |

sw reaction \
  "nuclease-free water ;          ; to 6.00 µL ; +" \
  "RNase H buffer      ;      10x ;    0.60 µL ; +" \
  "f89 + o224          ;   300 nM ;    3.20 µL ; -" \
  -t "−RNase H" -n 4 |

sw step "Incubate at 37°C for 20 min." |
sw note "https://tinyurl.com/y2deaw9k" |

sw gel urea 11 -v 5 -S |
sw gelgreen -I |
sw laser blue red

