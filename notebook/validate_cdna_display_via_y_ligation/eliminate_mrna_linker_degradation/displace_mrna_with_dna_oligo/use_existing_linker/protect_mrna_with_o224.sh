#!/usr/bin/env bash
set -euo pipefail

# Outline
# - Make mRNA
# - Add toehold
# - incubate
# - treat with RNase H
# - PAGE

# To-do
# - Scale RNase H reactions down to 6 µL
#   - Just enough to run gel
#   - Don't waste mRNA

sw zap |
sw cond \
  protect_mrna_with_o224_cond.xlsx \
  -k "dna/dna oligo values are molar ratios relative to mRNA+linker" |

# If I don't have any f89; I might have some.
sw cdna/make_mrna f85 o129 -A |

sw step "Label the product: f89" |

sw serial 9 '100µM' / 10 3 -0 -m o224 |

sw reaction \
  "PBS  ;     10x ; 0.30 µL ; +" \
  "f89  ;    1 µM ; 1.35 µL ; +" \
  "o224 ; 100-0µM ; 1.35 µL ; -" \
  -t "annealing" -n 4 |

sw step "Incubate at 95°C for 2 min, then let cool at room temperature." |

# 2 master mixes, +/- rnase H
# - Same volume/mRNA amount as expression reactions
# - excess of RNase H
sw reaction \
  "nuclease-free water;         ; to 6.00 µL; +" \
  "RNase H buffer     ;      10x;    0.60 µL; +" \
  "f89 + o224         ;     1 µM;    1.00 µL; -" \
  "RNase H            ; 5000U/mL;    0.12 µL; +" \
  -t "+RNase H" -n 4 |

sw note "mRNA volume and concentration derived from expt #18.  Expts #61,62,65 have a different ratio, but I'm not sure why." |
sw note "100 U/mL RNase H final; need at least 8 U/mL to digest all template, see expt #77." |

sw reaction \
  "nuclease-free water ;          ; to 6.00 µL ; +" \
  "RNase H buffer      ;      10x ;    0.60 µL ; +" \
  "f89 + o224          ;     1 µM ;    1.00 µL ; -" \
  -t "−RNase H" -n 4 |

sw step "Incubate at 37°C for 20 min." |
sw note "https://tinyurl.com/y2deaw9k"

sw gel urea 11 -v 5 -S |
sw gelgreen -I |
sw laser blue red

