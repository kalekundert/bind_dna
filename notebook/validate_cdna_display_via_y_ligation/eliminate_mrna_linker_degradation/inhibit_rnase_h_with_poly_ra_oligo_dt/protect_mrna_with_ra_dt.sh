#!/usr/bin/env bash
set -euo pipefail

sw zap |

# Test PURExpress or RNase H?
# - I know the basic idea will work with RNase H, so I'm really interested in 
#   the amount needed to eliminate activity in PURExpress.
# - It'll be nice to know the relative concentrations of RNase H in the 
#   different expression systems.  So I'm glad I spent the money on the Lee2018 
#   RNase H assay reagents.
# - Plus, an advantage of using pure RNase H is that I can stain using 
#   GelGreen, but that probably won't work very well with the poly-rA/oligo-dT 
#   reagent.

sw cond protect_mrna_with_ra_dt_cond.xlsx |

sw step "Dilute lyophilized poly-rA/oligo-dT (Midland P-4012) to 80 A260 U/mL:~5 A260 U poly-rA/oligo-dT~62.5 µL nuclease-free water"

sw serial \
  10 "8 U/mL" / 2 5 \
  -0 \
  -m "poly-rA/oligo-dT" \
  -d "nuclease-free water" \
  |

sw ivtt \
  f89 \
  -p 'purex/gfp' \
  -v 5 \
  -n 6 \
  -mr \
  -a 'poly-rA/oligo-dT;0-80 A260 U/mL;0.5 µL;-' \
  -t '30 min' \
  |
  
# Should I include the coupling reaction?  It's not really what I'm trying to 
# test here, but it is the end goal...
#
# I think I'll skip this for now.  If I do the coupling step, I'd also want to 
# exclude it just to make sure that it isn't causing any problems itself, and 
# that would mean a bunch more conditions.  Better to do this experiment to 
# work out a good poly-dA/oligo-dT concentration, then to do another experiment 
# to see if the coupling reaction works.

# Prepare controls:
# - 160 nM o194 (stock: 1 µM)
# - 160 nM f89  (stock: 1 µM)
sw reaction \
  -t "−mRNA control/s" \
  "nuclease-free water  ;              ;   to 5.0 µL ; +" \
  "o129                 ;         1 µM ;      0.8 µL ; +" \
  "poly-rA/oligo-dT     ; 80 A260 U/mL ;      0.5 µL ; -" \
  -n 2 \
  |

sw reaction \
  -t "−PURExpress control/s" \
  "nuclease-free water  ;              ;   to 5.0 µL ; +" \
  "f89                  ;         1 µM ;      0.8 µL ; +" \
  "poly-rA/oligo-dT     ; 80 A260 U/mL ;      0.5 µL ; -" \
  -n 2 \
  |

sw gel bolt/o194 10 -S |
sw laser blue red
