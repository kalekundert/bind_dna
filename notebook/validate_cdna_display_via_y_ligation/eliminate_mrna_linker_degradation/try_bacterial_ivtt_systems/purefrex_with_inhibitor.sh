#!/usr/bin/env bash
set -euo pipefail

IVTT_FLAGS='-v5 -r -t0'

sw zap |

sw cond purefrex_with_inhibitor_cond.xlsx |

# sw reaction \
#   -t "PUREfrex master mix/" \
#   'solution I   ;;    2.50 µL' \
#   'solution II  ;;    0.25 µL' \
#   'solution III ;;    0.50 µL' \
#   -V '4.4*1.1' |

# sw reaction \
#   -t "+inhibitor PUREfrex" \
#   'PUREFrex master mix      ;         ; 3.25 µL; +' \
#   'RNase inhibitor, murine  ; 40 U/µL ; 0.20 µL; +' \
#   'f89                      ; 1000 nM ; 1.55 µL; -' \
#   -n 2 |

# sw reaction \
#   -t "−inhibitor PUREfrex" \
#   'PUREFrex master mix      ;         ; 3.25 µL; +' \
#   'water                    ; 40 U/µL ; 0.20 µL; +' \
#   'f89                      ; 1000 nM ; 1.55 µL; -' \
#   -n 2 |

# Decided to just go ahead with pipetting small volumes of the inhibitor.  
# Making three master mixes doesn't dramatically increase the volumes being 
# pipetted, and I don't really care if the +/- f89 reactions get different 
# amounts of RNase inhibitor.  (The -f89 +inhibitor reaction is really not very 
# interesting.  I'm only including it for completeness.)

sw ivtt f89 -p frex/gfp -n4 -v5 -r -t30 -a 'RNase inhibitor, murine;40 U/µL;0.2 µL;-' |

sw gel bolt 4 -S |
sw laser blue red

