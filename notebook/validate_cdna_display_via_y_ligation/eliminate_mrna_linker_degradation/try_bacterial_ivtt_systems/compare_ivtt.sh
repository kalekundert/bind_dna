#!/usr/bin/env bash
set -euo pipefail

# Volume:
# - I only need 2.5 µL to run the gel.
# - But I'm gonna do 10 µL because I'm going to have to thaw an entire IVTT 
#   aliquot for each reaction anyways, and with no master mix, I want to make 
#   the volumes a little easier to pipet.
IVTT_FLAGS="-v 10 -r -C 1000 -t 0"

sw zap |
sw cond compare_ivtt_cond.xlsx |
sw ivtt f89 -p purex/gfp $IVTT_FLAGS |
sw ivtt f89 -p frex/gfp  $IVTT_FLAGS |
sw ivtt f89 -p nebex/gfp $IVTT_FLAGS |
sw ivtt f89 -p s30/gfp   $IVTT_FLAGS |
sw step "Incubate all 4 IVTT reactions at 37°C for 30 min." |
sw step "Prepare 200 nM f89 as a control." |
sw gel bolt 5 -S |
sw laser blue red
