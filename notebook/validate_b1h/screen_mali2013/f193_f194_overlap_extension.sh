#!/usr/bin/env bash
set -euo pipefail

# Reaction volume:
# - [Hilgarth2020] calls for this reaction to be 50 µL, but only 4 µL will be 
#   used in the next step.
# - In the future, I should make this reaction be 5 µL.  Maybe 10 µL if I have 
#   plans to run a gel or something.
sw future/reaction -f f193_f194_overlap_extension.xlsx |

# Touchdown PCR protocol:
# - This is a hybrid of [Hilgarth2020] and the Q5 protocol recommened by NEB: 
#   https://international.neb.com/protocols/2012/12/07/protocol-for-q5-high-fidelity-2x-master-mix-m0492
sw thermocycler 98/30 "14x 98/10 72/20 72/1m30" 72/2m |
sw sub 20s "20s; -0.5°C each cycle" |

sw future/pcr "previous reaction,o396,o397" \
  --reaction-volume 50 \
  --template-volume 4 \
  --primer-conc 1 \
  --skip-thermocycler |

# Annealing temperatures:
# - [Hilgarth2020] calls for 67-59°C.
# - My primers have a Tm of 66°C, though, so I want to start at 72°C (basically 
#   the highest possible temperature) and end a few degrees below 66°C.  I can 
#   do this by simply shifting the 67-59°C window up by 5°C.
sw thermocycler 98/30 "17x 98/10 72/20 72/1m30" "23x 98/10 64/20 72/1m30" 72/2m |
sw sub "72°C for 20s" "72°C for 20s; -0.5°C each cycle"
