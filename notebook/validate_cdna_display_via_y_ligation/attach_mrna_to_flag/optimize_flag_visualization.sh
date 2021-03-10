#!/usr/bin/env bash
set -euo pipefail

# Quick experiment to make sure I can visualize FLAG expression.

# Load volumes:
# - Want to try two volumes:
#   - 2.5 µL: What I normally do for these reactions.
#   - 1.0 µL: What the PUREfrex manual recommends with FluoroTect.
sw cond optimize_flag_visualization_cond.xlsx |

sw ivtt f110 -p frex1/lys -n 2 -v 10 -c 2 -C 75 |

# RNase A volume:
# - The PUREfrex 2.0 manual calls for 1 µL 1 mg/mL RNase A per 10 µL reaction.  
# - I don't have any RNase A on hand, but I do have RNase cocktail (Invitrogen 
#   AM2286).  This cocktail contains:
#   - 500 U/mL RNase A
#   - 20,000 U/mL RNase T1
# - The manual for the cocktail says to replace RNase A at equivalent RNase A 
#   concentrations.
# - I can't easily find a way to convert between mass and units, so I'm just 
#   going to use the same proportion as in the PUREfrex manual and hope that 
#   works.
sw step "Split each reaction into 2x 4 µL aliquots." |
sw step "Add 0.4 µL RNase cocktail (Invitrogen AM2286) to the +RNase A aliquots." |
sw step "Incubate at 37°C for 15 min." |

sw gel bolt/ivtt/blue 8


