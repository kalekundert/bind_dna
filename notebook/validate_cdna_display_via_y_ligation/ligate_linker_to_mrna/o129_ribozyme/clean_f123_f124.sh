#!/usr/bin/env bash
set -euo pipefail

# DNase treatment:
# - I haven't done this before, but it occurs to me that I don't want any 
#   template DNA in the protein expression reaction.  Without a DNase treatment 
#   step, I'm probably adding a lot of it.

sw step "Resuspend each vial of lyophilized DNase I (Zymo E1009-A, 250 U) in 275 µL nuclease-rfee water.  Mix by gentle inversion.  Store frozen aliquots." |

sw reaction \
  -s "DNase" \
  "water                ;         ; to 200 µL ; +" \
  "DNA digestion buffer ;     10x ;     20 µL ; +" \
  "DNase I              ;  1 U/µL ;     20 µL ; +" \
  "IVT reaction         ;  1 U/µL ;     20 µL ; -" \
  -n 2 |

sw note "From Zymo Clean & Concentrator manual" |

sw step "Incubate at room temperature for 15 min." |

sw spin_cleanup zymo/rna-clean-conc/100
