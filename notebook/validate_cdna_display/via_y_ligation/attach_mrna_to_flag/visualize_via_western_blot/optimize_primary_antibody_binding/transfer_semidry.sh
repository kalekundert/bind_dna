#!/usr/bin/env bash
set -euo pipefail

TRANS_BLOT_QUICK_REF=https://tinyurl.com/6v2cm88
TRANS_BLOT_MANUAL=https://tinyurl.com/3y8ent47

LC2002_MANUAL=https://tinyurl.com/hadb3v37
PVDF_MANUAL=https://tinyurl.com/8dx8ndys

# Membrane/filter products:
# - PVDF/Filter Paper Sandwich, 0.2 µm, 8.3 x 7.3 cm (Invitrogen LC2002)
#   - Pre-cut
#   - Prefer to use invitrogen for electrophoresis stuff
#   - $238
#
# - Low-Fluorescence PVDF Transfer Membranes, 0.2 μm 7 x 8.4 cm (Invitrogen 
#   22860)
#   - Membrane only
#   - $128
#   - I already have filter paper; don't need to pay $110 for it.

sw step "Equilibrate gel in transfer buffer for 10-15 min." |
sw note $PVDF_MANUAL |

sw step "Cut the membrane to the size of the gel (if necessary) and notch one corner." |
sw note $PVDF_MANUAL |

sw step "Wet the membrane in 100% ethanol for 15 s, then rinse with water for 2 min with gentle shaking." |
sw note $PVDF_MANUAL |

sw step "Soak the membrane in transfer buffer for 5 min." |
sw note $PVDF_MANUAL |

sw step "Cut the filter paper to size and briefly soak in transfer buffer."

se step "Assemble the transfer sandwich (start from the bottom):

~ filter paper
~ platinum anode" |

