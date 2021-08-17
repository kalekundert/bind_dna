#!/usr/bin/env bash
set -euo pipefail

# Standard concentrations:
# - [Knight2003] uses 50, 75, 100, 150, 200, 300, and 400 ng.
#   - Assuming 5 µL/lane, that's 10-80 µg/mL
#   - This isn't a very wide range; only 8-fold.
#   - If I were to use these bounds, I'd make a geometric series instead.

sw dilute 'BSA standard' -C 2mg/mL -c 80µg/mL -V 500 |
sw direct 20µL 80µg/mL to 10 7 -m 'BSA standard' |

sw step "Prepare the same dilutions of the concentrated PCV2-Zif268 stock." |

sw gel bolt/mes 14 -v 5 -S |

sw stain sypro-orange |

sw step "Quantify band intensities and calculate PCV2-Zif268 concentration by linear regression."



