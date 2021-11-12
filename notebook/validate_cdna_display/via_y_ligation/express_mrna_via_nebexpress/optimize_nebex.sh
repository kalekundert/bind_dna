#!/usr/bin/env bash
set -euo pipefail

sw zap |

sw serial 2.0 '5000 nM' to 50 7 -0 -m f15 -d 'nuclease-free water' |

# mRNA concentration:
# - In expt #18, the highest concentration I tested was 1600 nM.
#   - 1.6 µL of 10 µM stock in 10 µL reaction.
# - NEBexpress calls for half as much mRNA as PURExpress:
#   - NEBExpress: 1-5 µg in a 50 µL reaction
#   - PURExpress: 1-5 µg in a 25 µL reaction
# - So I'll use 800 nM as my highest concentration for this experiment.
# - This time, my template is only 5 µM.  Fortunately, that's enough.

sw ivtt f15 -p nebexpress -v 5 -n 8 -c 800 -C 5000  |

sw gel bolt 8 -S |
sw laser blue red

