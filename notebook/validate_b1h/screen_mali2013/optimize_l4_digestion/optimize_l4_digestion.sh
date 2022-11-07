#!/usr/bin/env bash
set -euo pipefail

# - 150 (6*25) µL PCR
# - PCR cleanup
# - Digest with 3 µL BbsI
# - Take aliquots at given time points, and heat denature.

sw pcr l2,sr4,sr5 -v 175 -y 10 -a 61 -x 15 |
sw step "Label the product: l4" |
sw spin_cleanup qiagen -s 175 -v 70 |

sw step "Reserve 10 µL as the −enzyme (0 min) control." |

sw digest l4 BbsI-HF -d 6 -D 100 |
sw skip |
sw step "Incubate at 37°C.  At the same time, pre-warm another thermocycler block to 65°C.  At the times indicated below, remove an 11.67 µL aliquot from the digestion reaction and incubate it at 65°C for 20 min.~5m~15m~30m~1h~2h" |

sw gel tbe f191 -n6 -v 11.67 -V 14.5875
