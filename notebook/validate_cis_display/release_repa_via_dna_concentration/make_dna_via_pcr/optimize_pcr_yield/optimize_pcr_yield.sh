#!/usr/bin/env bash
set -euo pipefail

# Primer conc:
# - Default primer concentration is 0.5 µM.
# - Want to try 4x and 2x that, i.e. 2 µM and 1 µM.
sw pcr \
  --product f66 \
  --num-reactions 6 \
  --primer-conc 2 \
  --primer-mix-volume 20 \
  --only-primer-mix \
  |

sw serial_dilution \
  10µL 10x / 2 3 \
  --material "primer mix" \
  |

sw pcr \
  --product f66 \
  --num-reactions 3 \
  --reaction-volume 24 \
  --master-mix dna \
  --skip-primer-mix \
  --skip-thermocycler \

sw step "Split each reaction into two 10 µL reactions, one for each of the 
following thermocyler protocols." |

sw pcr \
  --product f66 \
  --only-thermocycler \
  --num-cycles 35 \

sw pcr \
  --product f66 \
  --only-thermocycler \
  --num-cycles 45 \

sw spin_cleanup \

sw step "Measure yield by NanoDrop."
