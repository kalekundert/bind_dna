#!/usr/bin/env bash
set -euo pipefail

# - Prepare stock of 20x enzyme + 2x buffer
# - Dilute into 2x buffer
# - Prepare 2x stock of DNA
# - Mix enzyme with DNA
#   - 15 µL reactions: want to stay below 10% enzyme, will need ≈1.2 µL enzyme

# Enzyme volume:
# - pmol BsmBI site/unit: k = 14 * 1e6 / 3e7 = 0.467 pmol
# - 1x units: u = 15 * 2 * 1e3 / (50545.322 * k) = 1.27 U
# - 10x units per 8 µL: 10 * u / 8 = 1.5898 U
sw future/reaction \
  'water; ; ; to 11 µL' \
  'buffer; 10x; 2x;' \
  'enzyme; 10 U/µL; 1.5898 U/µL;' \
  -C 'enzyme, buffer' \
  -c 'BsmBI-v2, NEBuffer r3.1' \
  -c 'Esp3I, rCutSmart' \
  -S 'Setup {n:# enzyme mix/es}:' |

sw future/reaction \
  'water; ; ; to 10 µL' \
  'buffer; 10x; 2x;' \
  'enzyme mix; ; ; 1 µL' \
  -S 'Make a 10x dilution of each enzyme mix:' |

# Reaction volume:
# - The most important consideration is to keep the volume of enzyme below 10% 
#   of the reaction, as recommended by NEB.
# - The enzyme volume for the 10x reaction works out to 1.27 µL, so the total 
#   reaction volume need to be at least 12.7 µL.
# - I rounded up to 16 µL because that leads to nice round numbers when setting 
#   up the gel (4 µL loading due, 20 µL load).
sw future/reaction \
  'water; ; ; to 16 µL' \
  'f206; 15 ng/µL; ; 1 µL' \
  'enzyme mix; 2x; 1x;' \
  -C 'enzyme mix' \
  -c 'water' \
  -c '10x BsmBI' \
  -c '1x BsmBI' \
  -c '10x Esp3I' \
  -c '1x Esp3I' \
  -s 'digestion' |

sw step "Incubate as follows:~BsmBI:~~55°C for 1h~~80°C for 20m~Esp3I:~~37°C for 1h~~65°C for 20m" |

sw gel page/tbe/f206 5 -V 20
