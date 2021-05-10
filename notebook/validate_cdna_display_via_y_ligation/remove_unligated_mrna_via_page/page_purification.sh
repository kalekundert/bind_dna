#!/usr/bin/env bash
set -euo pipefail

cat <<EOF
The problem is that my RNA is so dilute after the ligation reaction that it's 
barely concentrated enough to matter.

Maybe I should use dialysis protocol
EOF

# Template
# - Candidates:
#   - FLAG (f115 = f111 + o239):
#     - not retained by spin filters
#     - should give the biggest separation
#   - Zif268 (f119 = f11 + o129):
#     - what I really care about
#   - mWasabi (f85?):
#     - I do already have this gel
# - I'm inclined to use Zif268 and FLAG.

# I already have three aliquots of f111, that might be enough.
# I have lots of f11.
#sw make f111 |

# mRNA volume:
# - Let's assume I want to end up with 5 µL ligated mRNA stock.
# - IVTT reactions require 1 µM mRNA.
# - Reaching that requires a stock concentration of ≈2 µM.
#   - 5 µL @ 2 µM = 10 pmol
# - Only about half the mRNA is ligated.
#   - 20 pmol
# - I'll only recover 10% of the ligated DNA.
#   - 200 pmol = 20 µL @ 10 µM
#   - f11: 25 µg
#   - f111: 8 µg
sw cdna/anneal f111,f11 o239,o119 -r 5 |
sw cdna/ligate 11.11 4.5 -n 2 |

# - Start with pre-cast gel.
# - Ask: How much can I load in a single lane?
#   - Ideally: 200 pmol / 12 lanes ≈ 15 pmol/lane
#     - f11: 2 µg
#     - f111: 0.6 µg
#   - 10 µL / lane = 1.5 µM
# - Lowest concentration to test:
#   - 200 ng / lane is usual recommendation
#     - f11: 10x dilution
#     - f111: 3x dilution
#   - 
#   - Max: 12 lanes → 1/12 reaction (≈10 µg)
#   - serial dilution by 2 → 2**12 ≈ 1000 → 10 ng
#
# - Ligation reaction:
#   - 450 nM mRNA
#     - 
#
sw gel urea/o194 10

# Purification protocol: follow [Petrov2013]_
# 
# Denaturing gels are pre-run because the tempature helps with denaturing...
#
# Calls for 1.5 mm spacers.  Keep that in mind when ordering.

#sw step "UV shadowing~In dark room" |

# - might buy a plastic tray for this, so I can cut directly in laser scanner
# - if shadowing works, this might not be necessary
# - but even so, the laser would probably be less damaging to the RNA
#sw "cut band"

# sw <<EOF
# - cut band
# - crush
# - add two volumes elution buffer
# - incubate 3h at RT or o/n at 4°C
# - centrifuge 3x

