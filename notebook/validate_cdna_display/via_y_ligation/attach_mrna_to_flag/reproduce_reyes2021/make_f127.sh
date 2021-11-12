#!/usr/bin/env bash
set -euo pipefail

# Linker quantity:
# - [Reyes2021] calls for an excess of linker, but I know that only about half 
#   of the mRNA will react.
# - I don't care about reproducing this part of the protocol exactly; I'll get 
#   clean, ligated mRNA after the gel purification step either way.  So it 
#   makes sense to conserve material.

# mRNA quantity:
# - 10 µM f111 is 399 ng/µL
# - I can load 20 µg/lane for gel purification.
# - 20 µg is 50 µL.
# - I'll have to check how much f111 I have available.  I might scale down the 
#   reaction so I odn't have to make more.
sw reaction \
  -s "annealing" \
  "water;; to 98.67 µL" \
  "T4 DNA ligase buffer; 10x; 10 µL" \
  "FLAG mRNA (f111); 10 µM; 50 µL" \
  "[Reyes2021] linker (o243); 100 µM; 5 µL" \
  |

sw step "Incubate as follows:~90°C for 30s~Cool to 25°C at 1°C/s" |

sw step "Add to the reaction:
  ~0.33 µL × 10 U/µL = 3U T4 PNK
  ~2 µL × 10 U/µL = 20U T4 RNA ligase" |

sw step "Incubate at 25°C for 30 min."
