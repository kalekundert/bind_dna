#!/usr/bin/env bash
set -euo pipefail

# Strain:
# - Haven't picked a strain yet.
# - My intention is to use whichever of the AmpR promoter strains grows the 
#   best.
sw step "Grow 3 mL szXXX in LB+Carb at 37°C for 16h" |

# Media:
# - I (want to) grow the cells in NM+His+Ura for my auxotrophy assays.
# - I grow the cells in LB for my RNAseq/GFP assays.
# - I should test both.
sw step "Inoculate 10 mL LB+Carb and NM+His+Ura+Carb with 100 µL saturated overnight culture.  Split each culture between 2 tubes." |

sw step "Incubate at 37°C for 1h" |
sw serial 1 / 10 -n 6 -v 1mL -m "day culture" -d "LB/NM" |
sw step "Spot 5 µL of each dilution in triplicate on a pre-warmed LB+Carb plate." |
sw step "Measure A600 for each dilution using:~The NanoDrop in cuvette mode~The nanodrop in pedestal mode~The plate reader" |
sw step "Repeat the above measurement steps every 1h until the cultures seem mostly saturated (≈6h)."
