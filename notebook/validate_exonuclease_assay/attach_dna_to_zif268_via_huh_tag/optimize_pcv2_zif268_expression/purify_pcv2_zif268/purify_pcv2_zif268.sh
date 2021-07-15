#!/usr/bin/env bash
set -euo pipefail

# Culture volume:
# - In "The QIAexpressionist", a "standard" culture is 100 mL.  That seems 
#   reasonable to me.
sw grow_pcv2_zif268_100mL.txt |
sw nta/lyse_native |
sw nta/purify_native

sw gel bolt/mes NI,I,CL,FT,W1,W2,E1,E2,E3,E4 -v 5 --mix-volume 20

# Future steps:
# - Measure protein concentration:
#   - Nanodrop
#   - Bradford
#
# - Measure protein activity
#   - EMSA
#
# - Determine best storage conditions:
#   - Flash freeze and store at −80°C
#   - −20°C in 50% glycerol
#
# This experiment is just about expressing the protein, though.
#
# I should've included steps to concentrate the protein (Amicon) and transfer 
# in to storage buffer.
