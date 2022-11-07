#!/usr/bin/env bash
set -euo pipefail

# Setup the digestions separately, to avoid diluting anything.
sw digest p241 BsaI-HFv2 -d 10 |
sw skip -- -1 |
sw digest p243 BsaI-HFv2 -d 10 |

# Minimum concentrations:
# - I want to choose a minimum concentration that I'm confident should resolve 
#   the bands nicely.
#
# - Alejandro Martin calls for 500 ng/band/cm² of the gel's cross-sectional 
#   area, as an empirical rule-of-thumb:
#   https://www.researchgate.net/post/How_much_DNA_can_you_load_per_well_in_a_2_acarbose_gel_without_overloading 
#
#   My wells are 7mm wide and 5 mm tall, for a total of 0.35 cm² per lane.  
#   That would suggest loading 175 ng/lane.
#
# - For the 1kb+ ladder, NEB calls for loading 1 µg of total DNA (120 ng of the 
#   brightest bands) in 6 µL for a 5mm well.  Assuming a well length of 1 mm, 
#   volume of 6 µL would have a height of 1.2 mm, for a total cross-sectional 
#   area of 0.06 cm².  That works out to 2000 ng/band/cm², which exceeds the 
#   above recommendation by 4x (despite being analytical rather than 
#   preparative) and suggests loading 700 ng/lane.
#
# - I think the Martin recommendation is too low, but I wantto go a little 
#   below the NEB recommendation, so I'm going to set my low-end concentration 
#   to 1000 ng/band/cm², or 350 ng/lane.
#
# - To get a height of 5mm, I'll need to load a total volume of 5mm × 7mm × 
#   1.5mm = 52.5 µL in each well.  Accounting for the 6x loading dye, that 
#   works out to 43.75 µL of DNA.
#
# - 350 ng / 43.75 µL = 8 ng/µL

# High concentrations:
# - This is basically determined by how concentrated the midiprepped DNA is.
# - I also want to try lyophilizing the DNA, to see if lyophilization 
#   itself has any ill-effect on the gel.
#
# - f202:
#   - 10 µg / 181.92 µL = 54.97 ng/µL
#
# - f202 after lyophilization:
#   - 10 µg × [(181.92 - 70.72) / 181.92] = 6.11 µg
#   - 6.46 µg / ≈60 µL = ≈100 ng/µL
#
# - f203:
#   - 10 µg / 50.36 µL = 198.57 ng/µL

# Dilutions:
# - The wide gels have 13 lanes.
#   - 4 dilutions with no gaps
#   - 3 dilutions with gaps
#
# - Care more about getting information than getting good purifications, so 
#   I would prefer doing 4 dilutions.  However, that ends up taking too much 
#   material, so I'll stick with 3.

sw serial '55 ng/µL' to 8 -n 3 -v 43.75 -m f202 -d water |
sw lyophilize -v 60 |
sw sub 'the sample.s.' 'the remaining f202 digestion' |
sw serial '100 ng/µL' to 8 -n 3 -v 43.75 -m 'f202 (lyo)' -d water |
sw serial '200 ng/µL' to 8 -n 3 -v 43.75 -m f203 -d water |

# Preset:
# - Relevant fragment sizes:
#
#   Fragment  Keep?  Size (kb)
#   --------------------------
#   f202        +          2.8
#               -          1.0
#   f203        +          4.1
#               -          0.2
#
# - The 4kb preset is probably about right.
# - I don't think I've actually optimized the 4 kb preset yet, so this will 
#   also be a chance to improve those parameters.
#
#   Update: 1%, 100V, 90m worked well.  4kb was right about in the middle of 
#   the gel.  This was the big gel, though, so it will be different for the 
#   smaller gels.

sw gel agarose/purify/4kb 'f202, f202 (lyo), f203' -n 9 -V 52.5 |
sw qiaex_ii '2.8 kb (f202) or 4.1 kb (f203)'
