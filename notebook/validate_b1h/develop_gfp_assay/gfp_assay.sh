#!/usr/bin/env bash
set -euo pipefail

# Would prefer flow cytometry, but don't think that's an option.

sw step "Grow 1 mL overnight cultures of sz221-sz228 in LB+Carb at 37°C for 16h." |
sw note "The small culture size is meant to reduce the chance of cross-contamination when using 24-well blocks." |

# Growth time:
# - There needs to be enough time to:
#   - Express rpoZ-Zif268 (induced with IPTG).
#   - GFP expression level to reach a steady state.
#
# - RFP is constitutively expressed, so it will already be at steady state.  If 
#   GFP is not also at steady state, the ratio of the two will not be 
#   reproducible.
#
# - Probably the best way of "degrading" GFP/RFP is for the cells to divide.  
#   The more divisions, the less the RFP from the overnight cultures will 
#   affect the assay.  So I'll want relatively thick cultures.

# - I probably want the cells to still be in log phase, because in stationary 
#   phase the rate of RFP/GFP "degradation" by dilution will go way down, and a 
#   new steady state will be reached.  Plus the cells will be more stressed, 
#   and in gnereal it seems like a good idea to stick to log phase.
#
# - I might want to measure several timepoints for the first experiment.
sw step "Inoculate day cultures:~1 mL LB + Carb + 10 µM IPTG~10 µL overnight culture~Reserve some media for making serial dilutions later." |
sw step "Grow day cultures at 37°C to late log phase." |
sw note "I don't know exactly what OD represents late log phase.  I'll have to do an experiment to figure that out." |

sw step "Wash 100 µL cells from each culture with PBS:~Repeat 2x:~~Spin 3500g, 3 min, 4°C~~Resuspend in 100 µL PBS" |
sw note "This is necessary because the yeast extract in LB is highly fluorescent: https://tinyurl.com/3nvmpu94"

# Is it necessary to mesure dilutions?
#
# - It's basically just replicates, but I think it's better because it also 
#   gives me some information on dynamic range.
#
# Number of dilutions:
#
# - 4 dilutions would allow me to fit two strains per column.  (Note that I 
#   don't need to worry about evaporation because this is an endpoint 
#   measurement.)
# 
# Volume:
# - I'm only going to use 50 µL of each culture, but 100 µL is convenient for 
#   10x serial dilutions.
#
# Dilution factor:
# - My goal is to determine a linear range, so I feel like a 100x difference 
#   between most dilute and most concentrated is reasonable.
sw serial 1x to 0.01x -n 4 -v 100µL -m "day culture" -d "LB+Carb+IPTG" |

sw step "Prepare a plate for fluorescence measurement:~Use a black, opaque-bottom plate.~Add:~~150 µL PBS~~50 µL washed/diluted cells~Pipet to mix."

sw step "Measure the following using a plate reader:~mWasabi: 493/509 nm ex/em~mCherry: 587/610 nm ex/em"
