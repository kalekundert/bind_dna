#!/usr/bin/env bash
set -euo pipefail

# Ni-NTA purification:
# - To measure expression, it would be most efficient to run only the following 
#   lanes for each condition:
#   - non-induced (NI)
#   - induced (I)
#   - cleared lysate (CL)
#
# - However, these lanes didn't work in the last experiment (#113).
# - So until I know that I can get the information I need from these lanes, I 
#   should continue to do the Ni-NTA purification.
# - It may also be informative to see if the amount of purified protein differs 
#   from the amount of protein in the clarified lysate.
# - In the meantime, I'll be a little limited in the number of conditions I can 
#   test (or I'll just have to run lots of gels).
#
# Rhamnose concentrations:
# - NEB recommends testing 7 rhamnose concetrations.
# - I want to start with fewer, since I'll be doing the Ni-NTA purification for 
#   each:
#
#   4 concentrations:
#   - won't fit on two gels
#
#   3 concentrations:
#   - will fit on two gels
#   - I'll do 0, 100, 250 to focus on the low concentrations.
#   - My thinking is that:
#     - The higher concentrations are for more difficult proteins, and my 
#       protein didn't seem that difficult in the first experiment.
#     - I want the lowest concentration that works, since that will be the 
#       highest expression.

sw cond optimize_inducer_cond.xlsx |
sw step "Thaw rhamnose and make 1 mL aliquots." |
sw grow_miniprep_rhamnose.txt |
sw ni_nta/miniprep_native |
sw gel bolt/mes 32
