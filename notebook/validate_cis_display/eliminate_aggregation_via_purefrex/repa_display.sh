#!/usr/bin/env bash
set -euo pipefail

# Time: 5h15
# - setup reaction: 30m
# - incubate: 2h
# - setup gel: 15m
# - run gel: 2h
# - image: 30m

sw conditions repa_display_cond.xlsx |

# Template concentration
# - As described in expt #97, I want the template to be as concentrated as 
#   possible.
# - f50-f66 were diluted to 75 nM, which is about 100 ng/µL.  This could 
#   probably be higher, but it's not unreasonably dilute.
# - 28.5 nM is the highest template concentration I can achieve in a PUREfex 
#   1.0 reaction with the 75 nM stocks.

# Reaction volume:
# - I want the biggest volume possible without using more than one aliquot of 
#   PUREfrex.
# - 1 aliquot: 12.1 µL solution I
# - 3 µL reactions use just a little less than that, so there's some extra.
sw ivtt f50 f61 f62 f63 f64 f65 f66 -p frex1 -v 3.0 -c 28.5 -C 75 |

sw gel natbr 14

