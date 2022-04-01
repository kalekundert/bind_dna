#!/usr/bin/env bash
set -euo pipefail

# Do I need to clone AAA targets?
#
# - I can't think of any reason that AmpR expression should have different 
#   effects for the different targets, but maybe there's some weird pressure on 
#   the cells.
#
# - In expt 158, I can see a difference in growth just in the selective 
#   condition.
#
# - Let me start with just the TGG targets.  That should be informative, and 
#   that'll be what I clone first anyways.  I'll probably follow up with the 
#   AAA targets for select constructs that look promising.
#
# Controls:
#
# - The goal of the experiment is to see if the constructs behave more like s13 
#   (low expression) or s16 (high expresison).
#
# Assay:
#
# - Right now I'm planning to use the colony assay, since that's the one I have 
#   most established.  If I decide that the OD600 is reliable by the time I 
#   have these constructs, I might do that instead.

sw b1h/plate_assay s13 s16 s{23..27}
