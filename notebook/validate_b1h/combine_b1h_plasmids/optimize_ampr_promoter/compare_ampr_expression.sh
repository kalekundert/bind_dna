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
#   (low expression) or s16 (high expression).
#
# - Maybe I should also grow the cells +/− antibiotic.  That runs the risk that 
#   the cells will lose the plasmid during the experiment.  I could use plain 
#   s3 as the no-antibiotic control.

sw b1h/od_kinetic_assay s13 s16 s{23..27}
