#!/usr/bin/env bash
set -euo pipefail

# I've read that supercoiled plasmid does not amplify consistently and should 
# not be used for qPCR, so linearize the plasmid here.
#sw digest p49 XmnI -D 389 |

# The DNA is 100 ng/µL after the digest (1 µg DNA in 10 µL).
#sw step "Dilute DNA from 100 ng/µL to 20 pg/µL." |
sw qpcr f108 o115 o116 -a 60 -g 10 -l 100 -n 8 -m dna,primers -T '20 pg/µL' |
sw note "The amount of template doesn't really matter in this reaction; it just needs to be reasonable.  It doesn't have to mimic an RT reaction."
