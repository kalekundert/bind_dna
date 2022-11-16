#!/usr/bin/env zsh
set -euo pipefail
zmodload zsh/mathfunc

# Digest with DpnI

# Copied and modified from `make_f196.sh` in expt 202.

# f201 quantity:
# - I expect to get ≈3 µg f210 (50 µL at ≈60 ng/µL) from the PCR reaction 
#   followed by the spin cleanup.
# - Update the variables below based on the actual yield.

f210_ug=4
f210_ng_uL=95

# Esp3I volume:
# - Using equations from `make_f196.sh`, expt #202.
# 
# - I can easily monitor the progress of this reaction and purify the complete 
#   product, so I use a 3.33x excess of enzyme instead of a 10x excess.  This 
#   also helps because otherwise (due to the relatively small size of f210) a 
#   very large volume of enzyme would be needed.

f210_uL=$((f210_ug * 1e3 / f210_ng_uL))
f210_mw=776552.1515

esp3i_sites_per_unit=$((14 * 1e6 / 3e7))
esp3i_uL=$((f210_ug * 1e6 * 2 / f210_mw / esp3i_sites_per_unit / 3))

# Reaction volume:
# - NEB recommends 50 µL per 1 µg template.

rxn_uL=$((25 * ceil(2 * (f210_ug - 1e-3))))
dpni_uL=$((rxn_uL / 100.))

sw future/reaction \
  "water            ;                    ;    ; to $rxn_uL µL" \
  "f210             ;  $f210_ng_uL ng/µL ;    ; $f210_uL µL" \
  "rCutSmart buffer ;                10x ; 1x ; " \
  "Esp3I            ;            10 U/µL ;    ; $esp3i_uL µL" \
  "DpnI             ;            20 U/µL ;    ; $dpni_uL µL" \
  -s digestion |

sw thermocycler 37/2h 80/20m

# Purification:
# - This could be PAGE purified...
