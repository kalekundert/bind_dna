#!/usr/bin/env zsh
set -euo pipefail
zmodload zsh/mathfunc

# Copied and modified from `make_f195.sh` in expt #202.

# p239 volume:
# - Ideally I'd use 5 µg, which is the amount I've een using for all of my 
#   digests so far.  But this is a low copy plasmid, and I've had trouble 
#   getting that much.  So instead I just use everything I have.
#
# - Update the variables below, and the amount of enzyme will be updated 
#   accordingly.

p239_uL=35
p239_ng_uL=58

# Esp3I volume:
# - Using equations from `make_f195.sh`, expt #202.

p239_ug=$((p239_uL * p239_ng_uL / 1e3))
p239_mw=$((2.12e6))

esp3i_sites_per_unit=$((14 * 1e6 / 3e7))
esp3i_uL=$((p239_ug * 1e6 * 2 / p239_mw / esp3i_sites_per_unit))

# Quick CIP volume:
# - Using equations from `make_f186.sh`, expt #199.

cip_uL=$((p239_ug * 1e6 * 2 / p239_mw))

# Reaction volume:
# - NEB recommends 50 µL per 1 µg template.

rxn_uL=$((25 * ceil(2 * (p239_ug - 1e-3))))

sw future/reaction \
  "water            ;                    ;    ; to $rxn_uL µL" \
  "p239             ;  $p239_ng_uL ng/µL ;    ; $p239_uL µL" \
  "rCutSmart buffer ;                10x ; 1x ; " \
  "Esp3I            ;            10 U/µL ;    ; $esp3i_uL µL" \
  "Quick CIP        ;             5 U/µL ;    ; $cip_uL µL" \
  -s digestion/phosphatase |

sw thermocycler 37/2h 80/20m |

sw step "Run an E-Gel to make sure that the digestion is complete, and has the expected bands.~The expected bands are 3.4 kb and 20 bp.~The Esp3I sites were in the insert of the previous reaction."

