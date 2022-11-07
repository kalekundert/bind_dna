#!/usr/bin/env bash

# DNA amount:
# - 5 µg seems like a good amount to aim for:
#   - I've seen recommendations around this for restriction cloning
#   - It's a lot, but not too much to run on a gel
#   - It's around the amount of material you'd have in a miniprep.
#
# - Right now, I have 50 µL of 200 ng/µL p236.
#   - Yield isn't super important here.  Worst case I use less DNA and get 
#     fewer transformants, but I'll get 100% library coverage regardless.
#
# Digestion reaction
# - I copied most of the values from:
#
#     sw digest p236 NotI-HF,HindIII-HF -d5 -D200
#
# Phosphatase:
# - Important to treat the backbone with phosphatase, to avoid re-ligation.
# - The digest protocol doesn't support this currently, so that's why I'm 
#   writing this as a custom protocol.
#
# - NEB calls for 1 µL Quick CIP for each 1 pmol of DNA ends.
#
#   - 5 µg 236 × (1 µmol / 1,505,770 µg) × (1e6 pmol / 1 µmol) × (2 pmol ends / 
#     1 pmol) × (1 µL CIP / 1 pmol ends) = 6.64 µL CIP
#

# Reaction volume:
# - NEB comments that "enzyme volume should not exceed 10% of the total 
#   reaction volume to prevent star activity due to excess glycerol".
# 
# - NEB also comments that "a 50 µl reaction volume is recommended for 
#   digestion of 1 µg of substrate".
#
# - Based on these recommendations, I'm going to scale the reaction to 250 µL.
#
# - I'll probably lyophilize to 50 µL before doing the gel purification.

# Enzymes:
# - With p236 as the substrate:
#   - Use 25 µL, 200 ng/µL
#   - Don't use DpnI

sw future/reaction \
  "water            ;           ;    ; to 250 µL" \
  "f207             ;           ;    ; 50 µL" \
  "rCutSmart buffer ;       10x ; 1x ; " \
  "NotI-HF          ;   20 U/µL ;    ; 2.50 µL" \
  "HindIII-HF       ;   20 U/µL ;    ; 2.50 µL" \
  "DpnI             ;   20 U/µL ;    ; 2.50 µL" \
  "Quick CIP        ;    5 U/µL ;    ; 6.64 µL" \
  -s digestion/phosphatase \
  -i "Split the reaction into 3 tubes with 83.3 µL each." |

sw thermocycler 37/15m 80/20m

