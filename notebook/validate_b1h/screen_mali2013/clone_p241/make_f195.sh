#!/usr/bin/env bash

# Copied and modified from `make_f186.sh` in expt #199.

# Esp3I vs BsmBI:
# - Prefer Esp3I because Quick CIP requires 37°C incubation.
# - Expt #206 seems to show that both enzymes work equivalently.
#
# Esp3I volume:
# - Unit definition: digest 1 µg λ DNA at 37°C in 50 µL.
#
#   1 µg λ DNA   1 µmol λ   14 µmol Esp3I sites   1e6 pmol
#              × ──────── × ─────────────────── × ──────── = 0.467 pmol
#                3.0e7 µg        1 µmol λ          1 µmol
#
#   1 unit is enough to digest 0.467 pmol of Esp3I sites.
#
# - p237:
#
#   5 µg p237    1 µmol    1e6 pmol   2 pmol Esp3I sites   1 U Esp3I
#             × ──────── × ──────── × ────────────────── × ────────── = 14.29 U
#               1.5e6 µg    1 µmol       1 pmol p237       0.467 pmol
#
#   NEB recommends a 10x excess of enzyme (and I found that this was necessary, 
#   at least for my expired enzyme).  Esp3I is 10 U/µL, so this works out to 
#   14.29 µL.
#
# Quick CIP volume:
# - In `make_f186.sh`, I calculated that 6.64 µL of Quick CIP was appropriate 
#   given the MW of p236.
# - p237 has nearly the same MW, and I'm digesting the same amount, so I just 
#   reused the number.

sw future/reaction \
  "water            ;           ;    ; to 250 µL" \
  "p237             ; 383 ng/µL ;    ; 13.05 µL" \
  "rCutSmart buffer ;       10x ; 1x ; " \
  "Esp3I            ;   10 U/µL ;    ; 14.29 µL" \
  "Quick CIP        ;    5 U/µL ;    ; 6.64 µL" \
  -s digestion/phosphatase \
  -i "Split the reaction into 3 tubes with 83.3 µL each." |

sw thermocycler 37/2h 80/20m

