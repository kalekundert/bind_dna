#!/usr/bin/env bash

# Copied and modified from `make_f195.sh`.

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
#   2.5 µg p237     1e6 pmol     2 pmol Esp3I sites   1 U Esp3I
#               × ──────────── × ────────────────── × ────────── = 4.94 U
#                 2169820.9 µg      1 pmol p237       0.467 pmol
#
#   NEB recommends a 10x excess of enzyme (and I found that this was necessary, 
#   at least for my expired enzyme).  Esp3I is 10 U/µL, so this works out to 
#   4.94 µL.

sw future/reaction \
  "water            ;           ;    ; to 125 µL" \
  "f209             ;  60 ng/µL ;    ; 41.67 µL" \
  "rCutSmart buffer ;       10x ; 1x ; " \
  "Esp3I            ;   10 U/µL ;    ; 4.94 µL" \
  "DpnI             ;   20 U/µL ;    ; 1.25 µL" \
  -s digestion \
  -i "Split the reaction into 2 tubes with 62.5 µL each." |

sw thermocycler 37/2h 80/20m

