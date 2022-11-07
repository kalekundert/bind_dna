#!/usr/bin/env bash

# This protocol is mostly copied from `make_f186.sh` in expt #199.
# - I recalculated the amount of CIP to used, based on the MW=2.21e6 of p238.

# 2022/11/04:
# - Double the amount of BbsI, because it's been hard getting digests to go 
#   to completion.
# - Update p238 concentration.
# - Increase incubation time.

sw future/reaction \
  "water            ;           ;    ; to 250 µL" \
  "p238             ; 133 ng/µL ;    ; 37.59 µL" \
  "rCutSmart buffer ;       10x ; 1x ; " \
  "BbsI-HF          ;   20 U/µL ;    ; 5.00 µL" \
  "Quick CIP        ;    5 U/µL ;    ; 4.52 µL" \
  -s digestion/phosphatase \
  -i "Split the reaction into 3 tubes with 83.3 µL each." |

# Need to be able to separate the following in gel purification:
#
# - target: 3202
# - 1 cut: 3576
# - circular plasmid
#
# Gonna be very hard.  Really need digest to go to completion, in which case 
# I'm really just purifying away the restriction enzymes.

sw thermocycler 37/2h 65/20m

