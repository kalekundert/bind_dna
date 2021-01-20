#!/usr/bin/env zsh
set -euo pipefail

# Annealing step:
# - On one hand, don't want to add extra salt to the system.
# - On the other, don't want to dilute the system either.
# - Best would be to anneal in 1x reaction buffer, but I can't do that for 
#   PURExpress
# - Either way, best is to use as little material as possible.

sw zap |
sw cond pick_rnase_oligos_conditions.xlsx |
sw anneal 'probe RNA (f11)' 'target oligos (−,−,o216-o223)' -c 2 -C 10 -n 10 -m 1 -v 2.5 |
sw reaction -t 'RNase H' 'nuclease-free water;;to 10 µL' 'RNase H buffer;10x;1 µL' 'RNase H;5000U/mL;0.2 µL' -n 10 |
sw note "100 U/mL RNase H" |
sw add_probe.sh |
sw smart_mmlv.py 'crude samples + probes' o214 -n 9 -m primer -t 0.5 -T 100 |
sw qpcr 'crude RT products' o214 o215 -a 'optimal Ta' -l 88 -n 11 -v 66 -m primers

# I might want to add more of the RT product to the qPCR reaction if 
# sensitivity is an issue, but I want to start off being conservative
