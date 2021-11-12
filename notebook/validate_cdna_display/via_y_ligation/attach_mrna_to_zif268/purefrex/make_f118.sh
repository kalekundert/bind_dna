#!/usr/bin/env bash
set -euo pipefail

# Linker oligo:
# - o128: Spacer18
# - o129: poly-A
#
# - It looks like I used o93 previously (expt #105), but that doesn't make 
#   sense, because o93 doesn't have puromycin
#   
# - I expect o128 to work better, but it's more expensive, so I've been using 
#   o129 for most of my experiments so far.  In particular I used o129 in expt 
#   #101, which is the exiperiment I'm most trying to follow here.
#
# - So I'll use o129 to start, and come back to test o128 once I have so 
#   coupling activity to optimize.
#
# Spin filter step:
#
# - MW of f11: 121.342 kDa
# - MW of o129: 19.8988 kDa
# - MW of f118: ≈140 kDa
#
# - From Amicon:
#   
#   - "Determine the molecular weight (MW) of your macromolecule of interest, 
#     and choose an ultrafiltration membrane with a molecular weight cut-off 
#     (MWCO) that is 2-3 times smaller than the MW."
#     https://tinyurl.com/4ffxu8zb
#
#   - That implies that I want a MWCO around 50-70 kDa.
sw cdna/make f11 o129 -m 50 -v 10
