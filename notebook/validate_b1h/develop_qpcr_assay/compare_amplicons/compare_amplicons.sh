#!/usr/bin/env bash
set -euo pipefail

# Template:
# - I used p224/p228 as my "standard" plasmids in expt #165, so I can just keep 
#   doing that here.
#
# −RT control:
# - I should include this, because I've been having problems with background 
#   DNA.
#
sw b1h/qpcr \
  sz224 sz228 \
  -P sr1/sr2   -P sr3/sr5 \
  -P o366/o367 -P o368/o369 \
  -P o370/o371 -P o372/o373 \
  -N -d2 \
  -Q \
  --save-extra-cdna

