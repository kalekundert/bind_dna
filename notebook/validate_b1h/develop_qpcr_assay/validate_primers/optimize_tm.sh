#!/usr/bin/env bash
set -euo pipefail

# - PCR in regular thermocycler (I don't trust the qPCR gradient)
# - E-gel

sw pcr \
  -p ssoadv \
  p221,sr1,sr2 \
  p221,sr3,sr5 \
  p221,o366,o367 \
  p221,o368,o369 \
  p221,o370,o371 \
  p221,o372,o373 \
  -v 5 -d 10 -g 15 |
sw e_gel
