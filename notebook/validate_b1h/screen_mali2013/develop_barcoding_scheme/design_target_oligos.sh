#!/usr/bin/env bash

dbp_expand_features \
  ../../../../sequences/libraries/l3.dna \
  library.xlsx \
  -f 'barcode:Barcode Seq' -f 'target:Sequence' \
  -s Targets

