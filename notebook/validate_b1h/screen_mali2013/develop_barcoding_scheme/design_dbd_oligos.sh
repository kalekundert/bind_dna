#!/usr/bin/env bash

dbp_expand_features \
  ../../../../sequences/libraries/l4.dna \
  library.xlsx \
  -f 'barcode:Barcode Seq' -f 'sZF:Variable CDS Seq' \
  -s 'DBDs'

