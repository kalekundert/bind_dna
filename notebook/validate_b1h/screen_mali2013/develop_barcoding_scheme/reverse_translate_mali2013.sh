#!/usr/bin/env bash
set -euo pipefail

klab_reverse_translate \
  mali2013.fa \
  -d mali2013_szf1_dna.fa \
  -R BsaI,BsmBI,BbsI
