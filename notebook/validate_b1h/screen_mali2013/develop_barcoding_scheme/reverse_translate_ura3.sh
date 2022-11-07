#!/usr/bin/env bash
set -euo pipefail

klab_reverse_translate \
  ura3_prot.fa \
  -d ura3_dna.fa \
  -R BsaI,BsmBI
