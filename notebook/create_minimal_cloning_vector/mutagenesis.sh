#!/usr/bin/env bash
set -euo pipefail

# IDT won't make primers shorter than 15 bp, hence the minimum overlaps.

klab_mutagenesis \
  PUC_BSAI \
  GGCCCCAGTGCTGCAATGATACCGCGAGACCCACGCTCACCGGCTCCAGATTTATCAGCA \
  GGCCCCAGTGCTGCAATGATACCGCGAGAGCCACGCTCACCGGCTCCAGATTTATCAGCA \
  --min-overlap 14 --tm 64 

klab_mutagenesis \
  PUC_XMNI \
  ATAGCAGAACTTTAAAAGTGCTCATCATTGGAAAACGTTCTTCGGGGCGAAAACTCTCAAGGATCTTACC \
  ATAGCAGAACTTTAAAAGTGCTCATCATTGGAAAACGCTCTTCGGGGCGAAAACTCTCAAGGATCTTACC \
  --min-overlap 14