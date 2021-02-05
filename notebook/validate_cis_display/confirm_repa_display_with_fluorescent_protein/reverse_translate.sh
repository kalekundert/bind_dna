#!/usr/bin/env bash
set -euo pipefail

klab_reverse_translate \
  fp_candidates_protein.fa \
  -o fp_candidates_dna.tsv \
  -R BamHI,BbsI,BbvCI,BsaI,BstNI,BtgZI,EcoRV,Esp3I,NdeI,NruI,SapI,XmnI \
  -5 gacgcggtctcatATGGCTAGCTGGAGCCACCCGCAGTTCGAAAAAGGCGCC \
  -3 GGGGagagacctacaag \
  --no-stop-codon
