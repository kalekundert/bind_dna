#!/usr/bin/env bash
set -euo pipefail

klab_mutagenesis \
  CURE_NOTI \
  gctattgctgaaggtcgtgcggccgcggactacaaggatgacgacgacaag \
  gctattgctgaaggtcgtgcggcAgcggactacaaggatgacgacgacaag \
  -T 64
