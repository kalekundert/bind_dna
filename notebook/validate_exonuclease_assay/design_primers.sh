#!/usr/bin/env bash
set -euo pipefail

# Install initial target sequences

klab_mutagenesis \
  DNASE_5 \
  TAATACGACTCACTATAGGGAGTTCCGGATCCCATATGCCTGGGATATCAACAACAACAACAACCGCGTGTAAAATCCGAGAACCCA \
  TAATACGACTCACTATAGGGAGTTCCGATATCTAAATTGCGTGGGCGATTATAGGATCCTTGTTGGTTGGTCGAAGAACAACAACAACAACCGCGTGTAAAATCCGAGAACCCA \

klab_mutagenesis \
  DNASE_3_CDNA \
  CATCTGCGTCAGAAAGATGGCGGCGGCAGCTAAGGATCCCATATGCCTGGGATATCAACAACAACAACAACCGCGTGTAAAATCCGAGAACCC \
  CATCTGCGTCAGAAAGATGGCGGCGGCAGCTAAAGTTCCGATATCTAAATTGCGTGGGCGATTATAGGATCCTTGTTGGTTGGTCGAAGAACAACAACAACAACCGCGTGTAAAATCCGAGAACCC \

klab_mutagenesis \
  DNASE_3_CIS \
  GGGTTCTCGGATTTTACACGCGGTTGTTGTTGTTGTTGATATCCCAGGCATATGGGATCCGGAACTGCTAGCAAAAGGCCAGCAAAAGG \
  GGGTTCTCGGATTTTACACGCGGTTGTTGTTGTTGTTCTTCGACCAACCAACAAGGATCCTATAATCGCCCACGCAATTTAGATATCGGAACTGCTAGCAAAAGGCCAGCAAAAGG \

# Replace target with non-target (same primers work for all plasmids)

klab_mutagenesis \
  DNASE_NONTARGET \
  TAATACGACTCACTATAGGGAGTTCCGATATCTAAATTGCGTGGGCGATTATAGGATCCTTGTTGGTTGGTCGAAGAACAACAACAACAACCGCGTGTAAAATCCGAGAACCCA \
  TAATACGACTCACTATAGGGAGTTCCGATATCTAAATTGCGAAAGCGATTATAGGATCCTTGTTGGTTGGTCGAAGAACAACAACAACAACCGCGTGTAAAATCCGAGAACCCA \
