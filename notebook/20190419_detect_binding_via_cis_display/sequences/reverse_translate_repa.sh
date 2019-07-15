#!/usr/bin/env bash

klab_reverse_translate repa_$1.prot \
  --restriction-sites EcoRV,NruI,BstNI,NdeI,BamHI,BbsI,BsaI,BsmBI,BtgZI,SapI,BbvCI,XmnI \
  --template-fasta repa_$1.fa \
  | tee repa_$1_no_digest.fa
  
