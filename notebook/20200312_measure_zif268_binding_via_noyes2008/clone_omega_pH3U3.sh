#!/usr/bin/env bash
set -euo pipefail

# The spacing between the target and the promoter on the pH3U3 plasmid is 
# different for ω- and α-based B1H.  Since I want to do ω-based B1H, but 
# ordered the α-based pH3U3, I need to fix the spacing.  See [Noyes2008], 
# Figure S6.

# Update: looking at these primers, I think this could be tough to clone.  I'm 
# just going to order the plasmid from addgene and save myself the trouble.

klab_mutagenesis \
  H3U3_10BP \
  cgtatcacgaggccctttcgtcttcaaacgcgtgtacacccgggcggccgctgcgtgggcggcactccggaggcgcgccgaattctttacactttatgcttccggctcgtatgttgtgtcgaccgagcggataacaatttcacacaggaaacagctatgaca \
  cgtatcacgaggccctttcgtcttcaaacgcgtgtacacccgggcggccgctgcgtgggcgggacgaattctttacactttatgcttccggctcgtatgttgtgtcgaccgagcggataacaatttcacacaggaaacagctatgaca

  
