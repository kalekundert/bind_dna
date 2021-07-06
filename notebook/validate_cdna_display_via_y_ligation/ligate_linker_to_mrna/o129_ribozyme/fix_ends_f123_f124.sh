#!/usr/bin/env bash
set -euo pipefail

# Strategy:
# - Gibson to make plasmid
# - PCR to get amplicon with desired end
# - Lyophilize to concentrate DNA
# - IVT
# - DNase treatment
# - Spin cleanup
# - PNK
# - gel purify
#   - product is likely to be messy
#   - even if 100% complete, better to be able to directly measure 
#     concentration.
#
# mRNA:
# - Zif
#   - Right size to easily purify.
#   - Easier cloning than FLAG, because Gibson should just work.
#   - Most relevant, ultimately.

# DNase treatment:
# - I haven't done this before, but it occurs to me that I don't want any 
#   template DNA in the protein expression reaction.  Without a DNase treatment 
#   step, I'm probably adding a lot of it.

# RNA yield
# - Not sure how to estimate this.  

# Buffer composition:
# - I'm just going to use the NEB T4 PNK buffer.  It's similar to the buffers 
#   used by all the reference protocols, but it's easier and more standardized.

# Enzyme concentration:
# - Reference protocols:
#   14:  6U /  50 pmol = 0.12 U/pmol
#   15:  1U / 100 pmol = 0.01 U/pmol
#   16: 10U / 300 pmol = 0.03 U/pmol
#
# - NEB T4 PNK: 10 U/µL
#
# - I have 100 µL × 4 µM ≈ 400 pmol RNA.
#   - 0.12 U/pmol = 48 U = 4.8 µL
#   - 0.05 U/pmol = 20 U = 2.0 µL
#   - 0.03 U/pmol = 12 U = 1.2 µL
#
# - I'll go with 0.05 U/pmol.
#   - Middle-of-the-road
#   - Not an unreasonable amount of enzyme.

sw reaction \
  -s "phosphatase" \
  "water;; to 110 µL; + " \
  "buffer; 10x; 11 µL; +" \
  "T4 PNK; 10 U/µL; 2 µL; +" \
  "f123,f124; 4 µM; 90 µL; -" \
  -n 2 \
  |

sw step "Incubate at 37°C for 6h" |

# Target concentration/volume:
# - Want to load 20 µg, 10 µL per lane
# - Loading buffer is 2x, so sample needs to be 20 µg/5 µL = 4 µg/µL
# - f123: 90 µL × 0.847 µg/µL ÷ 4 µg/µL = 19.1 µL
# - f124: 90 µL × 0.737 µg/µL ÷ 4 µg/µL = 16.6 µL
sw step "Concentrate samples to 4 µg/µL by lyophilization:~f123: 19.1 µL~f124: 16.6 µL" |

sw page_purify f123 f124 -b '404 nt' -c 4 -v 20 -C |

sw spin_cleanup zymo/rna-clean-conc/25 -s 400

