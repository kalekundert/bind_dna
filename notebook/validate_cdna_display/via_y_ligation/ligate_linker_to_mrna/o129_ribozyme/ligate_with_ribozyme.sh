#!/usr/bin/env zsh
set -euo pipefail

sw cond ligate_with_ribozyme_v2.xlsx |
sw zap |

# Control volume:
# - I only need 3.2 µL for the gel.
# - But I've been aiming to get 4.0 µL of each sample, to have some extra.
# - I also don't want to pipet volumes less than 0.5 µL.
sw step "Prepare the 750 nM linker-only control:~0.5 µL 10 µM o129~6.17 µL nuclease-free water" |
sw step "Prepare the 500 nM mRNA-only controls:~2 µL 10 µM f123,f124~2 µL nuclease-free water" |

# Excess linker:
# - I want an excess of linker, because I'm testing to see if the 
#   ribozyme-cleaved mRNA is more reactive, so I want to be sure that the 
#   linker itself isn't a limiting reagent.
# - I'm using a 1.5x excess based on [Reyes2021].
#
# mRNA concentration:
# - The default is 2 µM.
# - My stocks are only 1 µM though, so I had to reduce this.
# - I chose 0.5 µM, because it's a round number and still concentrated enough 
#   to load the full amount on the gel.
sw cdna/ligate f123,o129 f124,o129 -V 4 -c 0.5 -a 3 -l 1.5 |
sw step "Prepare the −ligase controls:~2.56 µL annealed mRNA/linker~1.44 µL nuclease-free water" |
sw gel urea/ligate/zif 7 -c 500

