#!/usr/bin/env bash
set -euo pipefail

sw zap |

sw cond optimize_nb_gfp_expression.xlsx |

# DNA concentrations:
# - NEB claims "typically, the optimal amount [of DNA template] will fall in a 
#   range of 25-1000 ng template [for a 25 µL reaction]." 
#   https://www.neb.com/protocols/0001/01/01/protein-synthesis-reaction-using-purexpress-e6800
# - For a 5 µL reaction, that's 5-200 ng.
# - I estimate that f144 will be 200 ng/µL.
# - My gels have 15 lanes:
#   - 1 ladder
#   - 1 −template control
#   - 7 RNA concentrations
#   - ≤6 DNA concentrations
#     - I don't care as much about finding the optimal DNA concentration, so 
#       I'm ok with having fewer steps in the DNA titration.
#     - 4 steps: 3 aliquots
#     - 5+ steps: 4 aliquots
#     - I should just do 6 steps.
#
# DNA volumes:
# - This DNA is less precious to me, so I use larger volume to make the 
#   pipeting easier.
sw serial 5 '200 ng/µL' to 5 6 -m f144 -d 'nuclease-free water' |

# RNA concentrations:
# - My stock is 10 µM.
# - In previous experiments I spanned 2 logs (see expt #99) in 7 steps.
# - So stock concentrations of 10 µM - 100 nM seem reasonable.

# RNA volumes:
# - Didn't really think about this.  Don't want to use too much, but also want 
#   pipeting to be accurate.
sw serial 3 '10000 nM' to 100 7 -0 -m f145 -d 'nuclease-free water' | 

# Hack:
# - This is a mixed RNA/DNA reaction.
# - `sw ivtt` doesn't like that, because it can't calculate concentration 
#   correctly.
# - I don't care, so I specify `-r` to stop `sw ivtt` from complaining.
# - Really, I need to make the stock conc system more flexible.  In this case, 
#   I'd really just like to ignore the stock conc entirely.
sw ivtt f144 f145 -p purex/lys -n 14 -v 5 -V 1.15 -r |
sw fluorotect_rnase_digest -V 5 |

sw gel sds/ivtt/nb 14
