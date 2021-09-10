#!/usr/bin/env bash
set -euo pipefail

# DNA template:
# - I should be able to directly express p107 in PURExpress.
#
#   - The PCV2-Zif268 gene is expressed with a T7 promoter, and the lac 
#     operator should not interfere with expression.
#   - AmpR might also be expressed, which isn't ideal, but also isn't really a 
#     problem.
#
# - If I want to PCR amplify the gene, I'd need to order suitable a forward 
#   primer.

sw cond huh_purexpress_cond_pcv.xlsx

sw ivtt p107 -p purexpress |

# PCV2-Zif268 concentration
# - The ribosome concentration in PURExpress is 2 µM, and NEB estimates that 
#   each ribosome is recycled 5 times (or the DHFR control).  That implies a 
#   final DHFR concentration of 10 µM.
#
# - I'll have to measure protein concentration empirically, but I can use 10 µM 
#   as a starting point.
sw huh |

sw gel bolt/mes |

sw stain sypro-orange -I |
sw laser blue red
