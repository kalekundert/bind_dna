#!/usr/bin/env bash
set -euo pipefail

# Confirm volumes...

sw zap |

sw cond optimize_purefrex1_zif_cond.xlsx |

#~ Volume:
#~ - f11 aliquots are 5 µL, so I tried to increase the volume as much 
#~   as possible (to make pipetting more accurate) without using more than one 
#~   aliquot.
sw serial 3.0 '10000 nM' to 100 6 -0 -m f11 -d 'nuclease-free water' |

# mRNA concentration:
# - f11 aliquots are 10 µM.
# - That corresponds to a maximum template concentration of 3.3 µM.
# - I want to use the highest concentration as possible, because I like being 
#   able to see a clear decrease in expression with too much template.
# - I'm going to make 33 nM the lowest concentration (100 nM template stock):
#   - This is the same dilution as before, just more volume per reaction.
#   - 20-40 nM has had pretty low expression in my other reactions.
#   - It doesn't really matter.
#
# DNA concentration:
# - GeneFrontier calls for DNA templates to have a 2 nM final concentration:
#   https://purefrex.genefrontier.com/resources/tech_notes/template-dna-preparation.html/
#   It's not totally clear which version they're talking about, but I got this 
#   link from the 1.0 product page.
# - 6 nM stock required to reach 2 nM final (based on the RNA stock 
#   concentrations in the following `ivtt` command)
#
# DNA template:
# - Would have preferred to use f2, but don't have any left.  So using p49 
#   instead.
# - p49 is miniprepped and not otherwise purified, so it is RNase-contaminated.  
#   GeneFrontier recommends purifying plasmid templates by phenol/chloroform 
#   extraction, but I'm just going to add RNase inhibitor to the reaction (also 
#   a recommended approach) and hope for the best.
sw dilute p49 -c 6nM -v2

sw ivtt p49 f11 -p frex1/zif -v 2.5 -n 8 -c 3300 -C 10000  |

sw step "Add 0.5 µL RNase cocktail (Invitrogen AM2286) to each reaction." |
sw step "Incubate at 37°C for 15 min." |

sw gel tricine/ivtt/zif 8

