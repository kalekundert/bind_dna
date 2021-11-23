#!/usr/bin/env python3

import stepwise
from stepwise_mol_bio import Ivtt

p = stepwise.Protocol()

# Template concentration:
# - NEB calls for 2 µL of 125 ng/µL supplied control template in 25 µL 
#   reaction.
# - The default reaction parameters in swmb are based on the DHFR template 
#   (converted into molar units), so I shouldn't have to change anything.

# Visualization:
# - Coomassie:
#   - This is what NEB does.
#   - The DHFR band does not overlap with any PURExpress components
#
# - FluoroTect GreenLys
#   - DHFR is 20 kDa
#   - I can use RNase.
#   - This should be more sensitive to low expression levels.
#   - PUREfrex2 manual shows gel with DHFR and FluoroTect, and the 
#     visualization is very good.
#
# - I think FluoroTect is the way to go.

# Reaction volume:
# - I wanted to do 5 µL, but the volumes were too small to pipet.
# - I created a master mix for all the reagents except PURExpress itself, and 
#   set the minimum volume to 0.5 µL.  That's enough still more than enough for 
#   10 µL reactions, but I think 10 µL should be comfortable.

ivtt = Ivtt.from_tags(['DHFR'])
ivtt.preset = 'purex/lys'
ivtt.num_reactions = 2
ivtt.volume_uL = 10
ivtt.reaction['water'].master_mix = True
ivtt.reaction['solution A'].master_mix = False
ivtt.reaction['solution B'].master_mix = False
ivtt.reaction['RNase inhibitor, murine'].master_mix = True
ivtt.reaction['FluoroTect GreenLys'].master_mix = True
ivtt.reaction['DNA'].master_mix = True
ivtt.reaction.extra_min_volume = 0.5, 'µL'
ivtt.setup_instructions = [
        "Setup the above reaction with both new and old PURExpress.",
        *ivtt.setup_instructions,
]

p += ivtt.protocol

p += stepwise.load('fluorotect_rnase_digest -V 10')

p += stepwise.load('gel bolt/ivtt/dna 2')

p.print()

