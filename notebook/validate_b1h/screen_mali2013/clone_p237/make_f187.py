#!/usr/bin/env python3

import stepwise
from stepwise import Reactions, AutoMix
from stepwise_mol_bio import RestrictionDigest

enzymes = ['NotI-HF', 'HindIII-HF']

digest = RestrictionDigest.from_tags(['l3'], enzymes)

# DNA quantity:
# - The goal is for this concentration to be reasonably well-matched to the 
#   concentration of f186, so that I don't need to pipet very small volumes 
#   when setting up the ligation reaction.
#
# - My current prep of f186 is very dilute, so 1 ng/µL is a reasonable target 
#   for f187.  If I had more f186, though, I'd want this reaction to be more 
#   concentrated.

digest.dna_ug = 0.02
digest.target_volume_uL = 20

rxn = digest.reaction.reaction

def init_mix(mix):
    mix.volume = rxn.volume / 2
    mix.name = 'enzyme'

mix = AutoMix({'water', 'buffer', *enzymes}, init=init_mix)

rxns = Reactions(rxn, mixes=[mix])
rxns.step_kind = 'restriction digestion'
rxns.extra.min_volume = '0.5 µL'

p = stepwise.Protocol()
p += rxns
p += digest.protocol.steps[-1]
p.print()


