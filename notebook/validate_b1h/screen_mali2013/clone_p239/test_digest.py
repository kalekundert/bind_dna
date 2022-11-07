#!/usr/bin/env python3

import stepwise
from stepwise import Reactions, AutoMix
from stepwise_mol_bio import RestrictionDigest

# Setup the reaction as a triple-digest, so that we get a buffer and incubation 
# parameters that are compatible with all the enzymes.  Later we'll modify the 
# reaction such that we only use one enzyme at a time.
enzymes = ['SapI', 'HaeII', 'NdeI']
digest = RestrictionDigest.from_tags(['p239'], enzymes)

# Just set this to a small number, so that we get a 10 µL reaction.
digest.dna_ug = 0.01

rxn = digest.reaction.reaction

for key in enzymes:
    del rxn[key]

# Whether or not these volume make sense will depend on the concentration of 
# the plasmid, but in most cases this will be a huge excess of enzyme.  For a 
# test digest, pipetting the "right" amount of material would be more effort 
# than it's worth.
rxn['DNA'].volume = '1.0 µL'
rxn['enzyme'].volume = '0.5 µL'

combos = [
        {'enzyme': 'SapI'},
        {'enzyme': 'HaeII'},
        {'enzyme': 'NdeI'},
]

rxns = Reactions(rxn, combos)
rxns.step_kind = 'restriction digestion'

p = stepwise.Protocol()
p += rxns
p += digest.protocol.steps[-1]
p += stepwise.load('e_gel')
p.print()


