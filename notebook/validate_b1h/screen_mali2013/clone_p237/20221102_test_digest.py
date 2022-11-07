#!/usr/bin/env python3

import stepwise
from stepwise import Reactions, AutoMix
from stepwise_mol_bio import RestrictionDigest

# Setup the reaction as if we'll use all the enzymes at once, so that we get 
# buffer and incubation parameters that are compatible with everything 
# (although I manually confirmed that all of these enzymes are recommended for 
# use with CutSmart buffer).  Later we'll modify the reaction such that we only 
# use one enzyme at a time.
plasmids = ['p237', 'p239']
enzymes = ['Esp3I', 'HaeII', 'BsrGI-HF', 'EcoRV-HF']

digest = RestrictionDigest.from_tags(plasmids, enzymes)
digest.time = '2 hr'

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
rxn['DNA'].name = 'DNA'
rxn['enzyme'].volume = '0.5 µL'

combos = [
        {'enzyme': enzyme, 'DNA': plasmid}
        for enzyme in ['water', *enzymes]
        for plasmid in plasmids
]

def init_mix(mix):
    mix.volume = rxn.volume / 2
    mix.name = 'enzyme'

mix = AutoMix({'water', 'buffer', 'enzyme'}, init=init_mix)

rxns = Reactions(rxn, combos, mixes=[mix])
rxns.step_kind = 'restriction digestion'

p = stepwise.Protocol()
p += rxns
p += digest.protocol.steps[-1]
p += stepwise.load('e_gel')
p.print()


