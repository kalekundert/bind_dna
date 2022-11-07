#!/usr/bin/env python3

import stepwise
from stepwise import Reactions, AutoMix, iter_all_mixes
from stepwise_mol_bio import RestrictionDigest

enzymes = ['NotI-HF', 'HindIII-HF', 'XmnI']

digest = RestrictionDigest.from_tags(['p236'], enzymes)
digest.dna_ug = 0.1

rxn = digest.reaction.reaction
combos = [
        {'NotI-HF': '−', 'HindIII-HF': '−'},
        {'NotI-HF': '+', 'HindIII-HF': '−'},
        {'NotI-HF': '−', 'HindIII-HF': '+'},
        {'NotI-HF': '+', 'HindIII-HF': '+'},
]

def init_mix(mix):
    mix.volume = rxn.volume / 2
    mix.name = 'enzyme'

    for child in iter_all_mixes(mix):
        if 'XmnI' in child.reagents:
            child.name = 'XmnI'

mix = AutoMix({'water', 'buffer', *enzymes}, init=init_mix)

rxns = Reactions(rxn, combos, mixes=[mix])
rxns.step_kind = 'restriction digestion'
rxns.extra.min_volume = '0.5 µL'

p = stepwise.Protocol()
p += rxns
p += digest.protocol.steps[-1]
p += stepwise.load('gel tbe −/−,+/−,−/+,+/+ -v 10 -V 12.5')
p.print()



