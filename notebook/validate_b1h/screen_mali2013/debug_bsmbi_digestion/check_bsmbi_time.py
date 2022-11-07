#!/usr/bin/env python3

import stepwise
from stepwise import Reactions, AutoMix, replace_text
from stepwise_mol_bio import RestrictionDigest

enzymes = ['BsmBI-v2']

digest = RestrictionDigest.from_tags(['f206'], enzymes)
digest.dna_ug = 0.015

rxn = digest.reaction.reaction

def init_mix(mix):
    mix.volume = rxn.volume / 2
    mix.name = 'BsmBI'

mix = AutoMix({'water', 'buffer', *enzymes}, init=init_mix)

rxns = Reactions(rxn, replicates=3, mixes=[mix])
rxns.step_kind = 'restriction digestion'
rxns.extra.min_volume = '0.5 µL'

p = stepwise.Protocol()
p += rxns
p += replace_text(digest.protocol.steps[-1], '5.15', '15, 60, and 120')
p += stepwise.load('gel tbe 15m,1h,2h -v 10 -V 12.5 -L "5 µL TriDye ultra low range ladder"')
p.print()



