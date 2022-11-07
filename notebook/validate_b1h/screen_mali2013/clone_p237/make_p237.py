#!/usr/bin/env python3

from stepwise import Protocol, Reactions
from stepwise_mol_bio.ligate import Ligate

# The `ligate` protocol doesn't have the ability to do a negative control, but 
# such a control is very important for this reaction, so I wrote this script to 
# explicitly include it.

assembly = Ligate.Fragment('f186'), Ligate.Fragment('f187')
ligate = Ligate([assembly])
rxn = ligate.reaction.reaction

combos = [
        {'f187': '+'},
        {'f187': '-'},
]
rxns = Reactions(rxn, combos)
rxns.step_kind = 'ligation'

p = Protocol()
p += rxns
p += ligate.protocol.steps[-1]
p.print()
