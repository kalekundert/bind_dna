#!/usr/bin/env python3

from stepwise import Protocol, Reactions, ul
from stepwise_mol_bio.ligate import Ligate

# The `ligate` protocol doesn't have the ability to do a negative control, but 
# such a control is very important for this reaction, so I wrote this script to 
# explicitly include it.

assembly = Ligate.Fragment('f195'), Ligate.Fragment('f196')
ligate = Ligate([assembly])
rxn = ligate.reaction.reaction
rxn.hold_ratios.volume = '10 µL'

combos = [
        {'f196': '+'},
        {'f196': '-'},
]
rxns = Reactions(rxn, combos)
rxns.step_kind = 'ligation'
rxns.instructions = ul("Mix all reagents on ice, so the reaction doesn't start before all of the fragments have been added.")

p = Protocol()
p += rxns
p += ligate.protocol.steps[-1]
p.print()

