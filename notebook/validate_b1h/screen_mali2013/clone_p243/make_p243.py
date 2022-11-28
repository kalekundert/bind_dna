#!/usr/bin/env python3

from stepwise import Protocol, Reactions, Extra, ul
from stepwise_mol_bio.ligate import Ligate

# The `ligate` protocol doesn't have the ability to do a negative control, but 
# such a control is very important for this reaction, so I wrote this script to 
# explicitly include it.

assembly = Ligate.Fragment('f200'), Ligate.Fragment('f201')
ligate = Ligate([assembly])
rxn = ligate.reaction.reaction

# Reaction volume:
# - In principle, I would want to make this large enough to use up most of my 
#   fragments.
#
# - In practice, I haven't even been plating all the cells I've transformed, 
#   and it's nice to have leftover DNA in case I need to try again.
#
# - I chose '5 µL' because that enough to use about half of the f200 that I 
#   have.

rxn.hold_ratios.volume = '5 µL'

combos = [
        {'f200': 'f200 (high)', 'f201': '+'},
        {'f200': 'f200 (high)', 'f201': '−'},
]
rxns = Reactions(rxn, combos)
rxns.step_kind = 'ligation'
rxns.instructions = ul("Mix all reagents on ice, to prevent the reaction from starting before all of the fragments have been added.")

p = Protocol()
p += rxns
p += ligate.protocol.steps[-1]
p.print()
