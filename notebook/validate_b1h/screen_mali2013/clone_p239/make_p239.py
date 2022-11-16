#!/usr/bin/env python3

from stepwise import Protocol, Reactions, Extra, ul
from stepwise_mol_bio.ligate import Ligate

# The `ligate` protocol doesn't have the ability to do a negative control, but 
# such a control is very important for this reaction, so I wrote this script to 
# explicitly include it.

assembly = Ligate.Fragment('f190'), Ligate.Fragment('f191')
ligate = Ligate([assembly])
rxn = ligate.reaction.reaction

# If I keep this, it messes up the "min volume" calculation.
del rxn['water']

# Reaction volume:
# - Ideally this would be large enough to use most of the ≈18 µL of f190 that I 
#   purified.
# - However, I only have 2.5 µL of f191 left, so I had to scale down the 
#   reaction to accommodate that.
rxn.hold_ratios.volume = '3 µL'


combos = [
        {'f191': '+'},
        {'f191': '-'},
]
rxns = Reactions(rxn, combos)
rxns.step_kind = 'ligation'
rxns.instructions = ul("Mix all reagents on ice, so the reaction doesn't start before all of the fragments have been added.")

# Extra:
# - Because I scaled down the reaction (see above), the master mix volumes were 
#   very small.
# - None of the master mix reagents are that precious, so I scaled up to more 
#   comfortable levels.
rxns.extra = Extra(min_volume='1 µL')

p = Protocol()
p += rxns
p += ligate.protocol.steps[-1]
p.print()
