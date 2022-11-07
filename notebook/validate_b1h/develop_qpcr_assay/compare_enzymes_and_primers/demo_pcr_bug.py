#!/usr/bin/env python3

from stepwise import Reaction, Reactions
from stepwise.reaction.mix import *

# rxn = Reaction.from_text("""\
# Reagent  Volume
# =======  ======
# a          1 µL
# b          2 µL
# c          3 µL
# d          4 µL
# e          5 µL
# """)

# combos = []
# for i in range(3):
#     combos += [
#             {'c': '1', 'd': '1', 'e': f'{i}'},
#             {'c': '2', 'd': '2', 'e': f'{i}'},
#     ]

# rxns = Reactions(rxn, combos)

# mix = Mix({
#     Mix({
#         Mix({'a', 'b'}),
#         'c', 'd',
#     }),
#     'e',
# })

# order_map = make_order_map(rxn.keys())
# score = score_mix(mix, combos, 0, [], order_map, {})
# debug(mix, score, order_map)

# rxns.protocol.print()

#plan_mixes(rxn, combos)


# The problem is that {{'a', 'b', 'c'}, 'd'} doesn't get generated.  And that's 
# because of the pairwise nature of the algorithm.
#
# What if levels should be a dict.  And the goal is to remove the top level of 
# the dict: i.e. pop the top level and merge those components into lower levels
#
# matching reagents must go together, but can go as mix or reagents.

x = iter_merged_top_levels(
        {'a', 'b', 'c'},
        {'d'},
        {'b': {'b', 'c'}, 'c': {'b', 'c'}},
        {},
)
debug(list(x))
