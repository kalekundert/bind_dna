#!/usr/bin/env python3

"""
Usage:
    debug_digest.py <product>
"""

import docopt
import stepwise

from stepwise import pl, ul
from stepwise_mol_bio import RestrictionDigest

args = docopt.docopt(__doc__)
p = stepwise.Protocol()

# 10x dilution of both enzymes, to keep the volumes reasonable.
p += pl(
        "Dilute EcoRI-HF and NotI-HF to 2 U/µL in 1x CutSmart:",
        ul(
            "8 µL water",
            "1 µL 10x CutSmart",
            "1 µL 20 U/µL EcoRI-HF/NotI-HF",
        ),
)

app = RestrictionDigest.from_product(args['<product>'])  # 'f79'
app.dna_ug = 0.1

rxn = app.reaction
rxn.num_reactions = 4

rxn['EcoRI-HF'].master_mix = False
rxn['NotI-HF'].master_mix = False

rxn['EcoRI-HF'].hold_conc.stock_conc = 2, 'U/µL'
rxn['NotI-HF'].hold_conc.stock_conc = 2, 'U/µL'

p += app.protocol

p.print()
