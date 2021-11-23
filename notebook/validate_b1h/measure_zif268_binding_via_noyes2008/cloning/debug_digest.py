#!/usr/bin/env python3

from stepwise_mol_bio import RestrictionDigest

app = RestrictionDigest.from_product('f79')
app.dna_ug = 0.1

rxn = app.reaction
rxn.num_reactions = 4

rxn['EcoRI-HF'].master_mix = False
rxn['NotI-HF'].master_mix = False

rxn['EcoRI-HF'].hold_conc.stock_conc = 2, 'U/µL'
rxn['NotI-HF'].hold_conc.stock_conc = 2, 'U/µL'

app.protocol.print()
