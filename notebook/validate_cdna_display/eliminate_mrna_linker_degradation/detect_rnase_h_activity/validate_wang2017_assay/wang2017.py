#!/usr/bin/env python3

import stepwise
import autoprop

@autoprop
class IncubateAndMeasure:

    def __init__(self, num_reactions):
        self.num_reactions = num_reactions

    def __iter__(self):
        yield from self.protocol.steps

    def get_reaction(self):
        rxn = stepwise.MasterMix()
        rxn.volume = '100 µL'
        rxn.num_reactions = self.num_reactions

        rxn['MgCl₂'].stock_conc = '2000 mM'
        rxn['MgCl₂'].hold_stock_conc.conc = '1 mM'
        rxn['MgCl₂'].master_mix = True

        rxn['o213'].stock_conc = '100 µM'
        rxn['o213'].hold_stock_conc.conc = '0.1 µM'
        rxn['o213'].master_mix = True

        rxn['samples'].volume = '10.5 µL'
        rxn['samples'].master_mix = False

        while rxn['MgCl₂'].volume * rxn.scale < '1 µL':
            rxn['MgCl₂'].hold_conc.stock_conc /= 2
        while rxn['o213'].volume * rxn.scale < '1 µL':
            rxn['o213'].hold_conc.stock_conc /= 2

        return rxn

    def get_protocol(self):
        p = stepwise.Protocol()

        p += stepwise.Step("""\
                Incubate at 30°C for 10 min.
        """)

        p += stepwise.Step("""\
                Setup the following reactions in a black,
                opaque-bottomed 96-well plate:
        """, self.reaction)

        p += stepwise.Step("""\
                Measure flourescence (ex: 450, em: 521) in a plate 
                reader for 30 min at 30°C without shaking.
        """)

        return p


