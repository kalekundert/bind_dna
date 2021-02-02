#!/usr/bin/env python3

import stepwise
import autoprop

@autoprop
class IncubateAndMeasure:

    def __init__(self, num_reactions, rnase_incubation_time_min=10, rnase_incubation_temp_C=30, fam_control=False):
        self.num_reactions = num_reactions
        self.rnase_incubation_time_min = rnase_incubation_time_min
        self.rnase_incubation_temp_C = rnase_incubation_temp_C
        self.fam_control = fam_control

    def __iter__(self):
        yield from self.protocol.steps


    def get_dnazyme_buffer(self):
        rxn = stepwise.MasterMix()
        rxn.volume = 20 * (self.num_reactions + 5), 'µL'
        rxn.extra_min_volume = '1 µL'

        rxn['Tris pH=8.5'].stock_conc = '1000 mM'
        rxn['Tris pH=8.5'].hold_stock_conc.conc = '250 mM'

        rxn['NaCl'].stock_conc = '5000 mM'
        rxn['NaCl'].hold_stock_conc.conc = '500 mM'

        rxn['MgCl₂'].stock_conc = '2000 mM'
        rxn['MgCl₂'].hold_stock_conc.conc = '5 mM'

        if rxn['MgCl₂'].volume < '1 µL':
            rxn.hold_ratios.volume *= '1 µL' / rxn['MgCl₂'].volume 

        return rxn

    def get_reaction(self):
        rxn = stepwise.MasterMix()
        rxn.volume = '100 µL'
        rxn.num_reactions = self.num_reactions

        rxn['water'].name = 'nuclease-free water'

        rxn['DNAzyme buffer'].stock_conc = '5x'
        rxn['DNAzyme buffer'].hold_stock_conc.conc = '1x'
        rxn['DNAzyme buffer'].master_mix = True

        rxn['o213'].stock_conc = '100 µM'
        rxn['o213'].hold_stock_conc.conc = '0.1 µM'
        rxn['o213'].master_mix = True

        rxn['samples'].volume = '10.5 µL'
        rxn['samples'].master_mix = False

        while rxn['o213'].volume * rxn.scale < '0.5 µL':
            rxn['o213'].hold_conc.stock_conc /= 2

        return rxn

    def get_protocol(self):
        p = stepwise.Protocol()

        p += stepwise.Step(f"""\
                Incubate at {self.rnase_incubation_temp_C}°C for {self.rnase_incubation_time_min} min.
        """)

        p += stepwise.Step("Prepare 5x DNAzyme buffer:", self.dnazyme_buffer)

        p += stepwise.Step("""\
                Setup the following reactions in a black,
                opaque-bottomed 96-well plate:
        """, self.reaction)

        if self.fam_control:
            p += """\
Prepare a FAM control:

10 µM o42:
- 9 µL nuclease-free water
- 1 µL 100 µM o42

FAM control:
- 68.50 µL nuclease-free water
- 20.00 µL 5x DNAzyme buffer
-  1.00 µL 10 µM o42
- 10.50 µL blank
            """

        p += stepwise.Step("""\
                Measure flourescence (ex: 450, em: 521) in a plate 
                reader for 30 min at 30°C without shaking.
        """)

        return p


