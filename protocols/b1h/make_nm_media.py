#!/usr/bin/env python3

import stepwise
import byoc
import autoprop

from stepwise import pl
from stepwise_mol_bio import Main
from byoc import DocoptConfig
from operator import not_


@autoprop
class MakeNmMedia:
    """\
Prepare permissive and selective liquid media for auxotrophy-based B1H assays.

Usage:
    make_nm_media.py <volume_mL> [-P] [-S] [-s <concs>] [-a <antibiotics>]

Options:
    -P --no-permissive
        Don't make permissive media (e.g. +His, +Ura).

    -S --no-selective
        Don't make selective media (e.g. −His, −Ura, +3AT).

    -s --3at-concs <mM,...>             [default: 10]
        What concentrations of 3-AT to include in the selective media.  
        Specify multiple concentrations separated by commas.

    -a --antibiotics <names,...>        [default: carb]
        What antibiotics to include in the media.  Must be either 'carb', 
        'kan', or 'carb,kan'.
"""
    __config__ = [DocoptConfig]

    def get_protocol(self):
        p = stepwise.Protocol()

        if self.include_permissive:
            p += pl(
                    f"Prepare {self.volume_mL:g} mL 8/7x permissive media:",
                    permissive_recipe(self.volume_mL, self.antibiotics),
            )

        if self.include_selective:
            for conc_mM in self.selective_concs_mM:
                p += pl(
                        f"Prepare {self.volume_mL:g} mL 8/7x selective media with {conc_mM:g} mM 3-AT:",
                        selective_recipe(self.volume_mL, conc_mM, self.antibiotics),
                )

        return p

    volume_mL = byoc.param(
            '<volume_mL>',
            cast=float,
    )
    selective_concs_mM = byoc.param(
            '--3at-concs',
            cast=lambda x: map(float, x.split(',')),
    )
    antibiotics = byoc.param(
            '--antibiotics',
            cast=lambda x: set(x.split(',')),
    )
    include_permissive = byoc.param(
            '--no-permissive',
            cast=not_,
            default=True,
    )
    include_selective = byoc.param(
            '--no-selective',
            cast=not_,
            default=True,
    )

def permissive_recipe(volume_mL, antibiotics):
    # Histidine concentration:
    # - [Noyes2008] calls for 0.1% His.
    # - I assume this is 0.1% (w/v) = 1 mg/mL
    # - MW=155.15 g/mol, 1 mg/mL = 6.44 mM
    # - My stock histidine is 100 mM, so it's more convenient to use a 
    #   molar concentration for this (e.g. 6 mM).
    #
    # Uracil concentration:
    # - [Noyes2008] consistently uses 200 µM uracil for permissive 
    #   conditions.

    nm = stepwise.MasterMix("""\
            Reagent   Stock    Conc       Volume
            =======  ======  ======  ===========
            NM                       to 70000 µL
            His      100 mM    6 mM
            Ura       20 mM  0.2 mM
    """)

    nm['His'].volume *= 8/7
    nm['Ura'].volume *= 8/7
    add_antibiotics(nm, antibiotics)
    nm.hold_ratios.volume = volume_mL * 1000, 'µL'
    nm.show_concs = True
    return nm

def selective_recipe(volume_mL, conc_mM, antibiotics):
    # 3-AT concentrations:
    # - From previous experience, the AAA-target strain (s5) cannot survive 
    #   with even 0 mM 3-AT.  So there isn't much point titrating low 
    #   concentrations of 3-AT (originally I did 1, 2, 4 mM).
    # - [Noyes2008] uses 5 and 10 mM 3-AT for selections.  Figure S5 
    #   furthermore shows that the TGG-target strain can survive 50 mM 3-AT.
    # - My goal is just to show that the strains and plasmids perform as 
    #   expected, and eventually to evaluate my single-plasmid constructs.  For 
    #   the latter, it might be useful to do a titration to see how much 
    #   more/less stringent the single plasmid is.  But I think it's reasonable 
    #   to use 10 mM 3-AT to just ask, "does these plasmids/strains mostly 
    #   work".

    nm = stepwise.MasterMix()
    nm.volume = volume_mL * 1000, 'µL'
    nm.solvent = 'NM'
    nm['3-AT'].stock_conc = '1000 mM'
    nm['3-AT'].hold_stock_conc.conc = 8/7 * conc_mM, 'mM'
    add_antibiotics(nm, antibiotics)
    nm.show_concs = True
    return nm

def add_antibiotics(rxn, antibiotics):
    if 'kan' in antibiotics:
        rxn['kanamycin'].stock_conc = '1000x'
        rxn['kanamycin'].hold_stock_conc.conc = 8/7 * 0.5, 'x'

    if 'carb' in antibiotics:
        rxn['carbenicillin'].stock_conc = '1000x'
        rxn['carbenicillin'].hold_stock_conc.conc = 8/7 * 1, 'x'


if __name__ == '__main__':
    app = MakeNmMedia()
    byoc.load(app)
    app.protocol.print()


