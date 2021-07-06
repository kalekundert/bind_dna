#!/usr/bin/env python3

import stepwise
import autoprop
import appcli

from stepwise import pl, ul
from stepwise_mol_bio import (
        Main, Argument, ShareConfigs,
        bind_arguments, int_or_expr, float_or_expr, merge_names,
)
from freezerbox import ReagentConfig
from appcli import Key, Method, DocoptConfig
from more_itertools import all_equal
from inform import plural

class MrnaConfig(ReagentConfig):
    tag_getter = lambda obj: obj.mrna

class LinkerConfig(ReagentConfig):
    tag_getter = lambda obj: obj.linker

def parse_mrna_linker_pairs(strs):
    return [parse_mrna_linker_pair(x) for x in strs]

def parse_mrna_linker_pair(str):
    mrna, linker = str.split(',')
    return LigateMrnaLinker.Reaction(mrna, linker)

@autoprop.cache
class LigateMrnaLinker(Main):
    """
Ligate a puromycin linker to an mRNA transcript in preparation for 
mRNA-display.

This protocol is based on [Naimudden2016], [Reyes2021], and my own 
experience.

Usage:
    ligate <mrna,linker>... [-n <int>] [-v <µL>] [-C <µM>] [-l <ratio>]
        [-L <µM>]

Options:
    -n --num-reactions <int>
        The number of reactions to setup up.  The default is the number of 
        mRNA/linker pairs specified on the command line.

    -v --mrna-volume <µL>
        The volume of mRNA to use in each annealing reaction, in µL.  This will 
        scale to total volume of the reaction proportionally.

    -C --mrna-stock <µM>
        The stock concentration of the mRNA, in µM.  The volume of mRNA will be 
        updated accordingly to keep the amount of material in the reaction 
        constant.  The default is read from the FreezerBox database, or 10 µM 
        if the given mRNA is not in the database.  Use '--mrna-volume' to 
        change the amount of mRNA in the reaction.

    -l --linker-ratio <float>           [default: ${app.linker_ratio}]
        The amount of linker to add to the reaction, relative to the amount 
        of mRNA in the reaction.

    -L --linker-stock <µM>
        The stock concentration of the linker, in µM.  The volume of linker be 
        updated accordingly, to keep the amount of material in the reaction 
        constant.  The default is read from the FreezerBox database, or 10 µM 
        if the given linker is not in the database.
    """
    __config__ = [
            DocoptConfig,
    ]

    class Reaction(ShareConfigs, Argument):
        __config__ = [
                MrnaConfig,
                LinkerConfig,
        ]

        mrna_stock_uM = appcli.param(
                Key(DocoptConfig, '--mrna-stock'),
                Key(MrnaConfig, 'stock_conc_uM'),
                default=10,
        )
        linker_stock_uM = appcli.param(
                Key(DocoptConfig, '--linker-stock'),
                Key(LinkerConfig, 'stock_conc_uM'),
                default=10,
        )

        def __init__(self, mrna, linker, **kwargs):
            self.mrna = mrna
            self.linker = linker
            self._set_known_attrs(kwargs)

    reactions = appcli.param(
            Key(DocoptConfig, '<mrna,linker>', cast=parse_mrna_linker_pairs),
            get=bind_arguments,
    )
    num_reactions = appcli.param(
            Key(DocoptConfig, '--num-reactions', cast=int_or_expr),
            Method(lambda self: len(self.reactions)),
    )
    mrna_volume_uL = appcli.param(
            Key(DocoptConfig, '--mrna-volume', cast=float_or_expr),
            default=None,
    )
    linker_ratio = appcli.param(
            Key(DocoptConfig, '--linker-ratio', cast=float),
            default=1,
    )

    def __init__(self, reactions):
        self.reactions = reactions

    def get_protocol(self):
        p = stepwise.Protocol()
        n = self.num_reactions

        f = "Using 0.6x linker reduces the amount of unligated linker, see expt #1"
        p += pl(
                f"Setup {plural(n):# annealing reaction/s}{p.add_footnotes(f)}:",
                self.annealing_reaction,
        )

        # This thermocycler protocol is from [Reyes2021].  From my own 
        # experimentation, I expect that this step may not even be necessary, 
        # but I haven't thoroughly tested that yet.
        p += pl(
                "Incubate as follows:",
                ul(
                    "90°C for 30s",
                    "Cool to 25°C at 1°C/s",
                ),
        )

        p += pl(
                f"Setup {plural(n):# ligation reaction/s}:",
                self.ligation_reaction,
        )
        p += pl(
                "Incubate as follows:",
                ul(
                    "25°C for 10 min",
                    "65°C for 10 min",
                ),
        )
        return p

    def get_combined_reaction(self):
        rxn = stepwise.MasterMix("""\
            Reagent                 Stock      Volume  MM?
            ====================  =======  ==========  ===
            nuclease-free water                3.8 µL   +
            T4 RNA ligase buffer      10x      1.0 µL   +
            ATP                     10 mM      1.0 µL   +
            mRNA                    10 µM      2.0 µL
            linker                  10 µM      2.0 µL
            T4 RNA ligase         10 U/µL      0.2 µL   +
        """)
        rxn.num_reactions = n = self.num_reactions

        rxn['mRNA'].hold_conc.stock_conc = \
                min(x.mrna_stock_uM for x in self.reactions), 'µM'
        rxn['linker'].hold_conc.stock_conc = \
                min(x.mrna_stock_uM for x in self.reactions), 'µM'

        rxn['linker'].volume *= self.linker_ratio

        if uL := self.mrna_volume_uL:
            rxn.hold_ratios.volume *= uL / rxn['mRNA'].volume.value

        rxn['mRNA'].master_mix = all_equal(x.mrna for x in self.reactions)
        rxn['linker'].master_mix = all_equal(x.linker for x in self.reactions)

        rxn['mRNA'].name = merge_names(x.mrna for x in self.reactions)
        rxn['linker'].name = merge_names(x.linker for x in self.reactions)

        rxn.fix_volumes('mRNA', 'nuclease-free water')
        return rxn

    def get_annealing_reaction(self):
        rxn = self.combined_reaction.copy()

        del rxn['nuclease-free water']
        del rxn['ATP']
        del rxn['T4 RNA ligase']

        return rxn

    def get_ligation_reaction(self):
        rxn = self.combined_reaction.copy()
        anneal_rxn = self.annealing_reaction

        for reagent in anneal_rxn:
            del rxn[reagent.key]

        rxn['annealed mRNA/linker'].volume = anneal_rxn.volume
        rxn['annealed mRNA/linker'].master_mix = \
                anneal_rxn['mRNA'].master_mix and \
                anneal_rxn['linker'].master_mix

        return rxn

if __name__ == '__main__':
    LigateMrnaLinker.main()
