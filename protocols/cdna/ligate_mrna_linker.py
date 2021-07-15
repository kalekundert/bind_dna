#!/usr/bin/env python3

import stepwise
import autoprop
import appcli

from stepwise import Quantity, pl, ul
from stepwise_mol_bio import (
        Main, Argument, ShareConfigs,
        bind_arguments, int_or_expr, float_or_expr, merge_names,
)
from freezerbox import ReagentConfig
from appcli import Key, Method, DocoptConfig
from more_itertools import all_equal
from inform import warn, plural

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
    """\
Ligate a puromycin linker to an mRNA transcript in preparation for 
mRNA-display.

This protocol is based on [Naimudden2016], [Reyes2021], and my own 
experience.

Usage:
    ligate <mrna,linker>... [-n <int>] [-v <µL> | -V <µL>] [-c <µM>] [-C <µM>]
        [-l <ratio>] [-L <µM>] [-a <µL>]

Options:
    -n --num-reactions <int>
        The number of reactions to setup up.  The default is the number of 
        mRNA/linker pairs specified on the command line.

    -v --mrna-volume <µL>
        The volume of mRNA to use in each annealing reaction, in µL.  This will 
        scale to total volume of the reaction proportionally.

    -V --total-volume <µL>
        The total volume of the reaction.  The volumes of all other reagents 
        will be scaled proportionally.

    -c --mrna-conc <µM>
        The final concentration of the mRNA, in µM.

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

    -a --extra-anneal-rxn <µL>
        How much extra annealing reaction to prepare, in µL.  By default, no 
        extra is prepared (it is assumed that the ligase master mix will be 
        directly added to the entire annealing reaction).
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
                Key(MrnaConfig, 'conc_uM'),
                default=10,
        )
        linker_stock_uM = appcli.param(
                Key(DocoptConfig, '--linker-stock'),
                Key(LinkerConfig, 'conc_uM'),
                default=10,
        )
        tag = None

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
    reaction_volume_uL = appcli.param(
            Key(DocoptConfig, '--total-volume', cast=float_or_expr),
            default=None,
    )
    mrna_volume_uL = appcli.param(
            Key(DocoptConfig, '--mrna-volume', cast=float_or_expr),
            default=None,
    )
    mrna_conc_uM = appcli.param(
            Key(DocoptConfig, '--mrna-conc', cast=float_or_expr),
            default=None,
    )
    linker_ratio = appcli.param(
            Key(DocoptConfig, '--linker-ratio', cast=float),
            default=1,
    )
    extra_anneal_uL = appcli.param(
            Key(DocoptConfig, '--extra-anneal-rxn', cast=float_or_expr),
            default=0,
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
        # If I load this from a preset, I'll want to make it so that the water 
        # is actually a solvent.  See stepwise issue #54.
        rxn = stepwise.MasterMix("""\
            Reagent                 Stock  Flags        Volume  MM?
            ====================  =======  =======  ==========  ===
            nuclease-free water            ligate       3.8 µL   +
            T4 RNA ligase buffer      10x  buffer       1.0 µL   +
            ATP                     10 mM  ligate       1.0 µL   +
            mRNA                    10 µM  anneal       2.0 µL
            linker                  10 µM  anneal       2.0 µL
            T4 RNA ligase         10 U/µL  ligate       0.2 µL   +
        """)
        rxn.num_reactions = n = self.num_reactions
        rxn.extra_min_volume = 0.5, 'µL'

        k = self.linker_ratio
        c1 = Quantity(self.mrna_conc_uM, 'µM') or rxn['mRNA'].conc
        s1 = Quantity(min(x.mrna_stock_uM for x in self.reactions), 'µM')
        s2 = Quantity(min(x.linker_stock_uM for x in self.reactions), 'µM')
        v_free = sum(
            rxn[k].volume 
            for k in ['nuclease-free water', 'mRNA', 'linker']
        )
        v1_ideal = rxn.volume * (c1 / s1)
        v1_max = v_free * (s2 / (k * s1 + s2))
        v1 = min(v1_ideal, v1_max)
        v2 = v1 * (k * s1 / s2)

        rxn['mRNA'].stock_conc = s1
        rxn['linker'].stock_conc = s2

        rxn['mRNA'].volume = v1
        rxn['linker'].volume = v2
        rxn['nuclease-free water'].volume = v_free - v1 - v2

        rxn['mRNA'].master_mix = all_equal(x.mrna for x in self.reactions)
        rxn['linker'].master_mix = all_equal(x.linker for x in self.reactions)

        rxn['mRNA'].name = merge_names(x.mrna for x in self.reactions)
        rxn['linker'].name = merge_names(x.linker for x in self.reactions)

        if rxn['mRNA'].volume < v1_ideal:
            warn(f"cannot reach {c1} {rxn['mRNA'].name}, reducing volume from {v1_ideal} to {v1}.")

        if uL := self.mrna_volume_uL:
            rxn.hold_ratios.volume *= (uL, 'µL') / v1
        if uL := self.reaction_volume_uL:
            rxn.hold_ratios.volume = uL, 'µL'

        return rxn

    def get_annealing_reaction_without_extra(self):
        rxn = self.combined_reaction.copy()
        v_reagents = 0

        for reagent in rxn:
            if 'ligate' in reagent.flags:
                del rxn[reagent.key]
            if 'anneal' in reagent.flags:
                v_reagents += reagent.volume

        rxn['T4 RNA ligase buffer'].volume = v_reagents / 9
        return rxn

    def get_annealing_reaction(self):
        rxn = self.annealing_reaction_without_extra.copy()
        rxn.hold_ratios.volume += self.extra_anneal_uL, 'µL'
        return rxn

    def get_ligation_reaction(self):
        rxn = self.combined_reaction.copy()
        anneal_rxn = self.annealing_reaction_without_extra

        for reagent in rxn:
            if 'anneal' in reagent.flags:
                del rxn[reagent.key]

        rxn['T4 RNA ligase buffer'].volume -= \
                anneal_rxn['T4 RNA ligase buffer'].volume
        rxn['annealed mRNA/linker'].volume = \
                anneal_rxn.volume
        rxn['annealed mRNA/linker'].master_mix = \
                anneal_rxn['mRNA'].master_mix and \
                anneal_rxn['linker'].master_mix

        return rxn

if __name__ == '__main__':
    LigateMrnaLinker.main()
