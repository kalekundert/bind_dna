#!/usr/bin/env python3

import stepwise
import docopt
import appcli
import autoprop
import po4

from appcli import DocoptConfig, Key
from inform import plural

@autoprop
class AnnealMrnaLinker(appcli.App):
    """\
Anneal linker-N and mRNA prior to ligation.

Usage:
    anneal [<mrna>] [<linker>] [-n <int>] [-v <µL>] [-m <reagents>] [-x <fold>] 
        [-R <µM>] [-L <conc>]

Arguments:
    <mrna>
        The name of the mRNA, e.g. f11.  Multiple comma-separated names may be 
        given.

    <linker>
        The name of the linker, e.g. o93.  Multiple comma-separated names may 
        be given.

Options:
    -n --num-reactions <int>
        The number of reactions to setup up.  The default is the number of 
        unique combinations of mRNA and linker.

    -v --volume <µL>                    [default: ${app.volume_uL}]
        The volume of each annealing reaction in µL.

    -m --master-mix <reagents>          [default: ${','.join(app.master_mix)}]
        A comma-separated list of reagents to include in the master mix.  This 
        flag is only relevant in <n> is more than 1.  The following reagents 
        are understood: mrna, link

    -x --excess-linker <fold>           [default: ${app.excess_linker}]
        The amount of linker to add to the reaction, relative to the amount of 
        mRNA in the reaction.

    -R --mrna-stock <µM>
        The stock concentration of the mRNA, in µM.  The volume of mRNA will be 
        updated accordingly to keep the amount of material in the reaction 
        constant.  The default is read from the PO₄ database, or 10 µM if the 
        given mRNA is not in the database.

    -L --linker-stock <conc>
        The stock concentration of the linker, in user-specified units.  The 
        volume of linker will not be updated, so this will change the relative 
        proportion of linker to mRNA.
"""
    __config__ = [
            DocoptConfig(),
    ]

    mrnas = appcli.param(
            '<mrna>',
            cast=lambda x: x.split(','),
            default_factory=list,
    )
    linkers = appcli.param(
            '<linker>',
            cast=lambda x: x.split(','),
            default_factory=list,
    )
    num_reactions = appcli.param(
            '--num-reactions',
            cast=int,
            default=0,
            get=lambda self, x: x or len(self.mrnas) * len(self.linkers),
    )
    volume_uL = appcli.param(
            '--volume',
            cast=eval,
            default=4,
    )
    master_mix = appcli.param(
            '--master-mix',
            cast=lambda x: set(x.split(',')),
            default_factory=set,
    )
    excess_linker = appcli.param(
            '--excess-linker',
            cast=float,
            default=1,
    )
    mrna_stock_uM = appcli.param(
            '--mrna-stock',
            default=None,
    )
    linker_stock_uM = appcli.param(
            '--mrna-stock',
            default=None,
    )

    def get_protocol(self):
        p = stepwise.Protocol()
        rxn = self.reaction
        n = rxn.num_reactions

        p += stepwise.Step(
                f"Setup {plural(n):# annealing reaction/s} [1]:",
                rxn,
        )

        p.footnotes[1] = stepwise.Footnote("""\
                Using 0.6x linker reduces the amount of unligated 
                linker, see expt #1."""
        )

        p += stepwise.Step(
                f"Perform the {plural(n):annealing reaction/s}:",
                substeps=[
                    "Incubate at 95°C for 2 min.",
                    "Cool at room temperature.",
                ],
                br='\n',
        )

        return p

    def get_reaction(self):
        rxn = stepwise.MasterMix.from_text("""\
                Reagent              Stock     Volume  MM?
                ===================  =====  =========  ===
                nuclease-free water         to 4.0 µL   +
                PBS                  10x       0.4 µL   +
                mRNA                 10 µM     0.5 µL   -
                linker               10 µM     0.5 µL   -
        """)

        rxn.num_reactions = n = self.num_reactions
        rxn.hold_ratios.volume = self.volume_uL, 'µL'
        rxn['mRNA'].master_mix = 'mrna' in self.master_mix
        rxn['linker'].master_mix = 'link' in self.master_mix
    
        if self.mrnas:
            rxn['mRNA'].name = ','.join(self.mrnas)
        if self.linkers:
            rxn['linker'].name = ','.join(self.linkers)

        db = po4.load_db()

        rxn['mRNA'].hold_conc.stock_conc = consensus(
                get_conc_uM(db, x, self.mrna_stock_uM)
                for x in self.mrnas
        )
        rxn['linker'].hold_conc.stock_conc = consensus(
                get_conc_uM(db, x, self.linker_stock_uM)
                for x in self.linkers
        )
        rxn['linker'].volume *= self.excess_linker

        return rxn

def get_conc_uM(db, tag, override):
    if override:
        return float(override), 'µM'
    try:
        return db[tag].conc_nM / 1000, 'µM'
    except po4.QueryError:
        return 10, 'µM'

def consensus(values):
    from itertools import groupby
    from more_itertools import one
    return one((x[0] for x in groupby(values)))


if __name__ == '__main__':
    app = AnnealMrnaLinker.from_params()
    appcli.load(app)
    app.protocol.print()


