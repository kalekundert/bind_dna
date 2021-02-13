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
    anneal [<mrna>] [<linker>] [-n <int>] [-m <reagents>] [-V <µL>]
            [-r <µL>] [-R <µM>] [-l <ratio>] [-L <µM>]

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

    -m --master-mix <reagents>          [default: ${','.join(app.master_mix)}]
        A comma-separated list of reagents to include in the master mix.  This 
        flag is only relevant in <n> is more than 1.  The following reagents 
        are understood: mrna, link

    -r --mrna-volume <µL>               [default: ${app.mrna_volume_uL}]
        The volume of mRNA to use in each annealing reaction, in µL.

    -R --mrna-stock <µM>
        The stock concentration of the mRNA, in µM.  The volume of mRNA will be 
        updated accordingly to keep the amount of material in the reaction 
        constant.  The default is read from the PO₄ database, or 10 µM if the 
        given mRNA is not in the database.  Use '--mrna-volume' to change the 
        amount of mRNA in the reaction.

    -l --linker-ratio <float>           [default: ${app.linker_ratio}]
        The amount of linker to add to the reaction, relative to the amount of 
        mRNA in the reaction.

    -L --linker-stock <µM>
        The stock concentration of the linker, in µM.  The volume of linker be 
        updated accordingly, to keep the amount of material in the reaction 
        constant.  The default is read from the PO₄ database, or 10 µM if the 
        given linker is not in the database.  
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
    mrna_volume_uL = appcli.param(
            '--mrna-volume',
            cast=eval,
            default=1,
    )
    master_mix = appcli.param(
            '--master-mix',
            cast=lambda x: set(x.split(',')),
            default_factory=set,
    )
    linker_ratio = appcli.param(
            '--linker-ratio',
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

    def __bareinit__(self):
        self.db = po4.load_db()

    def __init__(self, mrnas, linkers):
        self.mrnas = list_if_str(mrnas)
        self.linkers = list_if_str(linkers)

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
        rxn = stepwise.MasterMix()
        rxn.num_reactions = n = self.num_reactions
        rxn.solvent = None

        rxn['mRNA'].stock_conc = consensus(
                get_conc_uM(self.db, x, self.mrna_stock_uM)
                for x in self.mrnas
        )
        rxn['linker'].stock_conc = consensus(
                get_conc_uM(self.db, x, self.linker_stock_uM)
                for x in self.linkers
        )

        rxn['mRNA'].volume = self.mrna_volume_uL, 'µL'
        rxn['linker'].volume = (
                self.linker_ratio
                * rxn['mRNA'].volume
                * (rxn['mRNA'].stock_conc / rxn['linker'].stock_conc)
        )

        rxn['mRNA'].master_mix = 'mrna' in self.master_mix
        rxn['linker'].master_mix = 'link' in self.master_mix
    
        if self.mrnas:
            rxn['mRNA'].name = ','.join(self.mrnas)
        if self.linkers:
            rxn['linker'].name = ','.join(self.linkers)

        rxn['PBS'].volume = rxn.volume / 9
        rxn['PBS'].order = -1

        return rxn

def list_if_str(x):
    return [x] if isinstance(x, str) else x

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


