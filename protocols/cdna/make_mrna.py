#!/usr/bin/env python3

import stepwise
import appcli
import autoprop

from anneal_mrna_linker import AnnealMrnaLinker
from ligate_linker_n import LigateMrnaLinker
from wash_barendt import WashBarendt
from operator import not_

@autoprop
class MakeMrnaLinker(appcli.App):
    """\
Ligate a puromycin linker to mRNA.

Usage:
    make_mrna.py <mrna> <linker> [-n <rxns>] [-v <µL>] [options]

Arguments:
    <mrna>
        A comma-separated list of mRNA names (e.g. f85).  If possible, the 
        concentrations of these mRNAs will be read from the PO₄ database.

    <linker>
        A comma-separated list of linker names (e.g. o129).  If possible, the 
        concentrations of these linkers will be read from the PO₄ database.

Options:
    -n --num-reactions <int>        [default: ${app.num_reactions}]
        The number of separate reactions to set up.  By default, this is 
        inferred from the number of mRNAs and linkers specified.

    -c --target-conc <µM>           [default: ${app.target_conc_uM}]
        The concentration of purified mRNA you'd like to achieve.  Note that 
        you'll have to actually measure the concentration of the mRNA to know 
        how much to dilute it (see `sw dilute`).  This option just scales the 
        annealing and ligation reactions such that the result should be more 
        concentrated than the target.

    -v --mrna-volume <µL>
        The volume of mRNA to use, e.g. if you have an entire aliquot you want 
        to use.  This can only increase the amount of mRNA that will be used; 
        if you want to use less mRNA you may also need to lower the target 
        concentration (--target-conc).

    -m --mwco <kDa>                 [default: ${app.mwco_kDa}]
        The MWCO of the spin filter.  According to Millipore Sigma, this should 
        be 2-3x smaller than the molecular weight of the ligated product: 
        https://tinyurl.com/4ffxu8zb

    -y --expected-yield <percent>   [default: ${100 * app.expected_yield}]
        The percentage of the mRNA added to the reaction that you expect to 
        recover from the spin filter.  This is used to ensure that you start 
        with enough mRNA to reach the desired final concentration.

    -d --expected-dead-volume <µL>  [default: ${app.expected_dead_volume_uL}]
        The volume of purified mRNA that you expect to recover from the spin 
        column.  The 500 µL Amicon spin filters have a minimum dead volume of 
        15 µL, but I usually recover somewhat more than that.

    -l --label <name>
        The label that should be used for the ligated and purified product.
        
    -W --no-wash
        Skip the wash step.

    -A --no-aliquot
        Skip the aliquot step.
"""
    __config__ = [
            appcli.DocoptConfig(),
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
    mwco_kDa = appcli.param(
            '--mwco',
            cast=int,
            default=100,
    )
    target_conc_uM = appcli.param(
            '--target-conc',
            cast=float,
            default=1,
    )
    mrna_volume_uL = appcli.param(
            '--mrna-volume',
            cast=eval,
            default=0,
    )
    expected_yield = appcli.param(
            '--expected-yield',
            cast=lambda x: float(x) / 100,
            default=0.5,
    )
    expected_dead_volume_uL = appcli.param(
            '--expected-dead-volume',
            cast=float,
            default=25,
    )
    label = appcli.param(
            '--label',
            default=None,
    )
    wash = appcli.param(
            '--no-wash',
            cast=not_,
            default=True,
    )
    aliquot = appcli.param(
            '--no-aliquot',
            cast=not_,
            default=True,
    )

    def get_protocol(self):
        anneal = AnnealMrnaLinker(self.mrnas, self.linkers)
        anneal.linker_ratio = 0.6  # See expt #1

        # Scale the reactions such that we'll end up with enough mRNA:
        anneal_rxn = anneal.reaction

        mrna_pmol = (
                self.expected_dead_volume_uL * self.target_conc_uM /
                self.expected_yield
        )
        mrna_conc_uM = anneal_rxn['mRNA'].stock_conc.value

        anneal.mrna_volume_uL = max(
                mrna_pmol / mrna_conc_uM,
                self.mrna_volume_uL,
        )

        ligate = LigateMrnaLinker(anneal)
        wash = WashBarendt()
        wash.mwco_kDa = self.mwco_kDa

        p = stepwise.Protocol()
        p += anneal.protocol
        p += ligate.protocol

        if self.label:
            p += f"Label the product: {self.label}"
        if self.wash:
            p += wash.protocol
        if self.aliquot:
            p += stepwise.load('aliquot "4 µL" "1 µM"')

        return p

if __name__ == '__main__':
    app = MakeMrnaLinker.from_params()
    appcli.load(app)
    app.protocol.print()
