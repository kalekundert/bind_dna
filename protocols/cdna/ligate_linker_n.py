#!/usr/bin/env python3

import stepwise
import appcli
import autoprop

from appcli import DocoptConfig, Key
from stepwise import pl, ul
from anneal_mrna_linker import AnnealMrnaLinker
from inform import plural
from operator import not_

# I can't remember where exactly this protocol came from.  Some thoughts:
# 
# - 10x diluted PBS: my own optimizations of the annealing reaction.
# 
# - 20 U T4 RNA ligase: Probably [Naimudden2016]_, even though I'm using 10x 
#   less mRNA/linker (I probably assumed that more enzyme wouldn't hurt, and 
#   this is already about as little as I can pipet).
# 
# - BSA: Probably because Takara included BSA with the enzyme, although I'm 
#   using a different concentration that their online protocol suggests.
# 
# - 10 min incubation time at 25°C: [Naimudden2016]_.  I don't know where the 
#   65°C incubation/heat denaturation came from, though.

@autoprop
class LigateMrnaLinker(appcli.App):
    """\
Ligate linker-N to the mRNA.

Usage:
    ligate <mrna_µL> <mrna_µM> [-n <int>] [-v <µL>] [-x <percent>] [-i <time>] [options]

Arguments:
    <mrna_µL>
        The volume of the annealed mRNA, in µL.  The ligation reaction will 
        be 10x this volume, to dilute the salt from the annealing reaction.

    <mrna_µM>
        The concentration of the mRNA, in µM.  The amount of ligase will be 
        scaled relative to the quantity of mRNA to be ligated.

Options:
    -n --num-reactions <int>        [default: ${app.num_reactions}]
        The number of reactions to setup.

    -x --extra <percent>            [default: ${app.extra_percent}]
        How much extra master mix to prepare.

    -i --incubate <time>            [default: ${app.incubate_time}]
        How long to incubate the reaction at the temperature indicated by the 
        `-t` flag.  Include a unit.

    -t --incubate-temp <temp>       [default: ${app.incubate_temp}]
        What temperature to incubate the reaction at.  Include a unit.

    -m --master-mix <reagents>      [default: ${','.join(app.master_mix)}]
        Include the indicated reagents in the master mix.  The following 
        reagents are understood:

        pnk: T4 PNK
        lig: T4 RNA ligase
        rna: annealed mRNA/linker
        peg: PEG-8000

    -M --no-master-mix
        Exclude all optional reagents from the master mix, i.e. `-m ''`.

    -L --no-ligase
        Remove the ligase from the reaction, e.g. as a negative control.

    -k --kinase
        Add T4 PNK to the reaction, e.g. if using non-phosphorylated primers.

    -p --peg
        Include PEG-8000 in the reaction.  Many T4 RNA ligase protocols 
        recommend this, but in my hands it does not improve yield.  This may be 
        because my substrates are already annealed.

    -Q --no-quench
        Leave out the 65°C incubation to quench the reaction.  This is useful 
        if the reaction will be quenched by a downstream step anyways.

    -I --no-incubate
        Skip the entire incubation step.  This is useful when setting up 
        multiple ligate reaction in a row; only the last needs include the 
        incubation.
"""
    __config__ = [
            DocoptConfig(),
    ]

    num_reactions = appcli.param(
            '--num-reactions',
            cast=int,
            default=1,
    )
    mrna_volume_uL = appcli.param(
            '<mrna_µL>',
            cast=float,
    )
    mrna_conc_uM = appcli.param(
            '<mrna_µM>',
            cast=float,
    )
    extra_percent = appcli.param(
            '--extra',
            cast=float,
            default=10,
    )
    incubate_time = appcli.param(
            '--incubate',
            default='10 min',
    )
    incubate_temp = appcli.param(
            '--incubate-temp',
            default='25°C',
    )
    master_mix = appcli.param(
            Key(DocoptConfig, '--no-master-mix', cast=lambda x: set()),
            Key(DocoptConfig, '--master-mix', cast=lambda x: set(x.split(','))),
            default_factory=lambda: {'peg','lig'},
    )
    use_ligase = appcli.param(
            '--no-ligase',
            cast=not_,
            default=True,
    )
    use_kinase = appcli.param(
            '--kinase',
            default=False,
    )
    use_peg = appcli.param(
            '--peg',
            default=False,
    )
    quench = appcli.param(
            '--no-quench',
            cast=not_,
            default=True,
    )
    incubate = appcli.param(
            '--no-incubate',
            cast=not_,
            default=True,
    )

    def __init__(self, anneal: AnnealMrnaLinker):
        anneal_rxn = anneal.reaction
        self.mrna_volume_uL = anneal_rxn.volume.value
        self.mrna_conc_uM = anneal_rxn['mRNA'].conc.value

    def get_protocol(self):
        p = stepwise.Protocol()
        rxn = self.reaction
        rxn_name = 'ligation' if self.use_ligase else 'negative control'
        n = rxn.num_reactions

        p += stepwise.pl(
                f"Setup {plural(n):# {rxn_name} reaction/s}:",
                rxn,
        )
        if self.incubate:
            p += pl(
                    f"Incubate the {plural(n):ligation reaction/s} as follows:",
                    s := ul(
                        f"{self.incubate_temp} for {self.incubate_time}."
                    ),
            )
            if self.quench:
                s += "65°C for 10 min."

        return p

    def get_reaction(self):
        return self.reaction_neb

    def get_reaction_takara(self):
        # Unit definition:
        # - One unit is defined as the amount of enzyme that converts 1 pmol of 
        #   [5'-32P]pCp into its acid-insoluble form in 10 minutes at 5°C, 
        #   using oligo(A) as the substrate during 3' end labeling of RNA.
        # - 2x excess relative to above definition.
        # - https://tinyurl.com/3lhzf7a6
        rxn = stepwise.MasterMix.from_text("""\
                Reagent                  Stock       Volume  MM?
                =====================  =======  ===========  ===
                nuclease-free water             to 40.00 µL   +
                T4 DNA ligase buffer       10x       4.0 µL   +
                BSA                       0.1%       4.0 µL   +
                PEG 8000                   50%      20.0 µL   +
                T4 PNK                 10 U/µL      0.33 µL   -
                T4 RNA ligase          40 U/µL      0.25 µL   -
                annealed mRNA/linker   1.25 µM       4.0 µL   -
        """)
        rxn['T4 RNA ligase'].name = "T4 RNA ligase (Takara 2050)"
        return self._adjust_reaction(rxn)

    def get_reaction_neb(self):
        # Unit definition:
        # - One unit is defined as the amount of enzyme required to convert 1 
        #   nanomole of 5´-[32P]rA16 into a phosphatase-resistant form in 30 
        #   minutes at 37°C.
        #
        #   - I'm going to use the same volume of NEB enzyme as I did of Takara 
        #     enzyme.
        #   - This unit definition has 1000x more RNA, but a longer reaction 
        #     time and a warmer reaction temperature.
        #   - The NEB enzyme is also 10 U/µL, which the Takara one is 40 U/µL.
        #   - I'm not sure how these effects cancel out, so I'm just going to 
        #     assume that the same volume will be about right.
        #
        # - The NEB ligation protocol calls for:
        #   - 10 U T4 RNA ligase
        #   - 20 pmol RNA
        #   - 2h incubation at 25°C
        #   - https://www.neb.com/protocols/2018/10/17/protocol-ligation-of-an-oligo-to-the-3-end-of-rna-using-t4-rna-ligase-1m0204
        #
        # - My oligos are annealed, though, which should significantly speed up 
        #   the reaction.
        rxn = stepwise.MasterMix("""\
                Reagent                  Stock       Volume  MM?
                =====================  =======  ===========  ===
                nuclease-free water             to 40.00 µL   +
                T4 RNA ligase buffer       10x       4.0 µL   +
                ATP                      10 mM       4.0 µL   +
                PEG 8000                   50%      20.0 µL   +
                T4 PNK                 10 U/µL      0.33 µL   -
                T4 RNA ligase          10 U/µL      0.25 µL   -
                annealed mRNA/linker   1.25 µM       4.0 µL   -
        """)
        rxn['T4 RNA ligase'].name = "T4 RNA ligase (NEB M0204)"
        return self._adjust_reaction(rxn)

    def _adjust_reaction(self, rxn):
        v = self.mrna_volume_uL / rxn['annealed mRNA/linker'].volume.value
        c = self.mrna_conc_uM / rxn['annealed mRNA/linker'].stock_conc.value

        rxn.num_reactions = self.num_reactions
        rxn.extra_percent = self.extra_percent
        rxn.extra_min_volume = '0.5 µL'
        rxn.hold_ratios.volume *= v
        rxn['PEG 8000'].master_mix = 'peg' in self.master_mix
        rxn['T4 PNK'].volume *= c
        rxn['T4 PNK'].master_mix = 'pnk' in self.master_mix
        rxn['T4 RNA ligase'].volume *= c
        rxn['T4 RNA ligase'].master_mix = 'lig' in self.master_mix
        rxn['annealed mRNA/linker'].volume = self.mrna_volume_uL, 'µL'
        rxn['annealed mRNA/linker'].stock_conc = self.mrna_conc_uM, 'µM'
        rxn['annealed mRNA/linker'].master_mix = 'rna' in self.master_mix

        if not self.use_ligase:
            del rxn['T4 RNA ligase']
        if not self.use_kinase:
            del rxn['T4 PNK']
        if not self.use_peg:
            del rxn['PEG 8000']

        return rxn



if __name__ == '__main__':
    app = LigateMrnaLinker.from_params()
    appcli.load(app)
    app.protocol.print()
