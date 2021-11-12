#!/usr/bin/env python3

import stepwise, appcli, autoprop
from po4 import load_db
from appcli import Key, DocoptConfig
from inform import plural
from copy import deepcopy
from stepwise import pl, ul, pre
from stepwise_mol_bio import Main, comma_list, comma_set
from operator import not_

@autoprop
class SmartMmlv(Main):
    """\
Synthesize first-strand cDNA (up to 11.7 kb) using SMART MMLV reverse 
transcriptase.

Usage:
    smart_mmlv <templates> <primers> [-n <rxns>] [-m <reagents>] [-v <µL>]
        [-t <nM>] [-T <nM>] [-p <fold>] [-P <µM>] [-C] [-q <step>] 

Arguments:
    <templates>
        The names of the templates to reverse transcribe, comma-separated.

    <primers>
        The names of the primers to use, comma-separated.

Options:
    -n --num-reactions <rxns>
        The number of reactions to setup.  By default, this is the number of 
        templates specified.

    -m --master-mix <reagents>
        A comma-separated list of reagents to include in the annealing master 
        mix.  Understood reagents are 'rna' and 'primer'.  By default, each 
        reaction is assembled independently.

    -v --volume <µL>                    [default: ${app.volume_uL}]
        The volume of the reaction, in µL.

    -t --template-conc <nM>             [default: ${app.template_conc_nM}]
        The final concentration of the template in the reverse transcription 
        reaction.  Note that a higher concentration will be used in the 
        annealing reaction.

    -T --template-stock <nM>            [default: ${app.template_stock_nM}]
        The stock concentration of the template, in nM.

    -p --excess-primer <fold>           [default: ${app.primer_excess}]
        The amount of primer to use, relative to the template concentration.

    -P --primer-stock <µM>              [default: ${app.primer_stock_uM}]
        The stock concentration of the primer, in µM.

    -r --rna-min-volume <µL>            [default: ${app.rna_min_volume_uL}]
        The smallest volume of RNA to pipet in any step.  This is meant to help 
        avoid pipetting errors when trying to be quantitative, e.g. for qPCR.

    -a --extra-anneal <percent>         [default: ${app.extra_anneal_percent}]
        How much extra annealing reaction to prepare, e.g. to facilitate 
        multichannel pipetting.

    -C --no-control
        Exclude the control reaction with no reverse transcriptase.  This is a 
        standard control in RT-qPCR workflows.

    -q --quench <step>                  [default: ${app.quench}]
        How to quench the reaction:

        heat: Incubate at 70°C.
        edta: Add EDTA.
        none: Skip the quench step.

Ordering:
- SMART MMLV reverse transcriptase (Takara 639523)
- Advantage UltraPure PCR deoxynucleotide mix (Takara 639125)

References:
1. https://tinyurl.com/y4ash7dl
"""
    __config__ = [
            DocoptConfig(),
    ]

    templates = appcli.param(
            '<templates>',
            cast=comma_list,
    )
    primers = appcli.param(
            '<primers>',
            cast=comma_list,
    )
    num_reactions = appcli.param(
            '--num-reactions',
            cast=int,
            default=None,
            get=lambda self, n: n if n else len(self.templates),
    )
    master_mix = appcli.param(
            '--master-mix',
            cast=comma_set,
            default_factory=set,
    )
    volume_uL = appcli.param(
            '--volume',
            cast=float,
            default=20,
    )
    template_conc_nM = appcli.param(
            '--template-conc',
            cast=float,
            default=50,
    )
    template_stock_nM = appcli.param(
            '--template-stock',
            cast=float,
            default=1000,
    )
    primer_excess = appcli.param(
            '--excess-primer',
            cast=float,
            default=50,
    )
    primer_stock_uM = appcli.param(
            '--primer-conc',
            cast=float,
            default=100,
    )
    rna_min_volume_uL = appcli.param(
            '--rna-min-volume',
            cast=float,
            default=2,
    )
    extra_anneal_percent = appcli.param(
            '--extra-anneal',
            cast=float,
            default=20,
    )
    nrt_control = appcli.param(
            '--no-control',
            cast=not_,
            default=True,
    )
    quench = appcli.param(
            '--quench',
            default='heat',
    )

    def get_protocol(self):
        p = stepwise.Protocol()
        anneal, mmlv = self.reactions

        if self.nrt_control:
            mmlv_nrt = deepcopy(mmlv)
            mmlv_nrt.num_reactions = 1
            anneal.num_reactions += 1
            del mmlv_nrt['SMART MMLV RT']

        p += pl(
                f"Anneal the {plural(self.primers):RT primer/s} to the {plural(self.templates):RNA template/s} [1]:",
                anneal,
        )

        p.footnotes[1] = pre("""\
This protocol is based on the official Takara 
SMART MMLV reverse transcription protocol:

https://tinyurl.com/y4ash7dl

However, I made three modifications:

- I reduced the volume of the annealing step as 
  much as possible, the increase the 
  concentrations of the oligos.

- I added buffer to the annealing step, because 
  annealing works much better with some salt to 
  shield the backbone charge.  I don't actually 
  know how much salt is in the buffer, but some is 
  better than none.

- I included the MMLV in the transcription master 
  mix, because excluding it seemed too inaccurate.  
  The changes to the annealing step also mean that 
  the master mix has a 1x buffer concentration, so 
  I don't need to worry about the enzyme being 
  unhappy.
""")
        p += "Incubate at 70°C for 3 min, then immediately cool on ice."
        p += pl(
                "Setup {plural(mmlv.num_reactions):# reverse transcription reaction/s}:",
                mmlv,
        )
        if self.nrt_control:
            p += pl(
                    f"Setup a −reverse transcriptase (NRT) control:",
                    mmlv_nrt,
            )

        p += "Incubate at 42°C for 60 min [2]."
        p.footnotes[2] = "Samples can be incubated for 50-90 min if necessary."

        if self.quench == 'none':
            pass

        elif self.quench == 'heat':
            p += "Incubate at 70°C for 15 min."

        elif self.quench == 'edta':
            p += "Add 4 µL 60 mM EDTA."

        else:
            raise ConfigError(f"unexpected value for quench parameter: {self.quench}")

        return p
    
    def get_reactions(self):

        # Define the annealing reaction:

        anneal = stepwise.MasterMix()
        anneal.num_reactions = self.num_reactions
        anneal.solvent = None
        anneal.extra_min_volume = 0.5, 'µL'

        template_fmol = self.volume_uL * self.template_conc_nM

        anneal['template'].name = ','.join(self.templates)
        anneal['template'].stock_conc = self.template_stock_nM, 'nM'
        anneal['template'].volume = template_fmol / self.template_stock_nM, 'µL'
        anneal['template'].master_mix = 'rna' in self.master_mix
        anneal['template'].order = 2

        anneal['primer'].name = ','.join(self.primers)
        anneal['primer'].stock_conc = self.primer_stock_uM, 'µM'
        anneal['primer'].volume = self.primer_excess * template_fmol / self.primer_stock_uM / 1e3, 'µL'
        anneal['primer'].master_mix = 'primer' in self.master_mix
        anneal['primer'].order = 3

        # If necessary, dilute the annealing reaction such that it will be 
        # necessary to add the given volume to the RT reaction.  The purpose of 
        # this is to ensure that we can pipet this volume accurately, since we 
        # often want to be quantitative about how much RNA we have (e.g. qPCR).

        if anneal.volume < (self.rna_min_volume_uL, 'µL'):
            anneal.solvent = 'nuclease-free water'
            anneal.volume = self.rna_min_volume_uL, 'µL'

        anneal['first-strand buffer'].stock_conc = '5x'
        anneal['first-strand buffer'].volume = 0, 'µL'
        anneal['first-strand buffer'].volume = anneal.volume / 4
        anneal['first-strand buffer'].master_mix = bool(self.master_mix)
        anneal['first-strand buffer'].order = 1

        # Define the MMLV reaction:

        mmlv = stepwise.MasterMix("""\
                Reagent                      Stock      Volume  MM?
                ========================  ========  ==========  ===
                nuclease-free water                 to 20.0 µL  yes
                first-strand buffer             5x      4.0 µL  yes
                dNTP mix                     10 mM      2.0 µL  yes
                DTT                         100 mM      2.0 µL  yes
                SMART MMLV RT             200 U/µL      0.5 µL  yes
                annealed template/primer                0.0 µL
        """)
        mmlv.num_reactions = self.num_reactions
        mmlv.hold_ratios.volume = self.volume_uL, 'µL'
        mmlv['first-strand buffer'].volume -= anneal['first-strand buffer'].volume
        mmlv['annealed template/primer'].volume = anneal.volume
        mmlv['annealed template/primer'].stock_conc = \
                anneal['template'].stock_conc * (
                        anneal['template'].volume / anneal.volume)

        # Scale the volume of the annealing reaction to guarantee that none of 
        # the volume are too small to pipet accurately.

        min_pipet_volumes = {
                'template': (self.rna_min_volume_uL, 'µL'),
        }
        pipet_volumes = [
                (anneal.master_mix_volume, '0.5 µL'),
        ]
        pipet_volumes += [
                (anneal[k].volume, min_pipet_volumes.get(k, '0.5 µL'))
                for k in ['template', 'primer']
                if not anneal[k].master_mix
        ]
        anneal.hold_ratios.volume *= max(
                1 + self.extra_anneal_percent / 100,
                *(limit / curr for curr, limit in pipet_volumes),
        )

        return anneal, mmlv

if __name__ == '__main__':
    SmartMmlv.main()
