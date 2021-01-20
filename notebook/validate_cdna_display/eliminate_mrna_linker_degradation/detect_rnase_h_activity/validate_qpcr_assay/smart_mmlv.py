#!/usr/bin/env python3

import stepwise, appcli, autoprop
from po4 import load_db
from appcli import Key, DocoptConfig
from inform import plural
from copy import deepcopy
from stepwise_mol_bio import Main, comma_list, comma_set
from operator import not_

@autoprop
class SmartMmlv(Main):
    """\
Synthesize first-strand cDNA (up to 11.7 kb) using SMART MMLV reverse 
transcriptase.

Usage:
    smart_mmlv <templates> <primers> [-n <rxns>] [-m <reagents>] [-v <µL>]
        [-t <pmol>] [-T <nM>] [-p <fold>] [-P <µM>] [-C] [-q <step>] 

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

    -v --volume <µL>                [default: ${app.volume_uL}]
        The volume of the reaction, in µL.

    -t --template-pmol <pmol>       [default: ${app.template_pmol}]
        The amount of template to use, in pmol.  Note that this value is for a 
        20 µL reaction, and will be scaled if the `--volume` option is given.

    -T --template-conc <nM>         [default: ${app.template_conc_nM}]
        The stock concentration of the template, in nM.

    -p --excess-primer <fold>       [default: ${app.primer_excess}]
        The amount of primer to use, relative to the template concentration.

    -P --primer-conc <µM>           [default: ${app.primer_conc_uM}]
        The stock concentration of the primer, in µM.

    -C --no-control
        Exclude the control reaction with no reverse transcriptase.  This is a 
        standard control in RT-qPCR workflows.

    -q --quench <step>      [default: ${app.quench}]
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
    volume_uL = appcli.param(
            '--volume',
            cast=float,
            default=20,
    )
    template_pmol = appcli.param(
            '--template-pmol',
            cast=float,
            default=1,
    )
    template_conc_nM = appcli.param(
            '--template-conc',
            cast=float,
            default=1000,
    )
    primer_excess = appcli.param(
            '--excess-primer',
            cast=float,
            default=50,
    )
    primer_conc_uM = appcli.param(
            '--primer-conc',
            cast=float,
            default=100,
    )
    master_mix = appcli.param(
            '--master-mix',
            cast=comma_set,
            default_factory=set,
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

        p += f"""\
Anneal the {plural(self.primers):RT primer/s} to the {plural(self.templates):RNA template/s} [1]:

{anneal}
"""
        p.footnotes[1] = """\
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
"""
        p += f"""\
Incubate at 70°C for 3 min, then immediately
cool on ice.
"""
        p += f"""\
Setup {plural(mmlv.num_reactions):# reverse transcription reaction/s}:

{mmlv}
"""
        if self.nrt_control:
            p += f"""\
Setup a −reverse transcriptase (NRT) control:

{mmlv_nrt}
"""

        p += """\
Incubate at 42°C for 60 min [2].
"""
        p.footnotes[2] = """\
Samples can be incubated for 50-90 min if
necessary.
"""
        if self.quench == 'none':
            pass

        elif self.quench == 'heat':
            p += """\
Incubate at 70°C for 15 min.
"""
        elif self.quench == 'edta':
            p += """\
Add 4 µL 60 mM EDTA.
"""
        else:
            raise ConfigError(f"unexpected value for quench parameter: {self.quench}")

        return p
    
    def get_reactions(self):
        anneal = stepwise.MasterMix()
        anneal.num_reactions = self.num_reactions
        anneal.solvent = None
        anneal.extra_min_volume = 0.5, 'µL'

        anneal['template'].name = ','.join(self.templates)
        anneal['template'].stock_conc = self.template_conc_nM, 'nM'
        anneal['template'].volume = 1e3 * self.template_pmol / self.template_conc_nM, 'µL'
        anneal['template'].master_mix = 'rna' in self.master_mix
        anneal['template'].order = 2

        anneal['primer'].name = ','.join(self.primers)
        anneal['primer'].stock_conc = self.primer_conc_uM, 'µM'
        anneal['primer'].volume = self.primer_excess * self.template_pmol / self.primer_conc_uM, 'µL'
        anneal['primer'].master_mix = 'primer' in self.master_mix
        anneal['primer'].order = 3

        anneal['first-strand buffer'].stock_conc = '5x'
        anneal['first-strand buffer'].volume = 0, 'µL'
        anneal['first-strand buffer'].volume = anneal.volume / 4
        anneal['first-strand buffer'].master_mix = bool(self.master_mix)
        anneal['first-strand buffer'].order = 1

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
        mmlv['first-strand buffer'].volume -= anneal['first-strand buffer'].volume
        mmlv['annealed template/primer'].volume = anneal.volume
        mmlv['annealed template/primer'].stock_conc = \
                anneal['template'].stock_conc * (
                        anneal['template'].volume / anneal.volume)

        # Scale the volume of the RT reaction as requested, and scale the 
        # volume of the annealing reaction accordingly, but not so much that 
        # we'd need to pipet volumes smaller than 0.5 µL:

        mmlv.hold_ratios.volume = self.volume_uL, 'µL'
        anneal.hold_ratios.volume = mmlv['annealed template/primer'].volume * 1.2

        pipet_volumes = [
                anneal.master_mix_volume,
        ]
        pipet_volumes += [
                anneal[k].volume
                for k in ['template', 'primer']
                if not anneal[k].master_mix
        ]
        anneal.hold_ratios.volume *= max(
                1,
                *('0.5 µL' / x for x in pipet_volumes),
        )

        return anneal, mmlv

if __name__ == '__main__':
    SmartMmlv.main()
