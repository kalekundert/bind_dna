#!/usr/bin/env python3

import stepwise, appcli, autoprop
from inform import fatal, warn, plural

@autoprop
class NebExpress(appcli.App):
    """\
Express proteins from linear DNA templates using NEBExpress.

Usage:
    nebexpress.py <templates>... [-v <µL>] [-n <rxns>] [-c <nM>] [-C <nM>]
        [-t <time>] [-T <°C>] [-r]

Arguments:
    <templates>
        The templates to express.  The number of reactions will be inferred 
        from this list.

Options:
    -v --volume <µL>                    [default: ${app.volume_uL}]
        The volume of the reaction in µL.

    -n --num-reactions <int>
        The number of reactions to set up.  By default, this is inferred from
        the number of templates.

    -c --template-conc <nM>
        The desired final concentration of template in the reaction.

    -C --template-stock <nM>
        The stock concentration of the template DNA, in units of nM.  If not 
        specified, a concentration will be queried from the PO₄ database.  In 
        this case, all templates must be in the database and must have 
        identical concentrations.

    -t --incubation-time <time>         [default: ${app.incubation_time}]
        The amount of time to incubate the reactions.  No unit is assumed, so 
        be sure to include one.

    -T --incubation-temperature <°C>    [default: ${app.incubation_temp_C}]
        The temperature to incubate the reactions at, in °C.

    -r --mrna
        Use mRNA as the template instead of DNA.
"""
    __config__ = [
            appcli.DocoptConfig(),
    ]

    templates = appcli.param(
            '<templates>',
    )
    volume_uL = appcli.param(
            '--volume',
            cast=eval,
            default=10,
    )
    num_reactions = appcli.param(
            '--num-reactions',
            cast=eval,
            default=None,
    )
    incubation_time = appcli.param(
            '--incubation-time',
            default='2-4 hours',
    )
    incubation_temp_C = appcli.param(
            '--incubation-temp',
            cast=float,
            default=37,
    )
    use_mrna = appcli.param(
            '--mrna',
            default=False,
    )

    @appcli.param(
            '--template-conc',
            cast=float,
            default=None,
    )
    def template_conc_nM(self, x):
        if x is not None:
            return x

        if self.use_mrna:
            warn("mRNA template concentrations must be empirically optimized.  The default value is just a plausible starting point recommended by NEB (3 µg/50 µL reaction, for the DHFR control).")
            # NEB recommends 1–5 µg mRNA
            # Control DHFR mRNA: 760 bp, MW=469647 Da
            # 3 µg / 50 µL = 128 nM
            return 125

        else:
            warn("DNA template concentrations must be empirically optimized.  The default value is just a plausible starting point recommended by NEB (250 ng/50 µL reaction, for the DHFR control).")
            # NEB recommends 250 ng
            # Control DHFR plasmid: 2727 bp, MW=1.69e6 Da
            # 250 ng / 50 µL = 3 nM
            return 3

    @appcli.param(
            '--template-stock',
            cast=float,
            default=None,
    )
    def template_stock_nM(self, x):
        if x is not None:
            return x

        fatal("Must specify a template stock concentration")


    def get_protocol(self):
        p = stepwise.Protocol()
        rxn = self.reaction

        p += f"""\
Setup {plural(self.num_reactions):# NEBExpress reaction/s} [1]:

{rxn}

- Thaw all components on ice.
- Mix the S30 extract and protein synthesis buffer 
  by gently vortexing.
"""

        p += f"""\
Incubate at {self.incubation_temp_C}°C for {self.incubation_time} [4].
"""

        p.footnotes[1] = """\
During the experimental setup, it is recommended 
to add the linear DNA template in the last step to 
allow GamS to bind and inhibit RecBCD exonuclease 
before RecBCD has a chance to act on the DNA. 
"""

        p.footnotes[2] = """\
Aliquot to avoid multiple freeze/thaw cycles.
"""

        p.footnotes[3] = """\
Optimal concentration must be determined 
empirically for each template.
"""

        p.footnotes[4] = """\
Additional incubation time (maximum 10 hours) at 
37°C may increase yield.
"""

        return p

    def get_reaction(self):
        rxn = stepwise.MasterMix("""\
                Reagent                            Stock    Volume  MM?
                =============================  =========  ========  ===
                water                                     to 50 µL   +
                S30 extract [2]                              12 µL   +
                synthesis buffer [2]                  2x     25 µL   +
                T7 RNA polymerase               450 U/µL      1 µL   +
                RNase inhibitor (murine)         40 U/µL      1 µL   +
                GamS nuclease inhibitor [2,3]  1.5 µg/µL      1 µL   +
        """)
        rxn.hold_ratios.volume = self.volume_uL, 'µL'
        rxn.num_reactions = self.num_reactions or len(self.templates)

        rxn['template'].name = f"{','.join(self.templates)} [3]"
        rxn['template'].stock_conc = self.template_stock_nM, 'nM'
        rxn['template'].hold_stock_conc.conc = self.template_conc_nM, 'nM'

        rxn.hold_ratios.volume =  self.volume_uL, 'µL'

        if self.use_mrna:
            del rxn['T7 RNA polymerase']

        return rxn

if __name__ == '__main__':
    app = NebExpress.from_params()
    app.load()
    app.protocol.print()
