#!/usr/bin/env python3

"""\
Express proteins from linear DNA templates using NEBExpress.

Usage:
    nebexpress.py <templates> [options]

Arguments:
    <templates>
        Comma-separated list of templates.  The number of reactions will be 
        inferred from this list.

Options:
    -d --template-stock <nM>             [default: 75]
        The stock concentration of the template DNA, in units of nM.

    -t --incubation-time <time>         [default: 2-4 hours]
        The amount of time to incubate the reactions.

    -T --incubation-temperature <temp>  [default: 37°C]
        The temperature to incubate the reactions at.
"""

import docopt
import stepwise
from inform import plural

args = docopt.docopt(__doc__)
templates = args['<templates>']
num_templates = len(templates.split(','))

rxn = stepwise.MasterMix.from_text("""\
Reagent                            Stock    Volume  MM?
=============================  =========  ========  ===
water                                     to 50 µL  yes
S30 extract [2]                              12 µL  yes
synthesis buffer [2]                  2x     25 µL  yes
T7 RNA polymerase               450 U/µL      1 µL  yes
RNase inhibitor (murine)         40 U/µL      1 µL  yes
GamS nuclease inhibitor [2,3]  1.5 µg/µL      1 µL  yes
linear DNA template                75 nM      2 µL
""")
rxn.hold_ratios.volume = '10 µL'
rxn.num_reactions = num_templates
rxn['linear DNA template'].name = f'{templates} [3]'
rxn['linear DNA template'].hold_conc.stock_conc = (
        float(args['--template-stock']), 'nM')

protocol = stepwise.Protocol()

protocol += f"""\
Setup {plural(rxn.num_reactions):# NEBExpress reaction/s} [1]:

{rxn}

- Thaw all components on ice
- Mix the S30 extract and protein synthesis buffer 
  by gently vortexing.
"""

protocol += f"""\
Incubate at {args['--incubation-temperature']} for {args['--incubation-time']} [4].
"""

protocol.footnotes[1] = """\
During the experimental setup, it is recommended 
to add the linear DNA template in the last step to 
allow GamS to bind and inhibit RecBCD exonuclease 
before RecBCD has a chance to act on the DNA. 
"""

protocol.footnotes[2] = """\
Aliquot to avoid multiple freeze/thaw cycles.
"""

protocol.footnotes[3] = """\
Optimal concentration must be determined 
empirically for each template.
"""

protocol.footnotes[4] = """\
Additional incubation time (maximum 10 hours) at 
37°C may increase yield.
"""

print(protocol)

# vim: tw=50
