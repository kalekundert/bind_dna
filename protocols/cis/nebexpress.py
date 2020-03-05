#!/usr/bin/env python3

"""\
Express proteins from linear DNA templates using NEBExpress.

Usage:
    nebexpress.py <templates>

Arguments:
    <templates>
        Comma-separated list of templates.  The number of reactions will be 
        inferred from this list.
"""

import docopt
import stepwise
from inform import plural

args = docopt.docopt(__doc__)
templates = args['<templates>']
num_templates = len(templates.split(','))

rxn = stepwise.MasterMix.from_text("""\
Reagent                       Stock    Volume  MM?
========================  =========  ========  ===
water                                to 50 µL  yes
S30 extract                             12 µL  yes
synthesis buffer                 2x     25 µL  yes
T7 RNA polymerase                        1 µL  yes
RNase inhibitor (murine)                 1 µL  yes
linear DNA template       125 ng/µL      2 µL
GamS nuclease inhibitor                  1 µL  yes
""")
rxn.hold_ratios.volume = '10 µL'
rxn.num_reactions = num_templates
rxn['linear DNA template'].name = templates

protocol = stepwise.Protocol()

protocol += f"""\
Setup {plural(rxn.num_reactions):# NEBExpress reaction/s} [1]:

{rxn}

- Thaw all components on ice
- Mix the S30 extract and protein synthesis buffer 
  by gently vortexing.
"""

protocol += f"""\
Incubate at 37°C for 2-4 hours [2].
"""

protocol.footnotes[1] = """\
During the experimental setup, it is recommended 
to add the linear DNA template in the last step to 
allow GamS to bind and inhibit RecBCD exonuclease 
before RecBCD has a chance to act on the DNA. 
"""

protocol.footnotes[2] = """\
Additional incubation time (maximum 10 hours) at 
37°C may increase yield.
"""

print(protocol)

# vim: tw=50
