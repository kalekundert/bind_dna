#!/usr/bin/env python3

"""\
Setup in vitro transcription/translation reactions using the NEB PURExpress 
system (E6800).

Usage:
    purespress.py <num_rxns>
"""

import docopt
import dirty_water
from nonstdlib import plural

args = docopt.docopt(__doc__)
protocol = dirty_water.Protocol()
purexpress = dirty_water.Reaction('''\
Reagent               Conc  Each Rxn  Master Mix
===============  =========  ========  ==========
water                         1.5 μL         yes
A                             4.0 μL         yes
B                             3.0 μL         yes
RNase Inhibitor    40 U/μL    0.2 μL         yes
ZnOAc                 1 mM    0.5 μL         yes
DNA                  75 nM    0.8 μL
''')

purexpress.num_reactions = int(args['<num_rxns>'])
purexpress.show_master_mix = True

protocol += f"""\
Setup {plural(purexpress.num_reactions):? IVTT reaction/s}:

{purexpress}

- Keep on ice.
- Be sure to add A before B.
- The control template (125 ng/μL) is 75 nM."""

protocol += """\
Incubate at 37°C for 2h."""

print(protocol)
