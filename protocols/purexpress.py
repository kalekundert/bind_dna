#!/usr/bin/env python3

"""\
Setup in vitro transcription/translation reactions using the NEB PURExpress 
system (E6800).

Usage:
    purexpress.py <num_rxns> [-v <uL>] [-t]

Options:
    -v --rxn-volume <uL>  [default: 10]
        The volume of each individual reaction (in μL).  NEB recommends 25 μL, 
        but I typically use 10 μL and get enough yield for routine experiments.

    -t --add-target
        Add target DNA to the PURExpress reaction.
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
template DNA         75 nM    0.8 μL
''')

if args['--add-target']:
    purexpress['target DNA'].std_stock_conc = 750, 'nM'
    purexpress['target DNA'].std_volume = 0.8, 'μL'
    purexpress['target DNA'].master_mix = True
    purexpress['water'].std_volume = 0.7, 'μL'

purexpress.num_reactions = eval(args['<num_rxns>'])
purexpress.volume = eval(args['--rxn-volume'])
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
