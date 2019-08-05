#!/usr/bin/env python3

"""\
Setup in vitro transcription/translation reactions using the NEB PURExpress 
system (E6800).

Usage:
    purexpress.py <num_rxns> [-v <uL>] [-t] [-z] [-p]

Options:
    -v --rxn-volume <uL>  [default: 10]
        The volume of each individual reaction (in μL).  NEB recommends 25 μL, 
        but I typically use 10 μL and get enough yield for routine experiments.

    -t --add-target
        Add target DNA to the PURExpress reaction.

    -z --add-zinc
        Add ZnOAc to the PURExpress reaction.  This is necessary when 
        expressing Zn-finger proteins.

    -p --purify
        Purify the expressed protein using the reverse-His protocol recommended 
        by NEB.
"""

import docopt
import dirty_water
from nonstdlib import plural

args = docopt.docopt(__doc__)
protocol = dirty_water.Protocol()
purexpress = dirty_water.Reaction('''\
Reagent               Conc  Each Rxn  Master Mix
===============  =========  ========  ==========
water                         2.0 μL         yes
A                             4.0 μL         yes
B                             3.0 μL         yes
RNase Inhibitor    40 U/μL    0.2 μL         yes
''')

if args['--add-zinc']:
    purexpress['ZnOAc'].std_stock_conc = 1, 'mM'
    purexpress['ZnOAc'].std_volume = 0.5, 'μL'
    purexpress['ZnOAc'].master_mix = True
    purexpress['water'].std_volume = purexpress['water'].std_volume - 0.5, 'μL'

if args['--add-target']:
    purexpress['target DNA'].std_stock_conc = 750, 'nM'
    purexpress['target DNA'].std_volume = 0.8, 'μL'
    purexpress['target DNA'].master_mix = True
    purexpress['water'].std_volume = purexpress['water'].std_volume - 0.8, 'μL'

purexpress['template DNA'].std_stock_conc = 75, 'nM'
purexpress['template DNA'].std_volume = 0.8, 'μL'
purexpress['template DNA'].master_mix = False

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

if args['--purify']:

    protocol += """\
Dilute reaction to 100 μL with PBS [1]."""

    protocol.notes += """\
We recommend a minimum volume of 100 μl after 
dilution to minimize losses during purification. 
If the target will be too dilute after addition of 
the diluent, we suggest a larger reaction volume 
be used.

Use of concentrated NaCl to dilute the reaction 
may help dissociate complexes between the target 
protein and translation factors.  The NaCl will 
remain, however, after the final elution and 
downstream applications may require microdialysis.  
We suggest limiting the final concentration of 
NaCl to ≤ 0.4 M after dilution.  Additionally, 
magnesium acetate should be included to keep 
[Mg2+] close to 10 mM."""

    protocol += """\
Apply the diluted reaction mixture to a Amicon 
Ultracel 0.5 ml-100K spin concentrator."""

    protocol += """\
Spin 30 min, 15000g, 4°C."""

    protocol += """\
Add 0.25 volumes of Ni-NTA Agarose to the 
flow-through."""

    protocol += """\
Mix continuously for 30-45 min at 4°C to allow 
His-tagged components to bind the resin."""

    protocol += """\
Apply the reaction mixture slurry to an empty 
Bio-Rad micro-spin column."""

    protocol += """\
Spin 2 min, 1500g, 4°C."""

print(protocol)
