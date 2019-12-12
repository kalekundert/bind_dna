#!/usr/bin/env python3

"""\
Setup in vitro transcription/translation (IVTT) reactions using the NEB 
PURExpress system (E6800).

Usage:
    purexpress.py <num_rxns> [-v <uL>] [-D <nM>] [-t] [-z] [-p] [-n] [-s]

Options:
    -v --rxn-volume <uL>  [default: 10]
        The volume of each individual reaction (in μL).  NEB recommends 25 μL, 
        but I typically use 10 μL and get enough yield for routine experiments.

    -D --dna-stock-conc <nM>  [default: 75]
        The concentration of the DNA being added to the reaction.  If this 
        differs from the default, the volume of DNA to add will be adjusted 
        accordingly.

    -t --add-target
        Add target DNA to the PURExpress reaction.

    -z --add-zinc
        Add ZnOAc to the PURExpress reaction.  This is necessary when 
        expressing Zn-finger proteins.

    -p --purify
        Purify the expressed protein using the reverse-His protocol recommended 
        by NEB.

    -n --native-page
        Run the IVTT reaction on a native PAGE gel.

    -s --sds-page
        Run the IVTT reaction on an SDS-PAGE gel.
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
purexpress['template DNA'].stock_conc = args['--dna-stock-conc']

purexpress.num_reactions = eval(args['<num_rxns>'])
purexpress.volume = eval(args['--rxn-volume'])

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
Dilute reaction to 100 μL with PBS + 10 mM 
MgOAc[1].
- Save a 10 μL aliquot"""

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
Spin 30 min, 15000g, 4°C.
- Save a 10 μL aliquot of the flow-through.
- Dilute the retentate to 100 μL, then save a
  10 μL aliquot."""

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
Spin 2 min, 1500g, 4°C.
- Save a 10 μL aliquot of the eluate."""

if args['--sds-page']:
    protocol += """\
Setup an SDS-PAGE gel:

- Prepare samples:
  - 10.00 μL IVTT reaction
  -  3.85 μL 4x loading buffer
  -  1.54 μL 10x reducing agent
  - Incubate at 70°C for 10 min.

- Use a 4-12% gel (Invitrogen NW04120).
- Load 15.39 μL in each lane.
- Run at 165V for 42 min."""

if args['--native-page']:
    protocol += """\
Setup a native PAGE gel:

- DNA ladder:
  - 2.5 μL water
  - 5.0 μL 50 ng/μL ladder, i.e. 1kb+ (NEB N3232)
  - 2.5 μL 4x sample buffer (Invitrogen BN2003)

- IVTT reactions:
  - 6.25 μL water
  - 1.25 μL IVTT reaction
  - 2.50 μL 4x sample buffer (Invitrogen BN2003)

- Use a 3-12% gel (Invitrogen BN1003).
- Load 5 μL in each lane.
- Run at 150V for 115 min."""


print(protocol)
