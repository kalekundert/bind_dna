#!/usr/bin/env python3

"""\
Setup in vitro transcription/translation (IVTT) reactions using the NEB 
PURExpress system (E6800).

Usage:
    purexpress.py <num_rxns> [-v <uL>] [-D <nM>] [options]

Options:
    -v --rxn-volume <uL>                        [default: 10]
        The volume of each individual reaction (in µL).  NEB recommends 25 µL, 
        but I typically use 10 µL and get enough yield for routine experiments.

    -D --dna-stock-conc <nM>                    [default: 75]
        The concentration of the DNA being added to the reaction.  If this 
        differs from the default, the volume of DNA to add will be adjusted 
        accordingly.

    -r --mrna
        Use mRNA as the template instead of DNA.

    -g --add-target
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

    -t --time <time>                            [default: 2h]
        The amount of time to run the transcription/translation reaction for.  
        No unit is assumed, so be sure to specify one.
"""

import docopt
import stepwise
from inform import plural

args = docopt.docopt(__doc__)
protocol = stepwise.Protocol()
purexpress = stepwise.MasterMix.from_text('''\
Reagent              Stock      Volume  MM?
===============  =========  ==========  ===
water                       to 10.0 µL  yes
A                               4.0 µL  yes
B                               3.0 µL  yes
RNase Inhibitor    40 U/µL      0.2 µL  yes
ZnOAc                 1 mM      0.5 µL  yes
target DNA          750 nM      0.8 µL  yes
template DNA         75 nM      0.8 µL
template mRNA        10 µM      1.6 µL
''')

purexpress['water']

if not args['--add-zinc']:
    del purexpress['ZnOAc']

if not args['--add-target']:
    del purexpress['target DNA']

if args['--mrna']:
    del purexpress['template DNA']
else:
    del purexpress['template mRNA']

purexpress.num_reactions = eval(args['<num_rxns>'])
purexpress.hold_ratios.volume = eval(args['--rxn-volume']), 'µL'

protocol += f"""\
Setup {plural(purexpress.num_reactions):# PURExpress reaction/s}:

{purexpress}

- Keep on ice.
- Be sure to add A before B.
{'- The control template (125 ng/µL) is 75 nM.' if not args['--mrna'] else ''}
"""

protocol += f"""\
Incubate at 37°C for {args['--time']}."""

if args['--purify']:
    protocol += """\
Dilute reaction to 100 µL with PBS + 10 mM 
MgOAc[1].
- Save a 10 µL aliquot"""

    protocol.footnotes[1] = """\
We recommend a minimum volume of 100 µl after 
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
- Save a 10 µL aliquot of the flow-through.
- Dilute the retentate to 100 µL, then save a
  10 µL aliquot."""

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
- Save a 10 µL aliquot of the eluate."""

print(protocol)
