#!/usr/bin/env python3

"""\
Setup in vitro transcription/translation (IVTT) reactions using the NEB 
PURExpress system (E6800).

Usage:
    purexpress.py <templates>... [-v <uL>] [-D <nM>] [options]

Arguments:
    <templates>
        The names of the DNA/mRNA templates to express.

Options:
    -v --rxn-volume <uL>                        [default: 10]
        The volume of each individual reaction (in µL).  NEB recommends 25 µL, 
        but I typically use 10 µL and get enough yield for routine experiments.

    -t --template-conc <nM>
        The desired final concentration of template in the reaction.  The 
        default differs depending on whether the template is DNA or mRNA.  If 
        this differs from the default, the volume of the template will be 
        adjusted accordingly.

    -T --template-stock <nM>
        The concentration of the DNA/RNA being added to the reaction.  The 
        default depends on whether the template is DNA or mRNA.  If this 
        differs from the default, the volume of the template will be adjusted 
        accordingly.

    -w --time <time>                            [default: 2h]
        The amount of time to run the transcription/translation reaction for.  
        No unit is assumed, so be sure to specify one.

    -r --mrna
        Use mRNA as the template instead of DNA.

    -p --purify
        Purify the expressed protein using the reverse-His protocol recommended 
        by NEB.

    -g --add-target
        Add target DNA to the PURExpress reaction.

    -z --add-zinc
        Add ZnOAc to the PURExpress reaction.  This is necessary when 
        expressing Zn-finger proteins.

"""

import docopt
import stepwise
from inform import plural

args = docopt.docopt(__doc__)
protocol = stepwise.Protocol()
purexpress = stepwise.MasterMix.from_text('''\
Reagent                  Stock      Volume  MM?
===================  =========  ==========  ===
water                           to 10.0 µL  yes
A                                   4.0 µL  yes
B                                   3.0 µL  yes
RNase Inhibitor [1]    40 U/µL      0.2 µL  yes
ZnOAc                     1 mM      0.5 µL  yes
target DNA              750 nM      0.8 µL  yes
template DNA             75 nM      0.8 µL
template mRNA          1000 nM      1.6 µL
''')

if not args['--add-zinc']:
    del purexpress['ZnOAc']

if not args['--add-target']:
    del purexpress['target DNA']

if args['--mrna']:
    template = 'template mRNA'
    default_template_conc_nM = 160
    default_template_stock_nM = 1000
    del purexpress['template DNA']
else:
    template = 'template DNA'
    default_template_conc_nM = 6
    default_template_stock_nM = 75
    del purexpress['template mRNA']

def float_or_default(x, default):
    return float(x) if x else default

purexpress[template].name = ','.join(args['<templates>'])
purexpress[template].hold_stock_conc.conc = float_or_default(
        args['--template-conc'], default_template_conc_nM), 'nM'
purexpress[template].hold_conc.stock_conc = float_or_default(
        args['--template-stock'], default_template_stock_nM), 'nM'

purexpress.num_reactions = len(args['<templates>'])
purexpress.hold_ratios.volume = eval(args['--rxn-volume']), 'µL'

purexpress.fix_volumes(template)

protocol += f"""\
Setup {plural(purexpress.num_reactions):# PURExpress reaction/s}:

{purexpress}

- Keep on ice.
- Be sure to add A before B.
{'- The control template (125 ng/µL) is 75 nM.' if not args['--mrna'] else ''}
"""

protocol.footnotes[1] = """\
The PURExpress protocol recommends 0.8 U/µL (20 U 
per 25 µL reaction), while the product page for 
the inhibitor itself recommends 1 U/µL.  I'm using 
the former here because it's the recommendation I 
encountered first.

PURExpress protocol: https://tinyurl.com/y3m9lrcz
RNAse inhibitor FAQs: https://tinyurl.com/y3zabsoz
"""

protocol += f"""\
Incubate at 37°C for {args['--time']}."""

if args['--purify']:
    protocol += """\
Dilute reaction to 100 µL with PBS + 10 mM 
MgOAc [2].
- Save a 10 µL aliquot"""

    protocol.footnotes[2] = """\
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
