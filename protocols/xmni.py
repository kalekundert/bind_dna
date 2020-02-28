#!/usr/bin/env python3

"""\
Step up XmnI restriction digest reactions.

Usage:
    xmni [<N>] [options]

Options:
    -p --plasmid NAME
        The name of the plasmid being digested.

    -u --dna-ug AMOUNT
        The amount of DNA to digest, in µg.

    -D --dna-stock CONC
        The concentration of DNA.
"""

from docopt import docopt
from stepwise import Protocol, MasterMix
from nonstdlib import plural

args = docopt(__doc__)
xmni = MasterMix.from_text("""\
Reagent              Stock    Volume  MM?
===============  =========  ========  ===
water                       to 50 µL  yes
DNA              200 ng/µL      5 µL
CutSmart buffer        10x      5 µL  yes
XmnI               20 U/µL      1 µL  yes
""")

xmni.num_reactions = eval(args['<N>'] or '1')
xmni['XmnI'].name = 'XmnI [1]'

if x:= args['--plasmid']:
    xmni['DNA'].name = x

if x := args['--dna-stock']:
    xmni['DNA'].stock_conc = eval(x)

if x := args['--dna-ug']:
    xmni['DNA'].volume *= eval(x) / 1
    xmni['XmnI'].volume *= eval(x) / 1
    xmni.volume = max(
            xmni.volume,
            10 * xmni['XmnI'].volume,
    )

protocol = Protocol()

protocol += f"""\
Setup {plural(xmni.num_reactions):? XmnI digestion/s}:

{xmni}
"""

protocol += """\
Incubate at 37°C for 5–15 min.
"""

protocol.footnotes[1] = """\
NEB recommends 5–10 units of enzyme per µg DNA 
(10–20 units for genomic DNA) in a 1 hour digest. 
Enzyme volume should not exceed 10% of the total 
reaction volume to prevent star activity due to 
excess glycerol.
"""

if __name__ == '__main__':

    print(protocol)
