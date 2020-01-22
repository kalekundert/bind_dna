#!/usr/bin/env python3

"""\
Express protein using Promega S30 Extract for linear templates (L1030)

Usage:
    s30 <num_rxns>

Arguments:
    <num_reactions>
        The number of reactions to perform.
"""

import docopt
from math import ceil
from stepwise import Protocol, MasterMix
from inform import plural

args = docopt.docopt(__doc__)
n = float(args['<num_rxns>'])
n_aliquots = int(ceil(n/2))

s30 = MasterMix.from_text("""\
Reagent  Stock  Volume  MM?
=======  =====  ======  ===
water           1.2 µL  yes
A          10x    1 µL  yes
B         2.5x    4 µL  yes
C        3.33x    3 µL  yes
DNA      75 nM  0.8 µL
""")

protocol = Protocol()

protocol += f"""\
Thaw {plural(n_aliquots):# S30 extract aliquot/s} [1].
"""

protocol += f"""\
Setup the S30 extract reaction:

{s30}

- Vortex gently and spin to collect.
"""

protocol += f"""\
Incubate at 37°C for 1h [2].
"""

protocol += f"""\
Stop reaction by incubating at 4°C.
"""

protocol.footnotes[1] = """\
A: Complete amino acids
B: S30 premix without amino acids
C: S30 extract, linear
"""

protocol.footnotes[2] = """\
Enhanced expression at lower temperatures for 
longer times appears to be gene/protein-specific 
and may be tested if the standard reaction at 37°C 
for 1 hour does not produce the desired results.
"""

print(protocol)

# vim: tw=50
