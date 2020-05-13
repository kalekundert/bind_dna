#!/usr/bin/env python3
# vim: tw=50

"""\
Test small numbers of B1H constructs by plating on selective plates.

Usage:
    plate_assay <strains>...

Arguments:
    <strains>
        The strains to individually test.
"""

import docopt
import stepwise

args = docopt.docopt(__doc__)
p = stepwise.Protocol()

additives = [
        ['His',  '3-AT'],
        ['5 mM', '0 mM'],
        ['0 mM', '1 mM'],
        ['0 mM', '2 mM'],
        ['0 mM', '4 mM'],
]

p += f"""\
Make NM plates with the following additives:

{stepwise.tabulate(additives, True)}
"""

p += f"""\
Pick fresh colonies of the following strains and 
grow to OD≈1 in 3 mL LB+Carb+Kan: {','.join(args['<strains>'])}
"""

p += """\
Wash cells with minimal media:

- Spin 3500g, 3 min, 4°C.
- Resuspend in 1 mL NM.
- Spin 3500g, 3 min, 4°C.
- Resuspend in 1 mL NM.
- Dilute to OD=0.1
"""

p += stepwise.load('serial 90µL 0.1OD 6 -f 10 -m cells -d NM')

p += """\
Plate 5 µL of each dilution in triplicate on each 
plate.
"""

print(p)
