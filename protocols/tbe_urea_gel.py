#!/usr/bin/env python3

"""\
Load, run and stain a TBE/urea gel.

Usage:
    tbe_urea_gel <n>

Arguments:
    <n>
        The number of samples to prepare.
"""

import stepwise
import docopt

args = docopt.docopt(__doc__)

lanes = stepwise.MasterMix.from_text("""\
Reagent            Stock  Volume  MM?
=============  =========  ======  ===
water                       4 µL  yes
sample buffer         2x    5 µL  yes
RNA            200 ng/µL    1 µL
""")
lanes.num_reactions = eval(args['<n>'])

protocol = stepwise.Protocol()

protocol += f"""\
Prepare samples for TBE/urea PAGE:

- Aim for 200 ng/band

- IVT samples:
  - Dilute 10x (2000 ng/µL is a typical yield).
  - Load 1 µL.

- NEB Low Range ssRNA Ladder (N0364)
  - Load 1 µL.

{lanes}

- Incubate at 70°C for 3 min.
"""

protocol += """\
Run the gel:

- 6% TBE/urea PAGE
- Run at 180V for 50 min.
"""

protocol += """\
Stain in 1x PAGE GelRed for 30 min.
"""

print(protocol)

