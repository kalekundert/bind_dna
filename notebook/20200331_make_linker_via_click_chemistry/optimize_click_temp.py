#!/usr/bin/env python3

import stepwise

p = stepwise.Protocol()
rxn = stepwise.MasterMix("""\
        Reagent   Stock  Volume  MM?
        =======  ======  ======  ===
        o126     400 µM    1 µL  yes
        o127     400 µM    1 µL  yes
        PBS          2x    2 µL  yes
""")
rxn.num_reactions = 4

p += f"""\
Setup 4 click reactions:

{rxn}
"""

p += """\
Incubate the reactions as follows:

- 25°C for 4h
- 25°C for 3h30, 65°C for 30 min
- 37°C for 4h
- 37°C for 3h30, 65°C for 30 min
"""

p += stepwise.load("gel anneal 6 -c 1990")

p.steps[2] = p.steps[2].replace(':', ' [1]:')
p.footnotes[1] = """\
The product (o129) has MW = 19899 g/mol, so:
100 µM = 1990 ng/µL
"""

print(p)

