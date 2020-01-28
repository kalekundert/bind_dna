#!/usr/bin/env python3

import stepwise

anneal = stepwise.MasterMix.from_text("""\
Reagent       Stock   Volume  MM?
=======  ==========  =======  ===
water                 2.6 µL  yes
PBS             10x   0.4 µL  yes
o93      10 pmol/µL   0.5 µL  yes
f11      10 pmol/µL   0.5 µL  yes
""")
anneal.num_reactions = 7

ligate = stepwise.MasterMix.from_text("""\
Reagent                     Stock   Volume  MM?
====================  ===========  =======  ===
water                              2.75 µL  yes
T4 DNA ligase buffer          10x     1 µL  yes
BSA                          0.1%     1 µL  yes
T4 RNA ligase             40 U/µL  1.25 µL  yes
Annealed oligos       2.5 pmol/µL     4 µL  yes
""")
ligate.num_reactions = 4

protocol = stepwise.Protocol()
protocol += f"""\
Setup the annealing reactions:

{anneal}

Also setup individual reactions lacking:
- PBS
- o93
- f11
"""

protocol += """\
Run the following thermocyler protocol:

- 95°C for 3 min
- 95°C to 25°C in 700 2s steps of 0.1°C.

Save a "−anneal" control.
"""

protocol += """\
Dialyze the reaction:

- Resevoir of RNase-free water.
- 0.025 µm drop dialysis membrane.
- Incubate 1h at 25°C.

Save a "−dialysis" control.
"""

protocol += f"""\
Setup the ligation reactions:

{ligate}

Also setup individual reactions lacking:
- T4 ligase buffer
- BSA
- T4 RNA ligase

Save a "−ligation" control.
"""

print(protocol)
