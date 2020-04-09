#!/usr/bin/env python3

import stepwise

p = stepwise.Protocol()

thiol = stepwise.MasterMix.from_text("""\
Reagent                   Stock    Volume
========================  =====  ========
Na₂HPO₄, pH=9.0             1 M  to 35 µL
DTT                         1 M    2.5 µL
puromycin oligo            2 mM     10 µL
""")

emcs = stepwise.MasterMix.from_text("""\
Reagent                    Stock    Volume
========================  ======  ========
water                             to 70 µL
phosphate buffer, pH=7.2     1 M     14 µL
EMCS                      100 mM     10 µL
annealing oligo             1 mM     10 µL
""")

p += f"""\
Reduce the thiol in the puromycin oligo [1,2]:

{thiol}

- Incubate at room temperature for 1h.
- Vortex periodically.
"""

p += """\
Desalt using a NAP-5 column:

- Equilibrate with 9 mL buffer*.
- Load reaction onto the column.
- Load 450 µL buffer* onto the column.
- Elute with 500 µL buffer*.

*20 mM phosphate buffer (pH=7.2)
"""

p.footnotes[1] = """\
Arai2020: 10.1007/978-1-4939-9853-1_3
Mochizuki2011: 10.1021/co2000295
"""

p.footnotes[2] = """\
Use the reduced oligo as soon as possible, 
because the exposed thiols can form disulfide 
bonds with each other.
"""

p += f"""\
Couple EMCS to the amine in the annealing oligo:

{emcs}

- Incubate at 37°C for 30 min [3].
"""

p.footnotes[3] = """\
Do not allow the solution to react for over 30 
min because a longer reaction time causes 
hydrolysis of the maleimide group in EMCS.
"""

p += stepwise.load('ethanol_precipitation')

p += """\
Couple the puromycin and annealing oligos:

- Resuspend the precipitated annealing oligo (now 
  coupled to EMCS) in the eluted puromycin oligo.
- Incubate at 4°C overnight.
- Add 1/20 volume 1M DTT.
- Incubate at room termperature for 30 min.
- Vortex periodically.
"""

p += """\
Perfrom another ethanol precipitation.

- Resuspend in 50 µL water.
"""

p += """\
Purify by HPLC:

- Column: Waters Symmetry 300C18, 4.6x250 mm, 
  particle size 5 μm
- Solvent A: 100 mM TEAA
- Solvent B: acetonitrile/water (80:20, v/v)
- Gradient: B/A (15–36%, 30 min)
- Flow rate: 0.5 mL/min
- Detections: A260, A490
"""

print(p)
