#!/usr/bin/env python3

import stepwise
import wang2017

purex = stepwise.MasterMix("""\
Reagent                            Stock      Volume  MM?
===============================  =======  ==========  ===
water                                     to 10.5 µL   x
A                                             4.0 µL   x
B                                             3.0 µL   x
RNase inhibitor (murine)         40 U/µL      0.2 µL   x
o210,o211                           5 µM      0.5 µL
""")
nebex = stepwise.MasterMix("""\
Reagent                            Stock      Volume  MM?
=============================  =========  ==========  ===
water                                     to 10.5 µL   x
S30 extract                                   2.4 µL   x
synthesis buffer                      2x      5.0 µL   x
T7 RNA polymerase               450 U/µL      0.2 µL   x
RNase inhibitor (murine)         40 U/µL      0.2 µL   x
GamS nuclease inhibitor        1.5 µg/µL      0.2 µL   x
o210,o211                           5 µM      0.5 µL
""")
s30 = stepwise.MasterMix("""\
Reagent                            Stock      Volume  MM?
==============================  ========  ==========  ===
water                                     to 10.5 µL   x
A                                    10x      1.0 µL   x
B                                   2.5x      4.0 µL   x
C                                  3.33x      3.0 µL   x
DNA                                75 nM      0.8 µL   x
o210,o211                           5 µM      0.5 µL
""")

nebex.num_reactions = 2
purex.num_reactions = 2
s30.num_reactions = 2

p = stepwise.Protocol()

p += """\
Setup 3 controls (−,−,+):

- 10.0 µL water
-  0.5 µL 5 µM o210,o211,o212
"""

p += f"""\
Setup an RNase H assay for PURExpress:

{purex}
"""

p += f"""\
Setup an RNase H assay for NEBExpress:

{nebex}
"""

p += f"""\
Setup an RNase H assay for Promega S30 extract:

{s30}
"""

p += wang2017.IncubateAndMeasure(9)

p.print()




