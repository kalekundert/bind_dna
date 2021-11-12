#!/usr/bin/env python3

import docopt
import stepwise
import wang2017

controls = stepwise.MasterMix("""\
Reagent                            Stock      Volume  MM?
===============================  =======  ==========  ===
nuclease-free water                       to 10.5 µL   x
o210,o211,o212                      5 µM      0.5 µL
""")

purex = stepwise.MasterMix("""\
Reagent                            Stock      Volume  MM?
===============================  =======  ==========  ===
nuclease-free water                       to 10.5 µL   x
A                                             4.0 µL   x
B                                             3.0 µL   x
RNase inhibitor (murine)         40 U/µL      0.2 µL   x
o210,o211,o212                      5 µM      0.5 µL
""")

pfrex = stepwise.MasterMix("""\
Reagent                            Stock      Volume  MM?
===============================  =======  ==========  ===
nuclease-free water                       to 10.5 µL   x
Solution I                                    5.0 µL   x
Solution II                                   0.5 µL   x
Solution III                                  1.0 µL   x
o210,o211,o212                      5 µM      0.5 µL
""")

nebex = stepwise.MasterMix("""\
Reagent                            Stock      Volume  MM?
=============================  =========  ==========  ===
nuclease-free water                       to 10.5 µL   x
A: S30 extract                                2.4 µL   x
B: synthesis buffer                   2x      5.0 µL   x
RNase inhibitor (murine)         40 U/µL      0.2 µL   x
GamS nuclease inhibitor        1.5 µg/µL      0.2 µL   x
o210,o211,o212                      5 µM      0.5 µL
""")

s30 = stepwise.MasterMix("""\
Reagent                            Stock      Volume  MM?
==============================  ========  ==========  ===
nuclease-free water                       to 10.5 µL   x
A                                    10x      1.0 µL   x
B                                   2.5x      4.0 µL   x
C                                  3.33x      3.0 µL   x
o210,o211,o212                      5 µM      0.5 µL
""")

controls.num_reactions = 3
purex.num_reactions = 3
pfrex.num_reactions = 3
nebex.num_reactions = 3
s30.num_reactions = 3

p = stepwise.Protocol()

p += stepwise.load('zap')

p += """\
Set the plate reader temperature to 30°C.
"""

p += """\
Prepare 5 µM o210-o212:

- 19 µL nuclease-free water
-  1 µL 100 µM oligo
"""

p += """\
Setup two blank wells (one will become a FAM
positive control):

- 10.5 µL nuclease-free water
"""

p += f"""\
Setup control reactions:

{controls}
"""

p += f"""\
Setup an RNase H assay for PURExpress:

{purex}
"""

p += f"""\
Setup an RNase H assay for PUREfrex:

{pfrex}
"""

#p += f"""\
#Setup an RNase H assay for NEBExpress:
#
#{nebex}
#"""

#p += f"""\
#Setup an RNase H assay for Promega S30 extract:
#
#{s30}
#"""

# PURExpress: 37°C for 2 h
# PUREfrex: 37°C for 2-4h
# NEBExpress: 37°C for 2-4 h
# Promega S30 Extract (for Linear Templates): 37°C for 1-2 h

# Barendt2013: 37°C for 30 min.

# expt 65: 37°C for 2h
# expt 61: 37°C for 30 min, same result as expt 65 (all mRNA/DNA degraded).

p += wang2017.IncubateAndMeasure(
        num_reactions=9,
        rnase_incubation_time_min=30,
        rnase_incubation_temp_C=37,
        fam_control=True,
)

p.print()




