#!/usr/bin/env python3

import stepwise
import wang2017 

p = stepwise.Protocol()

p += stepwise.load('serial 12µL 10U/mL 0.01 7 -m "RNase H" -d "RNase H reaction buffer"')

p += """\
Setup 1 positive control:

- 10.0 µL water
-  0.5 µL 5 µM o212
"""

p += """\
Setup 8 negative controls:

- 10.0 µL RNase H dilution, including 0 U/mL.
-  0.5 µL 5 µM o211 [1]
"""

p += """\
Setup 8 experimental samples:

- 10.0 µL RNase H dilution, including 0 U/mL.
-  0.5 µL 5 µM o210 [1]
"""

p += wang2017.IncubateAndMeasure(17)

p.footnotes[1] = """\
If I could use an Echo, I would add 25 nL of 100 µM 
o210/o211, to minimize the amount of volume added to 
the sample.  But I don't think I can accurately 
pipet less than 0.5 µL by hand.
"""

p.print()
