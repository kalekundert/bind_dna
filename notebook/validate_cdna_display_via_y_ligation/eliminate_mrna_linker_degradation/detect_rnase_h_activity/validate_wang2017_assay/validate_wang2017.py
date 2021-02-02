#!/usr/bin/env python3

import stepwise
import wang2017 

p = stepwise.Protocol()

p += """\
Prepare 100 U/mL RNase H:

- 44 µL nuclease-free water
- 5 µL 10x RNase H buffer
- 1 µL 5000 U/mL RNase H
"""

p += stepwise.load('serial 12µL 100U/mL 0.1 7 -m "RNase H" -d "1x RNase H reaction buffer"')

p.footnotes[2] = """\
If I could use an Echo, I would add 25 nL of 100 µM 
o210/o211, to minimize the amount of volume added to 
the sample.  But I don't think I can accurately 
pipet less than 0.5 µL by hand.
"""

p += """\
Prepare 5 µM o210-o212:

- 1 µL 100 µM oligo
- 19 µL nuclease-free water
"""

p += """\
Setup 1 blank:

- 10.0 µL 1x RNase H reaction buffer
-  0.5 µL nuclease-free water [2]
"""

p += """\
Setup 1 positive control:

- 10.0 µL 1x RNase H reaction buffer
-  0.5 µL 5 µM o212 [2]
"""

p += """\
Setup 1 negative control:

- 10.0 µL 100 U/mL RNase H
-  0.5 µL 5 µM o211 [2]
"""

p += """\
Setup 8 experimental samples:

- 10.0 µL RNase H dilution, including 0 U/mL.
-  0.5 µL 5 µM o210 [2]
"""

p += wang2017.IncubateAndMeasure(
        num_reactions=11,
        rnase_incubation_time_min=120,
        rnase_incubation_temp_C=37,
)

p.print()
