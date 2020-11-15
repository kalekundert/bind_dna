#!/usr/bin/env python3

import stepwise
from inform import plural
from stepwise import Protocol, Step, Footnote, MasterMix

p = Protocol()
n = 6

loading_buf = MasterMix("""\
Reagent                 Stock     Volume MM?
======================  =====  ========= ===
water                          to 5.0 µL   y
Bolt LDS sample buffer     4x     2.5 µL   y
Bolt reducing agent       10x     1.0 µL   y
""")
loading_buf.num_reactions = 2 * n
loading_buf.extra_percent = 50

protein_samples = MasterMix("""\
Reagent                        Stock      Volume MM?
=========================  =========  ========== ===
water                                 to 10.0 µL   y
loading buffer master mix         2x      5.0 µL   y
BSA                        500 ng/µL      1.0 µL   y
""")
protein_samples.num_reactions = n
protein_samples.extra_percent = 10

dna_samples = MasterMix("""\
Reagent                        Stock      Volume MM?
=========================  =========  ========== ===
loading buffer master mix         2x      5.0 µL   y
1 kb Plus DNA ladder        50 ng/µL      5.0 µL   y
""")
dna_samples.num_reactions = n
dna_samples.extra_percent = 10

p += """\
Prepare 500 ng/µL BSA:
- 199 µL water
- 1 µL 0.1% BSA
"""

p += Step("Prepare known DNA and protein samples:")
p.s.body += loading_buf
p.s.body += protein_samples
p.s.body += dna_samples
p.s += "Incubate at 70°C for 10 min"

p += Step(f"""\
Load the gel such that it can be cut into {plural(n):# piece/s}, each with a 
protein lane and a DNA lane.
""")

p += """\
Run a gel:

- Use a 4−12% SDS PAGE gel.
- Load 10 µL of each sample.
- Run at 165V for 42 min.
"""

p += """\
Cut the gel as described above.  Identify each 
piece by notching the corners according to a 
binary code:

- Top left:     8 if notched
- Top right:    4 if notched
- Bottom left:  2 if notched
- Bottom right: 1 if notched
"""

p += """\
Prepare 90 mL 3x GelGreen, 100 mM NaCl.
"""

p += """\
Stain the gel pieces with the following protocols:

0: Coomassie only
1: GelGreen only (− wash)
2: GelGreen only (+ wash)
3: Coomassie then GelGreen
4: GelGreen then Coomassie (− wash)
5: GelGreen then Coomassie (+ wash)

0,2,3,5: Wash the gel

- Repeat 3 times:
  - Add 100 mL water.
  - Microwave until almost boiling (1 min).
  - Shake gently for 1 min.

1,4: Stain with GelGreen (− wash) [1]

- Rinse gel with water.
- Submerge gel in 3x GelGreen, 100 mM NaCl.
- Shake gently for 30 min.

2,5: Stain with GelGreen (+ wash) [1]

- Submerge gel in 3x GelGreen, 100 mM NaCl.
- Shake gently for 30 min.

0,3: Stain with Coomassie:

- Submerge gel in SimplyBlue SafeStain.
- Microwave until almost boiling (45-60s).
- Shake gently for 5 min.
- Discard stain and add 100 mL water.
- Shake gently for 10 min.

0: Wash with salt:

- Add 20 mL 5M NaCl.
- Shake gently for 5 min.

3: Stain with GelGreen:

- Submerge gel in 3x GelGreen, 100 mM NaCl.
- Shake gently for 30 min.

4,5: Stain with Coomassie

- Submerge gel in SimplyBlue SafeStain.
- Microwave until almost boiling (45-60s).
- Shake gently for 5 min.
- Discard stain and add 100 mL water.
- Shake gently for 10 min.
- Add 20 mL 5M NaCl.
- Shake gently for 5 min.
"""

p.footnotes[1] = Footnote("""\
Stain the washed and unwashed gel pieces 
separately, because it's possible that the excess 
SDS from the unwashed gels could interfere with 
staining.
""")

p += stepwise.load('laser blue red')

p.print()

