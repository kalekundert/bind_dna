#!/usr/bin/env python3

import stepwise
from stepwise import Protocol, Step, MasterMix

p = Protocol()

loading_buf = MasterMix("""\
Reagent                 Stock     Volume MM?
======================  =====  ========= ===
water                          to 5.0 µL   y
Bolt LDS sample buffer     4x     2.5 µL   y
Bolt reducing agent       10x     1.0 µL   y
""")
loading_buf.num_reactions = 8
loading_buf.extra_percent = 50

protein_samples = MasterMix("""\
Reagent                        Stock      Volume MM?
======================     =========  ========== ===
water                                 to 10.0 µL   y
loading buffer master mix         2x      5.0 µL   y
BSA                        500 ng/µL      1.0 µL   y
""")
protein_samples.num_reactions = 4
protein_samples.extra_percent = 50

dna_samples = MasterMix("""\
Reagent                        Stock      Volume MM?
======================     =========  ========== ===
loading buffer master mix         2x      5.0 µL   y
1 kb Plus DNA ladder        50 ng/µL      5.0 µL   y
""")
dna_samples.num_reactions = 4
dna_samples.extra_percent = 50

p += """\
Prepare 500 ng/µL BSA:
- 199 µL water
- 1 µL 0.1% BSA
"""

p += Step("Prepare known DNA and protein samples")
p.s.body += loading_buf
p.s.body += protein_samples
p.s.body += dna_samples

p += Step("""\
Load the gel such that it can be cut into four pieces, each with a protein lane 
and a DNA lane.
""")

p += """\
Run a gel:

- Use a 4−12% SDS PAGE gel.
- Load 10 µL of each sample.
- Run at 165V for 42 min.
"""

p += Step("Cut the gel as described above.")

p += """\
Stain the gel pieces with the following protocols:

A: Coomassie only

- Repeat 3 times:
    - Add 100 mL water.
    - Microwave until almost boiling (1 min).
    - Shake gently for 1 min.

- Submerge gel in SimplyBlue SafeStain.
- Microwave until almost boiling (45-60s).
- Shake gently for 5 min.
- Wash the gel with 100 mL water for 10 min.
- Add 20 mL 5M NaCl and wash for 5 min.

B: GelGreen only

- Submerge gel in 3x GelGreen, 100 mM NaCl.
- Shake gently for 30 min.

C: Coomassie then GelGreen

- Repeat 3 times:
  - Add 100 mL water.
  - Microwave until almost boiling (1 min).
  - Shake gently for 1 min.

- Submerge gel in SimplyBlue SafeStain.
- Microwave until almost boiling (45 sec).
- Shake gently for 5 min.

- Discard stain and add 100 mL water.
- Shake gently for 10 min.

- Submerge gel in 3x GelGreen, 100 mM NaCl.
- Shake gently for 30 min.

D: GelGreen then Coomassie

- Repeat 3 times:
  - Add 100 mL water.
  - Microwave until almost boiling (1 min).
  - Shake gently for 1 min.

- Submerge gel in 3x GelGreen, 100 mM NaCl.
- Shake gently for 30 min.

- Submerge gel in SimplyBlue SafeStain.
- Microwave until almost boiling (45 sec).
- Shake gently for 5 min.

- Discard stain and add 100 mL water.
- Shake gently for 10 min.
- Add 20 mL 5M NaCl.
- Shake gently for 5 min.
"""

p += stepwise.load('laser blue red')

p.print()

