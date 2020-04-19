#!/usr/bin/env python3
# vim: tw=50

from stepwise import Protocol, MasterMix, load
from attach_huh import huh
from inform import plural

huh['HUH-DNA'].name = "HUH-DNA (f12)"
huh['HUH-DNA'].master_mix = True
huh['sgRNA'].master_mix = False
del huh['EDTA']
huh.num_reactions = 2

p = Protocol()

p += f"""\
Use the following controls:

HUH-DNA (f16):  − + + + + + +
dCas9:          − − + + − + +
sgRNA:          − − − + − − +
BAL-31:         − − − − + + +
                          
Amplification   − + + + − − +
expected?
"""

p += f"""\
Attach HUH-tagged DNA to the dCas9 RNP.

{huh}
"""

bal31 = MasterMix("""\
Reagent            Stock    Volume  MM?
================  ======  ========  ===
water                     to 50 µL  yes
dCas9-PCV2-DNA    100 nM     10 µL
NEBuffer 3.1 [2]     10x      5 µL  yes
BAL-31            1 U/µL      1 µL     
""")
bal31.hold_ratios.volume = '20 µL'
bal31.num_reactions = 3

p += f"""\
Setup {plural(bal31.num_reactions):# BAL-31 digestion/s} [1]:

{bal31}

- Incubate at 30°C for 10 min.
- Add EGTA to 30 mM
- Incubate at 65°C for 10 min.
"""

p.footnotes[1] = """\
The NEB protocol calls for 3 pmol of DNA ends.  
This reaction has 1 pmol DNA, which is 2 pmol 
ends.  Maybe I should use less BAL-31 to 
compensate, but I think it's close enough.
"""

p.footnotes[2] = """\
This is the Cas9 reaction buffer recommended by 
NEB.  I might instead use whatever buffer Jorge 
has been using, or another similar buffer based on 
what I actually have.
"""

# 3x regular volume, for 3 technical replicates 
# per reaction.
p += load("pcr dCas9-PCV2-DNA o87 o88 7 -a 65 -x 10 -p ssoadv -v 66 -m primers")

print(p)
