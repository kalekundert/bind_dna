#!/usr/bin/env python3
# vim: tw=50

import stepwise

huh = stepwise.MasterMix.from_text("""\
Reagent    Stock    Volume  MM?
========  ======  ========  ===
water             to 10 µL  yes
dCas9       1 µM      1 µL  yes
sgRNA       1 µM      1 µL  yes
HUH-DNA   200 nM      5 µL
EDTA      500 mM      1 µL
CutSmart     10x      1 µL  yes
""")

if __name__ == '__main__':
    p = stepwise.Protocol()

    p += """\
Reactions:

dCas9:          +  +  +  +
DNA (f12,f16):  −  +  +  +
HUH-target:     −  −  +  +
EDTA:           −  −  −  +
"""

    huh['dCas9'].name = "dCas9 [1]"
    huh['sgRNA'].name = "sgRNA (f13)"
    huh['HUH-DNA'].name = "HUH-DNA (Ø,f16,f12) [2]"
    huh['EDTA'].name = "EDTA [3]"
    huh.num_reactions = 4

    p += f"""\
Attach HUH-tagged DNA to the Cas9 RNP:

{huh}

- Incubate at 37°C for 15 min.
"""

    p.footnotes[1] = """\
Invitrogen recommends loading no more than 250 
ng/band on Bolt SDS PAGE gels, and the detection 
limit for Coomassie (as an IR dye [Butt2013]) is 
at least 10 ng/band.  That corresponds to a range 
of ≈1.4-0.1 pmol Cas9-PCV2/band (MW: 176.8 kDa).  
I probably want to be on the high end of that.
"""

    p.footnotes[2] = """\
f16 and f12 are 414 bp.  At that length, 50 ng/µL 
(a typical PCR yield) corresponds to ≈150 nM.  
[Lovendahl2017] used a 10:1 DNA:protein ratio to
maximize the amount of coupled protein.  I'm going 
to use a 1:1 ratio instead, both because I don't 
want a lot of unbound DNA in my qPCR reactions and 
because a 10:1 ratio would use a lot of material.
"""

    p.footnotes[3] = """\
[VegaRocha2007] used 2.5 mM divalent metal and
30 mM EDTA.  This reaction has 10 mM Mg²⁺ and 50 
mM EDTA.  Not the same proportion, but still an 
excess.
"""

    p += stepwise.load('gel sdsmax 3 -S')
    p += stepwise.load('gelgreen -r')
    p += stepwise.load('simplyblue -r')

    print(p)
