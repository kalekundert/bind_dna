#!/usr/bin/env python3
# vim: tw=50

"""\
Usage:
    attach_huh.py [-S]

Arguments:
    -S --no-sgrna
        Leave the sgRNA out of the reaction.
"""

import docopt
import stepwise
import po4
from statistics import mean

huh = stepwise.MasterMix.from_text("""\
Reagent    Stock    Volume  MM?
========  ======  ========  ===
water             to 10 µL  yes
CutSmart     10x      1 µL  yes
dCas9       1 µM      1 µL  yes
sgRNA       1 µM      1 µL  yes
EDTA      500 mM      1 µL
HUH-DNA   200 nM      5 µL
""")

if __name__ == '__main__':
    args = docopt.docopt(__doc__)
    p = stepwise.Protocol()

    p += stepwise.load('zap')

    p += """\
In the following steps, setup these reactions:

dCas9:          −  +  +  +  +
DNA (f12,f16):  +  −  +  +  +
HUH-target:     +  −  −  +  +
EDTA:           −  −  −  −  +
"""
    huh.num_reactions = 4
    huh.extra_min_volume = '0.5 µL'

    huh['water'].name = "nuclease-free water"
    huh['dCas9'].name = "dCas9-PCV2 [1]"
    huh['sgRNA'].name = "sgRNA (f13)"
    huh['EDTA'].name = "EDTA [2]"
    huh['HUH-DNA'].name = "HUH-DNA (Ø,f16,f12) [3]"

    huh['dCas9'].hold_conc.stock_conc = 22, 'µM'
    huh['sgRNA'].hold_conc.stock_conc = 10, 'µM'

    db = po4.load_db()
    # Get an error when trying to get the MW for 
    # f12, because PO₄ doesn't know how to account 
    # for the nonstandard nucleotide correctly.  
    # So just use f16 and assume that f12 is the 
    f16_nM = round(db['f16'].conc_nM)
    huh['HUH-DNA'].hold_conc.stock_conc = f16_nM, 'nM'

    if args['--no-sgrna']:
        del huh['sgRNA']

    p += f"""\
Attach HUH-tagged DNA to the Cas9 RNP:

{huh}

- Add each reagent in order.
- Mix after adding the EDTA (and before adding the 
  HUH-DNA) [2].
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
The EDTA reaction is a negative control 
established in [VegaRocha2007].  They used 2.5 mM 
divalent metal and 30 mM EDTA to prevent coupling.  
This reaction has 10 mM Mg²⁺ and 50 mM EDTA.  Not 
the same proportion, but still an excess.
"""

    p.footnotes[3] = """\
f16 and f12 are 414 bp.  At that length, 50 ng/µL 
(a typical PCR yield) corresponds to ≈150 nM.  
[Lovendahl2017] used a 10:1 DNA:protein ratio to
maximize the amount of coupled protein.  I'm going 
to use a 1:1 ratio instead, both because I don't 
want a lot of unbound DNA in my qPCR reactions and 
because a 10:1 ratio would use a lot of material.
"""

    p += stepwise.load('gel sdsmax 5 -S')
    p += stepwise.load('gelgreen -r')
    p += stepwise.load('simplyblue -fr')

    print(p)
