#!/usr/bin/env python3

import stepwise
from stepwise import pl, ul

p = stepwise.Protocol()
p.footnotes[1] = 'https://tinyurl.com/5b9sujn2'

lb = stepwise.MasterMix("""\
  Reagent               Stock   Final     Volume
  ================  =========  ======  =========
  water                                to 200 µL
  tris HCl, pH=7.5    1000 mM  100 mM
  LiCl                7500 mM  500 mM
  LDS                     20%    0.5%
  EDTA                 500 mM    1 mM
  DTT                  500 mM    5 mM
""")
w1 = stepwise.MasterMix("""\
  Reagent               Stock   Final     Volume
  ================    =======  ======  =========
  water                                to 800 µL
  tris HCl, pH=7.5    1000 mM   20 mM
  LiCl                7500 mM  500 mM
  LDS                     20%    0.1%
  EDTA                 500 mM    1 mM
  DTT                  500 mM    5 mM
""")
w2 = stepwise.MasterMix("""\
  Reagent               Stock   Final     Volume
  ================    =======  ======  =========
  water                                to 800 µL
  tris HCl, pH=7.5    1000 mM   20 mM
  LiCl                7500 mM  500 mM
  EDTA                 500 mM    1 mM
""")
ls = stepwise.MasterMix("""\
  Reagent               Stock   Final     Volume
  ================    =======  ======  =========
  water                                to 400 µL
  tris HCl, pH=7.5    1000 mM   20 mM
  LiCl                7500 mM  200 mM
  EDTA                 500 mM    1 mM
""")

lb.show_concs = True
w1.show_concs = True
w2.show_concs = True
ls.show_concs = True

p += pl(f"Prepare {lb.volume} lysis/binding buffer [1]:", lb)
p += pl(f"Prepare {w1.volume} wash buffer I [1]:", w1)
p += pl(f"Prepare {w2.volume} wash buffer II [1]:", w2)
p += pl(f"Prepare {ls.volume} low-salt buffer [1]:", ls)

p.print()
