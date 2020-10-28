#!/usr/bin/env python3

import stepwise
import po4

p = stepwise.Protocol()

p += """\
In the following steps, setup these reactions:

Buffer:      Mn  Mn  Mn  Mn  Mn  Mn  Mg  Mg  Mg  Mg  Mg  Mg  
dCas9:       +   −   +   −   +   +   +   −   +   −   +   +
DNA:         −   +   +   +   +   +   −   +   +   +   +   +
HUH-target:  −   −   −   +   +   +   −   −   −   +   +   +
EDTA:        −   −   −   −   −   +   −   −   −   −   −   +
"""

p += stepwise.load('huh/huh_cas9 f99 f100 -n 4 -S -b "Mn²⁺/Mg²⁺ buffer" -B 2.5')
p += stepwise.load('gel sdsmax 12 -S')
p += stepwise.load('gelgreen -r')
p += stepwise.load('simplyblue -fr')

print(p)


