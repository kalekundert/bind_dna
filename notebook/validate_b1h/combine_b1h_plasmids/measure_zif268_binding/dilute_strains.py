#!/usr/bin/env python3

import numpy as np
import nestedtext as nt
import stepwise

from stepwise import pl, ul, table
from functools import partial

# c1: concentration of diluted cells
# c2: concentration of stock cells, c2 == 3 * c1
# c: desired target concentration
#
# v1: volume of diluted cells, v1 == 120 µL
# v2: volume of stock cells to add
#
# (c1 * v1) + (c2 * v2) = c * (v1 + v2)
# v2 = v1 * (c1 - c) / (c - c2)

strain_ods = nt.load('initial_concs.nt')

c1 = np.array([float(x) for x in strain_ods.values()])
c2 = 3 * c1
c = 0.1

v1 = 120
v2 = v1 * (c1 - c) / (c - c2)

header = [
        'Strain',
        'Conc (OD)',
        'Add (µL)',
]
rows = list(zip(
    strain_ods.keys(),
    c1,
    v2,
))

p = stepwise.Protocol()



p += pl(
        "Make the following dilutions:",
        table(rows=rows, header=header),
)

p.print()
