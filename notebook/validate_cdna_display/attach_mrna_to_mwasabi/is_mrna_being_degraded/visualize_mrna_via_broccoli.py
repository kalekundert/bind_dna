#!/usr/bin/env python3

import stepwise
import po4

db = po4.load_db()
mrna = 'f97'

p = stepwise.Protocol()

# Deciding how much to load:
# - I want one lane to be almost maximally concentrated, just to make sure I 
#   see something.
# - Normally I load 200 ng on urea gels, so I should have at least one lane 
#   close to that.
# - [Filonov2015] claims to be able to detect 100 pg of RNA (Figure 2), so that 
#   makes sense as a lower bound.
p += stepwise.load([
    'serial',
    '2 µL', f'{db[mrna].conc_ng_uL} ng/µL', 0.1, 8,
    '-m', 'f97',
    '-d', 'nuclease-free water',
])
p += stepwise.load([
    'gel', 'urea', mrna,
    '-S',
])
p += stepwise.load('stain_dfhbi')

p.print()

