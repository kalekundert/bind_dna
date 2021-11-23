#!/usr/bin/env python3

"""
Usage:
    fluorotext_rnase_digest.py [-V <µL>]

Options:
    -V --reaction-volume <µL>  [default: 10]
        The volume of the IVTT reaction to add the RNase to.
"""

import docopt
import stepwise
from stepwise import pl, ul

args = docopt.docopt(__doc__)
rxn_volume_uL = float(args['--reaction-volume'])

p = stepwise.Protocol()

# RNase volume:
# - The PUREfrex 2.0 manual calls for 1 µL 1 mg/mL RNase A per 10 µL reaction:
#   https://tinyurl.com/ml26pu09
#
# - I don't have any RNase A on hand, but I do have RNase cocktail (Invitrogen 
#   AM2286).  This cocktail contains:
#   - 500 U/mL RNase A
#   - 20,000 U/mL RNase T1
#
# - The manual for the cocktail says to replace RNase A at equivalent RNase A 
#   concentrations.
#
# - I can't easily find a way to convert between mass and units, so in expt 
#   #117 I just decided to use the same proportion as in the PUREfrex manual.   
#   That seems to work well.
#
p += f"Add {rxn_volume_uL / 10:g} µL RNase cocktail (Invitrogen AM2286) to each {rxn_volume_uL:g} µL IVTT reaction."

p += "Incubate at 37°C for 15 min."

p.print()

