#!/usr/bin/env python3

"""\
Attach the puromycin oligo to the RT/ligation oligo by click chemistry.

Usage:
    click <puro-oligo> <anneal-oligo> [-v <µL>] [-I]

Options:
    -v --volume <µL>
        The volume of the reaction, in µL.
    
    -I --no-incubation
        Skip the incubation step (e.g. so it can be specified differently in a 
        later step).
"""

import docopt
import stepwise

args = docopt.docopt(__doc__)

# PBS:
# - I don't remember why I decided to add PBS to this reaction.
# - The Glen Research protocol referenced by expt #57 doesn't call for any 
#   salt.
rxn = stepwise.MasterMix("""\
  Reagent                  Stock  Volume
  ======================  ======  ======
  puromycin oligo         400 µM    1 µL
  annealing oligo         400 µM    1 µL
  PBS                         2x    2 µL
""")

rxn['puromycin oligo'].name = args['<puro-oligo>']
rxn['annealing oligo'].name = args['<anneal-oligo>']

if v := args['--volume']:
    rxn.hold_ratios.volume = v, 'µL'

p = stepwise.Protocol()

p += stepwise.Step(
        "Couple the linker oligos:",
        rxn,
)

# −20°C incubation:
# - I haven't tested this yet.
# - In principle, it would work by increasing the concentration of the oligos 
#   or the salts (or both) as the solvent crystallizes.
if not args['--no-incubation']:
    p += stepwise.Step(
        "Incubate in the dark as follows:",
        substeps=[
            "room temperature for >4h",
            "−20°C overnight",
        ],
    )

    p += "Dilute 10x to 10 µM."

p.print()
