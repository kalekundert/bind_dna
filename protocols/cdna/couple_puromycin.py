#!/usr/bin/env python3

"""\
Describe how to setup the puromycin/peptide coupling reaction.

Usage:
    couple_puromycin [-n <rxns>] [-v <µL>] [-m <mM>] [-M <mM>] [-k <mM>]
        [-K <mM>] [-t <time>] [-T <°C>]

Options:
    -n --num-reactions <int>        [default: 1]
        The number of reactions to setup.

    -v --volume <µL>                [default: 5]
        The volume of the expression reaction.

    -m --mg-conc <mM>               [default: 65]
        The desired final concentration of Mg²⁺, in mM.  The default is from 
        [Naimudden2016].

    -M --mg-stock <mM>              [default: 1000]
        The stock concentration of Mg²⁺, in mM.

    -k --k-conc <mM>                [default: 750]
        The desired final concentration of K⁺, in mM.  The default is from 
        [Naimudden2016].

    -K --k-stock <mM>               [default: 3000]
        The stock concentration of K⁺, in mM.

    -t --incubate-time <time>       [default: 1h]
        How long to incubate the reaction for.  No unit is assumed.

    -T --incubate-temp <°C>         [default: 25]
        What temperature to incubate the reaction at, in °C.
"""

import docopt
import stepwise
import numpy as np

args = docopt.docopt(__doc__)
rxn_uL = float(args['--volume'])
mg_mM = float(args['--mg-conc'])
mg_stock_mM = float(args['--mg-stock'])
k_mM = float(args['--k-conc'])
k_stock_mM = float(args['--k-stock'])

M = np.array([
    [mg_mM, -mg_stock_mM,           0,      0],
    [ k_mM,            0, -k_stock_mM,      0],
    [    1,           -1,          -1, rxn_uL],
])
salt_uL = np.linalg.solve(M[:,:3], M[:,3])

rxn = stepwise.MasterMix()
rxn.num_reactions = int(args['--num-reactions'])
rxn.extra_min_volume = 5, 'µL'
rxn.solvent = None
rxn['PURExpress'].volume = rxn_uL, 'µL'
rxn['MgOAc'].stock_conc = mg_stock_mM, 'mM'
rxn['MgOAc'].volume = salt_uL[1], 'µL'
rxn['MgOAc'].master_mix = True
rxn['KCl'].stock_conc = k_stock_mM, 'mM'
rxn['KCl'].volume = salt_uL[2], 'µL'
rxn['KCl'].master_mix = True

p = stepwise.Protocol()
p += f"""\
Setup the coupling reaction:

{rxn}
"""
p += f"""\
Incubate at {args['--incubate-temp']}°C for {args['--incubate-time']}.
"""

print(p)
