#!/usr/bin/env python3

"""\
Ligate linker-N to the mRNA.

Usage:
    ligate <n> [-m <reagents>] [-P] [-Q]

Arguments:
    <n>
        The number of reactions to set up.

Options:
    -m --master-mix <reagents>  [default: pnk,lig]
        Include the indicated reagents in the master mix.  The following 
        reagents are understood:

        pnk: T4 PNK
        lig: T4 RNA ligase
        rna: annealed mRNA/linker

    -P --no-pnk
        Leave T4 PNK out of the reaction, e.g. if using phosphorylated linker.

    -Q --no-quench
        Leave out the 65°C incubation to quench the reaction.  This is useful 
        if the reaction will be quenched by a downstream step anyways.
        
"""

import stepwise, docopt
from inform import plural

args = docopt.docopt(__doc__)

ligate = stepwise.MasterMix.from_text("""\
Reagent                  Stock       Volume  MM?
=====================  =======  ===========  ===
water                           to 40.00 µL  yes
BSA                       0.1%       4.0 µL  yes
T4 DNA ligase buffer       10x       4.0 µL  yes
T4 PNK                 10 U/µL      0.33 µL
T4 RNA ligase          40 U/µL       0.5 µL
annealed mRNA/linker   1.25 µM       4.0 µL
""")

ligate.num_reactions = n = eval(args['<n>'])
ligate['T4 PNK'].master_mix = 'pnk' in args['--master-mix']
ligate['T4 RNA ligase'].master_mix = 'lig' in args['--master-mix']
ligate['annealed mRNA/linker'].master_mix = 'rna' in args['--master-mix']

if args['--no-pnk']:
    del ligate['T4 PNK']

protocol = stepwise.Protocol()

protocol += f"""\
Setup {plural(n):# ligation reaction/s}:

{ligate}
"""

protocol += f"""\
Incubate the {plural(n):ligation reaction/s} as follows:

- 25°C for 10 min.
{'' if args['--no-quench'] else '- 65°C for 10 min.'}
"""

print(protocol)
