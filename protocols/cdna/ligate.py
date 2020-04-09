#!/usr/bin/env python3

"""\
Ligate linker-N to the mRNA.

Usage:
    ligate <n> [-v <µL>] [-x <percent>] [-m <reagents>] [-M] [-P] [-Q]

Arguments:
    <n>
        The number of reactions to set up.

Options:
    -v --volume <µL>                [default: 40]
        The volume of each ligation reaction in µL.

    -x --extra <percent>            [default: 10]
        How much extra master mix to prepare.

    -m --master-mix <reagents>      [default: pnk,lig]
        Include the indicated reagents in the master mix.  The following 
        reagents are understood:

        pnk: T4 PNK
        lig: T4 RNA ligase
        rna: annealed mRNA/linker

    -M --no-master-mix
        Exclude all optional reagents from the master mix, i.e. `-m ''`.

    -P --no-pnk
        Leave T4 PNK out of the reaction, e.g. if using phosphorylated linker.

    -Q --no-quench
        Leave out the 65°C incubation to quench the reaction.  This is useful 
        if the reaction will be quenched by a downstream step anyways.
        
"""

import stepwise, docopt
from inform import plural

args = docopt.docopt(__doc__)
master_mix = '' if args['--no-master-mix'] else args['--master-mix']

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
ligate.extra_percent = eval(args['--extra'])
ligate.hold_ratios.volume = eval(args['--volume']), 'µL'
ligate['T4 PNK'].master_mix = 'pnk' in master_mix
ligate['T4 RNA ligase'].master_mix = 'lig' in master_mix
ligate['annealed mRNA/linker'].master_mix = 'rna' in master_mix

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
