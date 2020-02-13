#!/usr/bin/env python3

"""\
Ligate linker-N to the mRNA.

Usage:
    ligate <n> [-m <reagents>]

Arguments:
    <n>
        The number of reactions to set up.

    -m --master-mix <reagents>  [default: pnk,lig]
"""

import stepwise, docopt
from inform import plural

args = docopt.docopt(__doc__)

ligate = stepwise.MasterMix.from_text("""\
Reagent                  Stock    Volume  MM?
=====================  =======  ========  ===
water                           27.17 µL  yes
T4 DNA ligase buffer       10x    4.0 µL  yes
BSA                       0.1%    4.0 µL  yes
T4 PNK                 10 U/µL   0.33 µL
T4 RNA ligase          40 U/µL    0.5 µL
annealed mRNA/linker   1.25 µM    4.0 µL
""")

ligate.num_reactions = n = eval(args['<n>'])
ligate['T4 PNK'].master_mix = 'pnk' in args['--master-mix']
ligate['T4 RNA ligase'].master_mix = 'lig' in args['--master-mix']
ligate['annealed mRNA/linker'].master_mix = 'rna' in args['--master-mix']

protocol = stepwise.Protocol()

protocol += f"""\
Setup {plural(n):# ligation reaction/s}:

{ligate}
"""

protocol += f"""\
Incubate the {plural(n):ligation reaction/s} as follows:

- 25°C for 10 min.
- 65°C for 10 min.
"""

print(protocol)
