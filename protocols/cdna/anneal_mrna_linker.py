#!/usr/bin/env python3

"""\
Anneal linker-N and mRNA prior to ligation.

Usage:
    anneal <n> [<mrna>] [<linker>] [-v <µL>] [-m <reagents>] [-x <fold>] [-R <µM>]

Arguments:
    <n>
        The number of reactions to set up.

    <mrna>
        The name of the mRNA, e.g. f11.

    <linker>
        The name of the linker, e.g. o93.

Options:
    -v --volume <µL>                    [default: 4]
        The volume of each annealing reaction in µL.

    -m --master-mix <reagents>          [default: ]
        A comma-separated list of reagents to include in the master mix.  This 
        flag is only relevant in <n> is more than 1.  The following reagents 
        are understood: mrna, link

    -x --excess-linker <fold>           [default: 1]
        The amount of linker to add to the reaction, relative to the amount of 
        mRNA in the reaction.

    -R --mrna-stock <µM>                [default: 10]
        The stock concentration of the mRNA, in µM.  The volume of mRNA will be 
        updated accordingly to keep the amount of material in the reaction 
        constant.
"""

import stepwise, docopt
from inform import plural

args = docopt.docopt(__doc__)

anneal = stepwise.MasterMix.from_text("""\
Reagent  Stock     Volume  MM?
=======  =====  =========  ===
water           to 4.0 µL  yes
PBS      10x       0.4 µL  yes
mRNA     10 µM     0.5 µL
linker   10 µM     0.5 µL
""")

anneal.num_reactions = n = eval(args['<n>'])
anneal.hold_ratios.volume = eval(args['--volume']), 'µL'
anneal['mRNA'].master_mix = 'mrna' in args['--master-mix']
anneal['linker'].master_mix = 'link' in args['--master-mix']

if args['<mrna>']: anneal['mRNA'].name = args['<mrna>']
if args['<linker>']: anneal['linker'].name = args['<linker>']

anneal['mRNA'].hold_conc.stock_conc = int(args['--mrna-stock']), 'µM'
anneal['linker'].volume *= float(args['--excess-linker'])

protocol = stepwise.Protocol()

protocol += f"""\
Setup {plural(n):# annealing reaction/s}:

{anneal}
"""

protocol += f"""\
Perform the {plural(n):annealing reaction/s}:

- Incubate at 95°C for 2 min.
- Cool at room temperature.
"""

print(protocol)


