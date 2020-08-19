#!/usr/bin/env python3

"""\
Anneal linker-N and mRNA prior to ligation.

Usage:
    anneal <n> [<mrna>] [<linker>] [-v <µL>] [-m <reagents>] [-L <conc>]

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
        The reagents to include in the master mix.  The following reagents are 
        understood: 'mrna' and 'link'.  To specify both reagents, separate the 
        two names with a comma.

    -L --linker-stock <conc>
        The stock concentration of the linker.  Note that this setting does not 
        change the volume of the reaction, so it does change the amount of 
        linker in the reaction.  Not unit is assumed; the unit should be 
        specified.
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

if args['--linker-stock']:
    anneal['linker'].stock_conc = args['--linker-stock']

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


