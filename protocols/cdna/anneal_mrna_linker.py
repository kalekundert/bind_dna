#!/usr/bin/env python3

"""\
Anneal linker-N and mRNA prior to ligation.

Usage:
    anneal [<mrna>] [<linker>] [-n <int>] [-v <µL>] [-m <reagents>] [-x <fold>] 
        [-R <µM>] [-L <conc>]

Arguments:
    <n>
        The number of reactions to set up.

    <mrna>
        The name of the mRNA, e.g. f11.  Multiple comma-separated names may be 
        given.

    <linker>
        The name of the linker, e.g. o93.  Multiple comma-separated names may 
        be given.

Options:
    -n --num-reactions <int>
        The number of reactions to setup up.  The default is the number of 
        unique combinations of mRNA and linker.

    -v --volume <µL>                    [default: 4]
        The volume of each annealing reaction in µL.

    -m --master-mix <reagents>          [default: ]
        A comma-separated list of reagents to include in the master mix.  This 
        flag is only relevant in <n> is more than 1.  The following reagents 
        are understood: mrna, link

    -x --excess-linker <fold>           [default: 1]
        The amount of linker to add to the reaction, relative to the amount of 
        mRNA in the reaction.

    -R --mrna-stock <µM>
        The stock concentration of the mRNA, in µM.  The volume of mRNA will be 
        updated accordingly to keep the amount of material in the reaction 
        constant.  The default is read from the PO₄ database, or 10 µM if the 
        given mRNA is not in the database.

    -L --linker-stock <conc>
        The stock concentration of the linker, in user-specified units.  The 
        volume of linker will not be updated, so this will change the relative 
        proportion of linker to mRNA.
"""

import stepwise
import docopt
from inform import plural
import po4

args = docopt.docopt(__doc__)
db = po4.load_db()
mrnas = args['<mrna>'].split(',') if args['<mrna>'] else []
linkers = args['<linker>'].split(',') if args['<linker>'] else []

anneal = stepwise.MasterMix.from_text("""\
Reagent              Stock     Volume  MM?
===================  =====  =========  ===
nuclease-free water         to 4.0 µL  yes
PBS                  10x       0.4 µL  yes
mRNA                 10 µM     0.5 µL
linker               10 µM     0.5 µL
""")

anneal.num_reactions = n = int(
        args['--num-reactions'] or 
        len(mrnas) * len(linkers)
)
anneal.hold_ratios.volume = eval(args['--volume']), 'µL'
anneal['mRNA'].master_mix = 'mrna' in args['--master-mix']
anneal['linker'].master_mix = 'link' in args['--master-mix']

if mrnas: anneal['mRNA'].name = ','.join(mrnas)
if linkers: anneal['linker'].name = ','.join(linkers)

def get_conc_uM(tag, override):
    if override:
        return float(override), 'µM'
    try:
        return db[tag].conc_nM / 1000, 'µM'
    except po4.QueryError:
        return 10, 'µM'

def consensus(values):
    from itertools import groupby
    from more_itertools import one
    return one((x[0] for x in groupby(values)))

anneal['mRNA'].hold_conc.stock_conc = consensus(
        get_conc_uM(x, args['--mrna-stock']) for x in mrnas
)
anneal['linker'].hold_conc.stock_conc = consensus(
        get_conc_uM(x, args['--linker-stock']) for x in linkers
)
anneal['linker'].volume *= float(args['--excess-linker'])

protocol = stepwise.Protocol()

protocol += f"""\
Setup {plural(n):# annealing reaction/s} [1]:

{anneal}
"""

protocol.footnotes[1] = """\
Using 0.6x linker reduces the amount of unligated 
linker, see expt #1.
"""

protocol += f"""\
Perform the {plural(n):annealing reaction/s}:

- Incubate at 95°C for 2 min.
- Cool at room temperature.
"""

print(protocol)


