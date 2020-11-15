#!/usr/bin/env python3

"""\
Ligate a puromycin linker to mRNA.

Usage:
    make_mrna.py <mrna> <linker> [-n <rxns>]

Options:
    -n --num-reactions <int>    [default: 1]
        The number of separate reactions to set up.

"""

import docopt
import stepwise
import shlex
import po4

args = docopt.docopt(__doc__)
db = po4.load_db()

mrna = args['<mrna>']
linker = args['<linker>']
n = args['--num-reactions']

try:
    mrna_conc = db[mrna].conc_nM / 1000
except po4.QueryError:
    mrna_conc = 10

anneal = [
        'cdna/anneal', n,
        mrna,
        linker,
        '-v', 12,
        '-x', 0.6,
        '-R', mrna_conc,
]
ligate = [
        'cdna/ligate', n,
        '-v', 120,
]

p = stepwise.Protocol()
p += stepwise.load(anneal)
p += stepwise.load(ligate)
p += stepwise.load('cdna/wash')
p += stepwise.load('aliquot "4 µL" "1 µM"')
p.print()
