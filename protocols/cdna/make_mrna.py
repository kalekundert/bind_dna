#!/usr/bin/env python3

"""\
Ligate a puromycin linker to mRNA.

Usage:
    make_mrna.py <mrna> <linker> [-n <rxns>]

Arguments:
    <mrna>
        A comma-separated list of mRNA names (e.g. f85).

    <linker>
        A comma-separated list of linker names (e.g. o129).

Options:
    -n --num-reactions <int>
        The number of separate reactions to set up.

    -W --no-wash
        Skip the wash step.
"""

import docopt
import stepwise

args = docopt.docopt(__doc__)
n = args['--num-reactions']

anneal = [
        'cdna/anneal',
        args['<mrna>'],
        args['<linker>'],
        '-n', n,
        '-v', 12,
        '-x', 0.6,
]
ligate = [
        'cdna/ligate',
        '-n', n,
        '-v', 120,
]

p = stepwise.Protocol()
p += stepwise.load(anneal)
p += stepwise.load(ligate)

if not args['--no-wash']:
    p += stepwise.load('cdna/wash')

p += stepwise.load('aliquot "4 µL" "1 µM"')

p.print()
