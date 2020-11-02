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

args = docopt.docopt(__doc__)

anneal = [
        'cdna/anneal',
        args['--num-reactions'],
        args['<mrna>'],
        args['<linker>'],
        '-v', 8,
        '-x', 0.6,
        '-R', 5,
]
ligate = [
        'cdna/ligate',
        args['--num-reactions'],
        '-v', 80,
]

p = stepwise.Protocol()
p += stepwise.load(shlex.join(str(x) for x in anneal))
p += stepwise.load(shlex.join(str(x) for x in ligate))
p += stepwise.load('cdna/wash')
p += stepwise.load('aliquot "4 µL" "333 nM"')
print(p)

