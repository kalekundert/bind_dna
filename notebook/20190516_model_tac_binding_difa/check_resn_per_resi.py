#!/usr/bin/env python3

"""\
Usage:
    check_pdb.py <pdb>
"""

import docopt
from pathlib import Path

args = docopt.docopt(__doc__)
pdb_path = Path(args['<pdb>'])
pdb_lines = pdb_path.read_text().split('\n')

resis = {}

for line in pdb_lines:
    if not line.startswith('ATOM'):
        continue

    chain = line[21:22]
    resi = int(line[23:26])
    resn = line[17:20].strip()

    k = chain, resi
    resis.setdefault(k, set())
    resis[k].add(resn)

pprint(resis)


