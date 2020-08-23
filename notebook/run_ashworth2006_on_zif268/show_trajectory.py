#!/usr/bin/env python3

"""
Usage:
	show_optectory.py <workspace>
"""

import docopt
import json
from pathlib import Path
from pprint import pprint

args = docopt.docopt(__doc__)
work = args['<workspace>']

opt_path = Path(work) / 'optimization.json'
with opt_path.open() as f:
	opt = json.load(f)

for k in range(1, len(opt)+1):
    x = opt[str(k)]
    print(f"{k}\t{x['cst_stdev']}\t{x['objective']}\t{x['score_reu']}")




