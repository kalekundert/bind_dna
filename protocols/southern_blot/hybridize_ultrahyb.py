#1/usr/bin/env python3

"""\
Hybridize a probe to a Southern blot membrane.

Usage:
  hybridize_ultrahyb <probe>
"""

import docopt
from stepwise import Protocol

args = docopt.docopt(__doc__)
probe = args['<probe>']

p = Protocol()

p += f"""\
Hybridize the probe to the blotted DNA/RNA [1]:

- Preheat ULTRAhyb-Oligo to 68°C.  Make sure that 
  any precipitated material has redissolved.

- Add enough ULTRAhyb-Oligo to keep the membrane 
  uniformly wet (≈1 mL per 10 cm² membrane).

- Incubate at 42°C for 30 min.

- Add {args['<probe>']} to a concentration of 1 nM [2,3].

- Incubate at 42°C for 14-24 h.

- Discard probe buffer.

- Repeat twice:
  - Add 50 mL wash buffer: 2x SSC, 0.5% SDS
  - Incubate at 42°C for 30 min.

- Allow to dry.
"""
p.footnotes[1] = """\
ULTRAhyb-Oligo manual:
https://tinyurl.com/yxd8tpgt

This product is meant for hybridizing DNA probes 
to RNA, not DNA.  But I think it's better than 
using either Denhardt's solution (no nucleic acid 
blocking agents, e.g. ssDNA) or ULTRAhyb (not 
recommended for oligonucleotide probes).
"""
p.footnotes[2] = """\
Miller et al (2018). doi:10.1261/rna.068213.118
"""
p.footnotes[3] = """\
Make sure that the undiluted probe never touches 
the membrane directly.  If necessary, pour the 
blocking solution into a clean tube, add the 
oligo, mix, then pour back over the membrane.
"""
p.print()
