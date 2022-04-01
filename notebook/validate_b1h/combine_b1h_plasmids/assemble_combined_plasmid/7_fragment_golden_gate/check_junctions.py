#!/usr/bin/env python3

"""
Check the p170 assembly by colony PCR

Usage:
    check_junctions.py <n>

Arguments:
    <n>
        The number of colonies to test.
"""

import stepwise
import docopt

from stepwise import pl, dl
from stepwise_mol_bio import Pcr, Gel
from more_itertools import flatten
from inform import plural

if __name__ == '__main__':
    args = docopt.docopt(__doc__)

    n = int(args['<n>'])
    fwd_primer = 'o2'
    rev_primers = 'o88', 'o185', 'o188'

    amplicons = list(flatten(
            [
                Pcr.Amplicon.from_tags('p170', fwd_primer, rev_primer)
                for rev_primer in rev_primers
            ]
            for i in range(n)
    ))
    pcr = Pcr(amplicons)
    pcr.master_mix = {'fwd'}
    pcr.reaction[0]['template DNA'].name = 'p170 colony'
    pcr.reaction[0]['template DNA'].stock_conc = None

    gel = Gel('agarose/1', n * len(rev_primers))

    p = stepwise.Protocol()
    p += f"Pick {plural(n):# colon/y/ies} and resuspend each in 20 µL water."
    p += pcr.protocol
    p += gel.protocol

    p += pl(
            "Check for bands of the following sizes:",
            dl(
                ('o88', '2.1 kb'),
                ('o185', '2.2 kb'),
                ('o188', '3.2 kb'),
            ),
    )

    p.print()
