#!/usr/bin/env python3

"""\
Display a protocol for making the given fragment.

Usage:
    make <fragment_tag>

The fragment tags refer to a row in `sequences/fragments.xlsx`.  The protocol 
is parsed from the "Construction" column of that table.  The syntax for a protocol is:

    <method>: [<key>=<value> ]...

The following methods are currently understood:

PCR: Polymerase chain reaction
    template=<name of template (required)>
    primers=<forward and reverse primers, separated by a comma (required)>
    Ta=<annealing temperature, will derive from primers if not given>

RE: Restriction enzyme digest
    template=<name of template (required)>
    enzyme=<name of enzyme (required)>

IVT: In vitro transcription of DNA.
    template=<name of DNA template (required)>
"""

import docopt
import bind_dna as dbp
from subprocess import run

args = docopt.docopt(__doc__)
p = dbp.get_fragment_protocol(args['<fragment_tag>'])

stepwise_cmd = ['stepwise', '-x']

if p.method == 'PCR':
    tx = int(p.extension_time_seconds)

    stepwise_cmd += [
            'pcr',
            p.template_tag,
            *p.primer_tags,
            '1',
            str(p.annealing_temp_celsius),
            f'{tx}s' if tx < 60 else f'{tx//60}m{tx%60}'
    ]

elif p.method == 'IVT':
    # TODO: Modernize this protocol to take more arguments.
    stepwise_cmd += ['ivt', '1']

elif p.method == 'RE' and p.enzyme_name == 'XmnI':
    # TODO: Write a general restriction digest protocol.
    stepwise_cmd += ['xmni']

else:
    inform.terminate("{p.method!r} protocols are not yet supported.")

run(stepwise_cmd)



