#!/usr/bin/env python3

"""\
Display a protocol for making the given fragment.

Usage:
    make <fragment_tag>...

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
import shlex
import bind_dna as dbp
from itertools import groupby
from more_itertools import one
from subprocess import run

def format_tx(pcr):
    tx = int(pcr.extension_time_seconds)
    return f'{tx}s' if tx < 60 else f'{tx//60}m{tx%60}'

def join(items):
    items = list(items)
    if len(set(items)) == 1:
        return items[0]
    else:
        return ','.join(items)

args = docopt.docopt(__doc__)
protocols = [
        dbp.get_fragment_protocol(x)
        for x in args['<fragment_tag>']
]

stepwise_cmds = []

for key, group in groupby(protocols, key=lambda x: x.method):
    group = list(group)
    stepwise_cmd = ['stepwise']

    if key == 'PCR':
        stepwise_cmd += [
                'pcr',
                join(x.template_tag for x in group),
                join(x.primer_tags[0] for x in group),
                join(x.primer_tags[1] for x in group),
                str(len(group)),
                join(str(x.annealing_temp_celsius) for x in group),
                join(format_tx(x) for x in group),
        ]


    elif p.method == 'IVT':
        # TODO: Modernize this protocol to take more arguments.
        stepwise_cmd += [
                'ivt',
                str(len(group)),
        ]

    elif p.method == 'RE' and p.enzyme_name == 'XmnI':
        # TODO: Write a general restriction digest protocol.
        stepwise_cmd += [
                'xmni',
                str(len(group)),
        ]

    else:
        inform.warning("{p.method!r} protocols are not yet supported.")
        continue

    stepwise_cmds.append(stepwise_cmd)

if not stepwise_cmds:
    inform.terminate("no protocols found.")

# Hack to deal with stepwise #17:
stepwise_cmds[-1].insert(1, '-x')
stepwise_pipeline = ' | '.join(
        shlex.join(x)
        for x in stepwise_cmds
)

run(stepwise_pipeline, shell=True)

