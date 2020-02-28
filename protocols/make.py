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

def str_strip_insig(x):
    return str(x).rstrip('0').rstrip('.')

def pcr_scale(group):
    custom_scale = one(
            {x.scale for x in group},
            too_long=ValueError(f"PCR reactions have different scales: {','.join(repr(x.product_tag) for x in group)}"),
    )
    # Default to 50 µL, because usually if we're trying to *make* something, we 
    # need it at >10 µL scale.
    return custom_scale or 50

def pcr_master_mix(group):
    master_mix = {'dna', 'primers'}
    num_templates = len({x.template_tag for x in group})
    num_primers = len({x.primer_tags for x in group})

    if num_templates > 1:
        master_mix.discard('dna')
    if num_primers > 1:
        master_mix.discard('primers')

    return master_mix

args = docopt.docopt(__doc__)

protocols = []
for tag in args['<fragment_tag>']:
    protocol = dbp.get_fragment_protocol(tag)
    protocol.product_tag = tag
    protocols.append(protocol)

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
                join(str_strip_insig(x.annealing_temp_celsius) for x in group),
                join(format_tx(x) for x in group),
                '-m', ','.join(pcr_master_mix(group)),
                '-v', str(pcr_scale(group)),
        ]

    elif key == 'IVT':
        # TODO: Modernize this protocol to take more arguments.
        stepwise_cmd += [
                'ivt',
                str(len(group)),
        ]

    elif key == 'RE' and join(x.enzyme_name for x in group) == 'XmnI':
        # TODO: Write a general restriction digest protocol.
        stepwise_cmd += [
                'xmni',
                str(len(group)),
                '-p', join(x.template_tag for x in group),
        ]

    else:
        inform.warning("{key!r} protocols are not yet supported.")
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

