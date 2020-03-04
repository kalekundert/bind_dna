#!/usr/bin/env python3

"""\
Display a protocol for making the given fragment.

Usage:
    make <tag>... [--] [<options>...]

Arguments:
    <tag>
        The name of a DNA/RNA construct listed in one of the following tables:

            sequences/plasmids.xlsx
            sequences/fragments.xlsx

    <options>
        Options that will be passed directly to all protocols invoked.

The protocol is parsed from the "Construction" column of these tables.  The 
syntax for a protocol is:

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

GG: Golden gate assembly
    bb=<backbone plasmid (required)>
    ins=<comma-separated list of inserts (required)>
    enzyme=<Type IIS enzyme, will default to BsaI if not given>
"""

import sys, re
import shlex
import docopt
import bind_dna as dbp
from itertools import groupby
from more_itertools import one
from subprocess import run
from statistics import mean

def make_pcr_command(group):
    def get_tx(pcr):
        tx = int(pcr.extension_time_seconds)
        return f'{tx}s' if tx < 60 else f'{tx//60}m{tx%60}'

    def get_scale(group):
        custom_scale = one(
                {x.scale for x in group},
                too_long=ValueError(f"PCR reactions have different scales: {','.join(repr(x.product_tag) for x in group)}"),
        )
        # Default to 50 µL, because usually if we're trying to *make* 
        # something, we need it at >10 µL scale.
        return custom_scale or 50

    def get_master_mix(group):
        master_mix = {'dna', 'primers'}
        num_templates = len({x.template_tag for x in group})
        num_primers = len({x.primer_tags for x in group})

        if num_templates > 1:
            master_mix.discard('dna')
        if num_primers > 1:
            master_mix.discard('primers')

        return master_mix

    return [
            'pcr',
            join(x.template_tag for x in group),
            join(x.primer_tags[0] for x in group),
            join(x.primer_tags[1] for x in group),
            str(len(group)),
            join(str_strip_insig(x.annealing_temp_celsius) for x in group),
            join(get_tx(x) for x in group),
            '-m', ','.join(get_master_mix(group)),
            '-v', str(get_scale(group)),
    ]

def make_ivt_command(group):
    return [
            'ivt',
            str(len(group)),
    ]

def make_digest_command(group):
    return [
            'xmni',
            str(len(group)),
            '-p', join(x.template_tag for x in group),
    ]

def make_gg_command(group):
    def fragment_from_tag(tags):
        def get_conc(tag):
            try:
                conc = dbp.get_conc_str(tag)
                return re.sub(r'\s*ng/[µu]L$', '', conc)
            except ValueError:
                return '50'

        def tabulate(x):
            from textwrap import indent
            from tabulate import tabulate
            return indent(
                    tabulate(x.items(), tablefmt='plain'),
                    '    ',
            )

        name = join(tags)

        concs = {
                x: get_conc(x)
                for x in tags
        }
        conc = one(
                set(concs.values()),
                too_long=ValueError(f"inserts have different concentrations:\n{tabulate(concs)}"),
        )

        lengths = {
                x: len(dbp.get_seq(x))
                for x in tags
        }
        if min(lengths.values()) > 0.5 * max(lengths.values()):
            length = int(mean(lengths.values()))
        else:
            raise ValueError(f"inserts have lengths that differ by more than 50%:\n{tabulate(lengths)}")

        return f'{name}:{conc}:{length}'

    def fragment_names():
        from itertools import count
        yield 'bb'
        for i in count(1):
            yield str(i)

    num_inserts = one(
            {len(x.insert_tags) for x in group},
            too_long=ValueError("All golden gate assemblies must have the same number of inserts"),
    )
    tag_groups = [
            [x.backbone_tag for x in group]
    ] + [
            [x.insert_tags[i] for x in group]
            for i in range(num_inserts)
    ]
    fragments = [
            fragment_from_tag(x)
            for x in tag_groups
    ]

    master_mix = []
    for name, tags in zip(fragment_names(), tag_groups):
        if len(set(tags)) == 1:
            master_mix.append(name)

    stepwise_cmd = [
            'golden_gate',
            *fragments,
            '-e', join(x.enzyme_name for x in group),
    ]
    if (n := len(group)) > 1:
        stepwise_cmd += [
            '-n', str(n),
            '-m', ','.join(master_mix),
        ]

    return stepwise_cmd


def join(items):
    items = list(items)
    if len(set(items)) == 1:
        return items[0]
    else:
        return ','.join(items)

def str_strip_insig(x):
    return str(x).rstrip('0').rstrip('.')

try: 
    i = sys.argv.index('--')
    args = docopt.docopt(__doc__, sys.argv[1:i])
    args['<options>'] = sys.argv[i:][1:]
except ValueError:
    i = -1
    args = docopt.docopt(__doc__)
    args['<options>'] = []


# Instead of looping through protocol objects, I should be using "Construct" 
# object that represent all the columns in the database.  That would let me do 
# things like scale down for cloning reaction.

protocols = []
for tag in args['<tag>']:
    protocol = dbp.get_protocol(tag)
    protocol.product_tag = tag
    protocols.append(protocol)

stepwise_cmds = []
for key, group in groupby(protocols, key=lambda x: x.method):
    group = list(group)
    stepwise_cmd = ['stepwise']

    if key == 'PCR':
        stepwise_cmd += make_pcr_command(group)

    elif key == 'IVT':
        # TODO: Modernize this protocol to take more arguments.
        stepwise_cmd += make_ivt_command(group)

    elif key == 'RE' and join(x.enzyme_name for x in group) == 'XmnI':
        # TODO: Write a general restriction digest protocol.
        stepwise_cmd += make_digest_command(group)

    elif key == 'GG':
        stepwise_cmd += make_gg_command(group)

    else:
        inform.warning("{key!r} protocols are not yet supported.")
        continue

    stepwise_cmd += args['<options>']
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

