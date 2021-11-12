#!/usr/bin/env python3

"""
Use PCR to verify the insertion of the rpoZ-zif268 fragment.

Usage:
    check_rpoZ_zif268_insert <plasmids>... [-p] [-r <primers>]

Arguments:
    <plasmids>
        The names of the plasmids to check.  It's ok to specify the same 
        plasmid multiple times, e.g. to check multiple colonies.  
        Alternatively, you can specify "N*plasmid" to repeat the given plasmid 
        the given number of times (you may have to make sure that the '*' isn't 
        interpreted as a glob).

Options:
    -p --plasmid
        Use purified plasmid as the template (instead of bacterial colonies).

    -r --primers <names>        [default: o262,o185,o188]
        The names of the reverse primers to use, separated by commas.  Each 
        primer must be present in the freezerbox database.  Below are the 
        primers that I've designed for this experiment:

        o262: between the backbone and the barcode (SR045)
        o185: between the barcode and Zif268 (SR022)
        o188: between Zif268 and the terminator (SR151)
        o266: inside Zif268

"""

import stepwise
import docopt

from stepwise import pl, ul
from stepwise_mol_bio import Pcr
from more_itertools import one, flatten, unique_everseen as unique
from collections import Counter
from inform import plural

def parse_plasmids(plasmid_strs):
    return list(flatten(parse_plasmid(x) for x in plasmid_strs))

def parse_plasmid(plasmid_str):
    fields = plasmid_str.split('*')

    if len(fields) == 1:
        return fields
    if len(fields) == 2:
        return int(fields[0]) * fields[1:]

    raise ValueError(f"expected <plasmid> or <N*plasmid>, got: {plasmid_str}")

def parse_primers(primers_str):
    return primers_str.split(',')

if __name__ == '__main__':
    args = docopt.docopt(__doc__)
    use_plasmid = args['--plasmid']

    plasmids = parse_plasmids(args['<plasmids>'])
    fwd_primer = fwd = 'o2'
    rev_primers = parse_primers(args['--primers'])

    n = len(plasmids) + 1  # +1: negative control
    m = len(rev_primers)

    p = stepwise.Protocol()

    if not use_plasmid:
        p += "Resuspend each colony in 20 µL EB."

    if m == 1:
        # If there is only one reverse primer, the default PCR protocol does a 
        # good job making a master mix.
        rev = one(rev_primers)
        amplicons = [
                Pcr.Amplicon.from_tags(plasmid, fwd, rev)
                for plasmid in plasmids
        ]
        pcr = Pcr(amplicons)

        counts = Counter(plasmids)
        names = [
                f'{xn}×{x}' if (xn := counts[x]) > 1 else x
                for x in counts
        ]
        mm, primer_mix = pcr.reaction
        mm['template DNA'].name = ','.join(names)

        p += pcr.protocol

    else:
        # If there are multiple reverse primers, explicitly describe how to 
        # setup a master mix for each one.
        amplicons_setup = [
                Pcr.Amplicon.from_tags(plasmid, fwd, 'rev')
                for plasmid in plasmids
        ]
        amplicons_thermo = [
                Pcr.Amplicon.from_tags(plasmid, fwd, rev)
                for plasmid in plasmids
                for rev in rev_primers
        ]

        pcr_setup = Pcr(amplicons_setup)
        pcr_setup.only_primer_mix = True
        pcr_setup.force_primer_mix = True
        mm, primer_mix = pcr_setup.reaction

        primer_mix['reverse primer'].name = ','.join(rev_primers)
        primer_mix['forward primer'].master_mix = True
        primer_mix.num_reactions = m

        vol_rxn = mm.volume
        vol_template = mm['template DNA'].volume

        mm.num_reactions = m
        mm.hold_ratios.volume = vol_rxn * n * 1.2, 'µL'
        mm.volume -= vol_template
        mm['primer mix'].master_mix = False
        del mm['template DNA']

        p += pcr_setup.protocol
        p += pl(f"Setup {plural(m):# PCR master mix/es}:", mm)

        rxn = stepwise.MasterMix()
        rxn.volume = vol_rxn, 'µL'
        rxn['PCR master mix'].volume = vol_rxn - vol_template, 'µL'
        rxn['template'].volume = vol_template, 'µL'
        rxn['template'].name = ','.join(['Ø', *unique(plasmids)])
        if use_plasmid:
            rxn['template'].stock_conc = 20, 'pg/µL'

        p += pl(f"Setup {n*m} PCR reactions (every primer/template combo, plus a negative control):", rxn)

        pcr_thermo = Pcr(amplicons_thermo)
        pcr_thermo.only_thermocycler = True

        p += pcr_thermo.protocol

    if n*m <= 10:
        p += pl("Run a 1% E-gel EX.")
    else:
        egel_quickstart = 'https://tinyurl.com/38nmnmk5'
        p += pl(f"Run a 1% E-gel 48{p.add_footnotes(egel_quickstart)}:", ul(
            pl("Prepare the samples:", ul(
                "5 µL PCR reaction",
                "15 µL water"
            ), br='\n'),
            pl("Prepare the ladder:", ul(
                "10 µL water",
                "10 µL 50 ng/µL 1kb+ ladder",
            ), br='\n'),
            "Load 15 µL of each sample.",
            "Run for 20 min.",
            "Visualize with a blue-light transilluminator.",
        ))

    p.print()


