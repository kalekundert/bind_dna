#!/usr/bin/env python3

"""
Use PCR to verify the insertion of the rpoZ-zif268 fragment.

Usage:
    check_rpoZ_zif268_insert <plasmids>...

Arguments:
    <plasmids>
        The names of the plasmids to check.  It's ok to specify the same 
        plasmid multiple times, e.g. to check multiple colonies.
"""

import stepwise
import docopt

from stepwise import pl, ul
from stepwise_mol_bio import Pcr

if __name__ == '__main__':
    args = docopt.docopt(__doc__)

    plasmids = args['<plasmids>']
    fwd_primer = fwd = 'o2'
    rev_primers = 'o262', 'o185', 'o188'

    n = len(plasmids)
    m = len(rev_primers)

    amplicons_setup = [
            Pcr.Amplicon.from_tags(plasmid, fwd, 'rev')
            for plasmid in plasmids
    ]
    amplicons_thermo = [
            Pcr.Amplicon.from_tags(plasmid, fwd, rev)
            for plasmid in plasmids
            for rev in rev_primers
    ]

    p = stepwise.Protocol()

    p += "Resuspend each colony in 20 µL EB."

    pcr_setup = Pcr(amplicons_setup)
    pcr_setup.only_primer_mix = True
    pcr_setup.force_primer_mix = True
    mm, primer_mix = pcr_setup.reaction

    primer_mix['reverse primer'].name = ','.join(rev_primers)
    primer_mix['forward primer'].master_mix = True
    primer_mix.num_reactions = 3

    vol_rxn = mm.volume
    vol_template = mm['template DNA'].volume

    mm.num_reactions = n
    mm.hold_ratios.volume = vol_rxn * n * 1.2, 'µL'
    mm.volume -= vol_template
    mm['primer mix'].master_mix = False
    del mm['template DNA']

    p += pcr_setup.protocol
    p += pl("Setup 3 PCR master mixes:", mm)

    rxn = stepwise.MasterMix()
    rxn.volume = vol_rxn, 'µL'
    rxn['PCR master mix'].volume = vol_rxn - vol_template, 'µL'
    rxn['template'].volume = vol_template, 'µL'
    rxn['template'].name = ','.join(set(plasmids))

    p += pl(f"Setup {n*m} PCR reactions (every primer/template combo):", rxn)

    pcr_thermo = Pcr(amplicons_thermo)
    pcr_thermo.only_thermocycler = True

    p += pcr_thermo.protocol

    p += pl("Run a 1% E-gel 48:", ul(
        "Add 15 µL water to each reaction.",
        "Load 15 µL in each lane.",
        pl("Ladder:", ul(
            "10 µL water",
            "5 µL 50 ng/µL 1kb+ ladder",
        ), br='\n'),
    ))
    p.print()


