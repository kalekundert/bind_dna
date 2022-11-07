#!/usr/bin/env python3

# Use antisense primers as control.
#
# Conditions:
# - primers: sense, antisense
# - target: TGG, AAA
# - dpni: +, -
# 
# grow, rna miniprep, RT, qPCR

import stepwise

from stepwise_mol_bio import Grow, ReverseTranscribe, Dnase
from stepwise_mol_bio.dnase import dnase_digest
from stepwise_mol_bio.reverse_transcribe import reverse_transcribe
from stepwise_mol_bio.future.pcr import Pcr, pcr
from stepwise import pl, ul
from freezerbox import load_db
from itertools import product

db = load_db()
pm = '−', '+'

p = stepwise.Protocol()

p += stepwise.load('grow -p "b1h/qpcr" sz234 sz235')
p += stepwise.load('zymo_quick_rna')

# DNase:

dpni_samples = []
dpni_templates = 'sz234', 'sz235'
dpni_presets = 'thermo/ezdnase', 'neb/dpni'
dpni_preset_names = {
        'thermo/ezdnase': 'ezDNase',
        'neb/dpni': 'DpnI',
}

for template, preset in product(dpni_templates, dpni_presets):
    sample = Dnase.Sample(
            template,
            preset=preset,
            db=db,
    )
    rxn = sample.reaction_prototype
    rxn['RNA'].volume = rxn.volume / 2
    rxn.hold_ratios.volume = 10 * 4.4, 'µL'
    dpni_samples.append(sample)

p += dnase_digest(dpni_samples)

# RT:

# I wanted to include the DNase reaction in the RT reaction, but the code 
# doesn't understand that DNase reactions with different parameters are in fact 
# different and can't be mixed together.  This probably wouldn't be too hard to 
# fix, but for now this is the path of least resistance.

rt_samples = []
rt_primers = [
        {'primer': 's1+s3', 'gene_specific_primer': True},
        {'primer': 's2+s5', 'gene_specific_primer': True},
]

for dpni, primer in product(dpni_samples, rt_primers):
    Sample = ReverseTranscribe.partial(
            f'{dpni.rna_name}, {dpni_preset_names[dpni.preset]}',
            **primer,
            preset='thermo/superscript-iv',
            db=db,
    )
    rt_samples += [
            Sample(include_rt=True),
            Sample(include_rt=False),
    ]

p += reverse_transcribe(rt_samples)

# qPCR:

p += pl(
        "Dilute the reverse-transcribed DNA 10x:",
        ul(
            "18 µL water",
            "2 µL RT reaction",
        ),
)

pcr_samples = []
pcr_primers = ('s1','s2'), ('s3','s5')

for rt_sample, (fwd, rev) in product(rt_samples, pcr_primers):
    pcr_sample = Pcr.from_tags(
            f'{rt_sample.template}, {rt_sample.primer}, {pm[rt_sample.include_rt]}RT',
            fwd, rev,
            preset='ssoadv',
            product_length_bp=74,
            db=db,
    )
    pcr_samples.append(pcr_sample)

p += pcr(3 * pcr_samples)

p.print()
