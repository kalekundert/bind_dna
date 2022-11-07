#!/usr/bin/env python3

import stepwise

from stepwise_mol_bio import Grow, ReverseTranscribe
from stepwise_mol_bio.reverse_transcribe import reverse_transcribe
from stepwise_mol_bio.future.pcr import Pcr, pcr
from stepwise import pl, ul
from freezerbox import load_db

# 2022/07/21:
# I've just compared the RT enzymes, and saw that SuperScript IV performed 
# markedly better than any of the others.  With this experiment, I want to (i) 
# replicate that result and (ii) see if the choice of primer is important.

db = load_db()

p = stepwise.Protocol()

p += stepwise.load('grow -p "b1h/qpcr" sz224 sz228')
p += stepwise.load('zymo_quick_rna')

# RT:

# Primers:
# - random hexamers
# - oligo dT
# - custom: s2+s5
# - custom, but wrong direction: s1+s3

# TODO:
# - Fix ezDNase volume bug

rt_primers = [
        {'primer': 'dt'},
        {'primer': 'hex'},
        {'primer': 's1+s3', 'gene_specific_primer': True},
        {'primer': 's2+s5', 'gene_specific_primer': True},
]

rt_templates = 'sz224', 'sz228'
rt_samples = []

for template in rt_templates:
    for primer in rt_primers:
        Sample = ReverseTranscribe.partial(
                template,
                **primer,
                preset='thermo/superscript-iv/ezdnase',
                db=db,
        )
        rt_samples += [
                Sample(include_rt=True),
                Sample(include_rt=False),
        ]

p += reverse_transcribe(rt_samples)

# PCR:

p += pl(
        "Dilute the reverse-transcribed DNA 10x:",
        ul(
            "18 µL water",
            "2 µL RT reaction",
        ),
)

pcr_samples = []
pcr_primers = ('s1','s2'), ('s3','s5')
rt_controls = ('−RT', '+RT')

for rt_sample in rt_samples:
    for fwd, rev in pcr_primers:
        pcr_sample = Pcr.from_tags(
                f'{rt_sample.template}, {rt_sample.primer}, {rt_controls[rt_sample.include_rt]}',
                fwd, rev,
                preset='luna',
                product_length_bp=74,
                db=db,
        )
        pcr_samples.append(pcr_sample)

p += pcr(3 * pcr_samples)

p.print()


