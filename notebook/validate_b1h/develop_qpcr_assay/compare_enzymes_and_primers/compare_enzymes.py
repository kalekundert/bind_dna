#!/usr/bin/env python3

import stepwise

from stepwise_mol_bio import Grow, ReverseTranscribe
from stepwise_mol_bio.reverse_transcribe import reverse_transcribe
from stepwise_mol_bio.future.pcr import Pcr, pcr
from freezerbox import load_db

db = load_db()

p = stepwise.Protocol()

# Strains:
# - I've been using sz224 and sz228 as my "default" strains, so I can just 
#   continue that here.
p += stepwise.load('grow -p "b1h/qpcr" sz224 sz228')
p += stepwise.load('zymo_quick_rna')

# RT:

gene_specific_primers = dict(
        primer='s2+s5',
        gene_specific_primer=True
)
rt_params = {
        'SuperScript IV': dict(
            preset='thermo/superscript-iv/ezdnase',
            **gene_specific_primers,
        ),
        'SuperScript IV VILO': dict(
            preset='thermo/superscript-iv/vilo',
        ),
        'ProtoScript II': dict(
            preset='neb/protoscript-ii',
            **gene_specific_primers,
        ),
        'AMV': dict(
            preset='neb/amv',
        ),
}
rt_templates = ['sz224', 'sz228']
rt_samples = []

for rt in rt_params:
    for template in rt_templates:
        rt_sample = ReverseTranscribe(
                template,
                **rt_params[rt],
                no_rt_control=True,
                db=db,
        )
        rt_sample.name = rt
        rt_samples.append(rt_sample)

p += reverse_transcribe(rt_samples)

# PCR:

pcr_samples = []
pcr_primers = ('s1','s2'), ('s3','s5')

for rt_sample in rt_samples:
    for rt_control in ('+RT', '−RT'):
        for fwd, rev in pcr_primers:
            pcr_sample = Pcr.from_tags(
                    f'{rt_sample.template}, {rt_sample.name}, {rt_control}',
                    fwd, rev,
                    preset='luna',
                    product_length_bp=74,

                    # Temporary hack b/c I don't have enough master mix.
                    reaction_volume_uL=18,

                    db=db,
            )
            pcr_samples.append(pcr_sample)

p += pcr(3 * pcr_samples)

p.print()
