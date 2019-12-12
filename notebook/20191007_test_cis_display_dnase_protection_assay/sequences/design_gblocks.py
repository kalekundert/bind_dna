#!/usr/bin/env python3

import autosnapgene as snap
import pandas as pd
from pathlib import Path
from klab.cloning import reverse_translate
import random

random.seed(0)

class Block:
    sequence = None

    def __init__(self, sequence, **feature_kwargs):
        self.sequence = sequence
        self.feature_kwargs = feature_kwargs

    def get_sequence(self, key):
        return try_key(self.sequence, key).upper()

    def get_features(self, key):
        if not self.feature_kwargs:
            return []

        kwargs = {
                k: try_key(v, key)
                for k, v in self.feature_kwargs.items()
                if k != 'qualifiers'
        }
        seq = self.get_sequence(key)
        feature = snap.Feature.from_segment(**kwargs)
        return [(seq, feature)]

class RbsBlock(Block):

    def get_features(self, key):
        rbs_kwargs = {
                k.split('_', 1)[1]: v
                for k, v in self.feature_kwargs.items()
                if k.startswith('rbs_')
        }
        sd_kwargs = {
                k.split('_', 1)[1]: v
                for k, v in self.feature_kwargs.items()
                if k.startswith('sd_')
        }
        return [
                (self.sequence, snap.Feature.from_segment(**rbs_kwargs)),
                ('TAAGGAGGA',   snap.Feature.from_segment(**sd_kwargs)),
        ]

def assemble_gblock(key):
    return ''.join([
        x.get_sequence(key)
        for x in blocks
    ])

def make_gblock_xlsx(keys, path):
    rows = []
    for key in keys:
        rows.append({
            'name': names[key],
            'sequence': assemble_gblock(key),
        })

    df = pd.DataFrame(rows)
    df.to_excel(path)

def make_gblock_snapgene(key, path):
    dna = snap.SnapGene()
    dna.sequence = assemble_gblock(key)

    for block in blocks:
        for seq, feat in block.get_features(key):
            dna.add_feature(feat, seq)

    dna.write(path)

def random_barcode(n=40):
    """
    Generate a random barcode with perfectly even usage of A/C/T/G.
    """
    import itertools as it
    nucs = list(it.islice(it.cycle('GACT'), n))

    while True:
        barcode = ''.join(random.sample(nucs, n))
        if not has_restriction_site(barcode):
            return barcode

def has_restriction_site(seq):
    from Bio.Seq import Seq
    from Bio.Restriction import RestrictionBatch

    mix = RestrictionBatch(restriction_sites)
    hits = mix.search(Seq(seq))

    return any(hits.values())

def try_key(maybe_dict, key):
    return maybe_dict[key] if isinstance(maybe_dict, dict) else maybe_dict

def paragraph(p):
    from textwrap import dedent, fill
    p = dedent(p)
    return p.replace(' \n', ' ')

restriction_sites = 'EcoRV,NruI,BstNI,NdeI,BamHI,BbsI,BsaI,BsmBI,BtgZI,SapI,BbvCI,XmnI'.split(',')

names = {
        'wt': '059_insert',
        'null': '063_insert',
        'redv': '065_insert',
        'rgpd': '067_insert',
        'lrhn': '069_insert',
}
targets = {
        'wt': 'TGG',
        'null': 'AAA',
        'redv': 'GCG',
        'rgpd': 'GCG',
        'lrhn': 'TAT',
}
colors = {
        'yellow': '#ffff99',
        'orange': '#ffcc99',
        'blue': '#99ccff',
}
blocks = [
        # 5' BsaI site for Golden Gate cloning.
        Block('tGGTCTCa'),

        # EcoRV site to create starting point for exonuclease digestion.
        Block(
            name="DNase start",
            sequence='TCCAgtGATATC',
            color=colors['yellow'],
            qualifiers=dict(
                note=paragraph("""\
                    According to NEB, the EcoRV site needs to be at least 5bp 
                    from the end of the DNA. To be a little safe, I use 6bp 
                    here. (Note that T7 RNAP often adds a few G's to the 
                    beginning of transcripts, too, but this cassette won't 
                    necessarily be used with T7 RNAP.) I chose the 6bp sequence 
                    ACTAGC. This begins with one of the validated Golden Gate 
                    sites from Potapov2018 (ACTA) and has relatively even GC 
                    content.

                    https://www.neb.com/tools-and-resources/usage-guidelines/cleavage-close-to-the-end-of-dna-fragments
                """),
            ),
        ),

        # `sr08` primer binding site (forward strand).
        Block('GTATGTCGGCTCTCGTATCG'),

        # Zif268 target sequence.
        Block(
            name={
                k: f"{v.upper()} target"
                for k, v in targets.items()
            },
            sequence={
                k: 'taaattGCG{NNN}GCGattata'.format(NNN=v)
                for k, v in targets.items()
            },
            color=colors['yellow'],
        ),

        # `sr13` primer binding site (forward strand).
        Block('GCAGCGTTTTAGCCTACAAG'),

        # Random barcode.  I won't need it for the qPCR assay (unless I want to 
        # use probes), but including it means less will need to change when I 
        # switch to NGS.
        Block(
            name={
                k: f"barcode {i+1}"
                for i, k in enumerate(names)
            },
            sequence={
                    k: random_barcode()
                    for k in names
            },
            color=colors['yellow'],
        ),

        # `sr09` primer-binding site (reverse strand).
        Block('GTGTATAATCCGCTCCCGAA'),

        # RBS from NEB
        RbsBlock(
            sequence='CTTAAGTATAAGGAGGAAAAAAT',

            rbs_name="RBS",
            rbs_type='RBS',
            rbs_color=colors['orange'],

            sd_name="SD",
            sd_color=colors['orange'],
            sd_qualifiers=dict(
                note=paragraph("""\
                NEB recommends the following 5’-UTR genes being amplified by 
                PCR: gcgaatTAATACGACTCACTATAgggcttaagtaTAAGGAGGaaaaaat (the T7 
                promoter and SD sequence are upper case). Note that this SD has 
                8 consecutive bp of complementarity with the 16S rRNA, compared 
                with just 6 for g10-L (although NEB doesn’t consider the 3’ G 
                to be part of the SD sequence—this seems to be a mistake).

                I don’t know where the cttaagta 5’ of the SD sequence came 
                from. It’s different than the sequence in the DHFR control 
                plasmid, and an internet search doesn’t produce any obvious 
                origins. It’s also notably longer than the 5 bp minimum 
                suggested by NEB above, which is interesting because this 
                sequence is presumably optimized to be short (as it’s already 
                on the long-side for a PCR primer).

                It’s also worth noting that the aaaaaat 3’ of the SD sequence 
                creates a 6 bp spacing, 1 bp longer than ideal. Again, it’s not 
                clear what the rationale for this is."""),
            ),
        ),

        Block(
            name={
                k: " ".join(["Zif268"] + ([] if k == 'wt' else [k.upper()]))
                for k in names
            },
            sequence={
                k: reverse_translate(
                    snap.parse(f'zif268_{k}.prot').protein_sequence,
                    forbidden_seqs=restriction_sites,
                    include_stop=False,
                )
                for k in names
            },
            type='CDS',
            color=colors['blue'],
            directionality='forward',
            is_translated=True,
        ),

        # 3' BsaI site for Golden Gate cloning.
        Block('tGAGACCa'),
]

if __name__ == '__main__':
    make_gblock_xlsx(names, 'gblocks.xlsx')
    for key in names:
        make_gblock(key, f'{names[key]}.dna')
