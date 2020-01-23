#!/usr/bin/env python3

import re
import pandas as pd
import autosnapgene as snap
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqUtils import molecular_weight
from Bio.Restriction import AllEnzymes
from pathlib import Path
from exmemo import Workspace
from functools import lru_cache
from more_itertools import one

work = Workspace.from_path(__file__)
work.sequence_dir = work.root_dir / 'sequences'
work.plasmid_dir  = work.sequence_dir / 'plasmids'
work.plasmid_db   = work.sequence_dir / 'plasmids.xlsx'
work.fragment_db  = work.sequence_dir / 'fragments.xlsx'
work.oligo_db     = work.sequence_dir / 'oligos.xlsx'

DnaSeq = lambda x: Seq(x.upper(), generic_dna)

def get_seq(tag):
    return dispatch_to_tag(
            tag,
            p=get_plasmid_seq,
            f=get_fragment_seq,
            o=get_oligo_seq,
    )

def get_mw(tag):
    return molecular_weight(
            seq=get_seq(tag),
            double_stranded=is_double_stranded(tag),
            circular=is_circular(tag),
    )

def get_plasmid_path(tag):
    id = int(str(tag).strip('p'))
    return work.plasmid_dir / f'{id:>03}.dna'

def get_plasmid_seq(tag):
    dna = snap.parse(get_plasmid_path(tag))
    return DnaSeq(dna.sequence)

def get_fragment_seq(tag):
    id = int(str(tag).strip('f'))
    fragment_db = read_fragment_db()
    construction = fragment_db.at[id,'Construction']
    method = construction.split(':')[0]
    seq_from_construction = {
            'PCR': get_fragment_seq_pcr,
            'RE':  get_fragment_seq_digest,
            'IVT': get_fragment_seq_ivt,
    }
    seq_getter = seq_from_construction[method]
    return seq_getter(construction)

def get_fragment_seq_pcr(construction):
    template_tag, = parse_param(r'template=([fp]?\d+)', construction)
    primer_tags = parse_param(r'primers=(\d+),(\d+)', construction)

    seq = get_seq(template_tag)
    primers = p = [get_oligo_seq(x) for x in primer_tags]
    primer_pairs = [
            (p[0], p[1].reverse_complement()),
            (p[1], p[0].reverse_complement()),
    ]

    for fwd, rev in primer_pairs:
        # Assume perfect complementarity in the last 15 bases.  This is a bit 
        # of a hack...
        primer_ends = fwd[-15:], rev[:15]
        i, j = sorted(seq.find(x) for x in primer_ends)
        if i > 0 and j > 0:
            break
    else:
        raise ValueError(f"{primer_ids[0]!r} and {primer_ids[1]!r} not found in {template_id!r}")


    return fwd[:-15] + seq[i:j] + rev

def get_fragment_seq_digest(construction):
    template_tag, = parse_param(r'template=(p?\d+)', construction)
    enzyme_name, = parse_param(r'enzyme=(\w+)', construction)

    seq = get_plasmid_seq(template_tag)
    enzyme = AllEnzymes.get(enzyme_name)
    sites = enzyme.search(seq)
    site = -1 + one(
            sites,
            ValueError(f"{enzyme_name!r} does not cut {template_tag!r}."),
            ValueError(f"{enzyme_name!r} cuts {template_tag!r} {len(sites)} times."),
    )
    return seq[site:] + seq[:site]

def get_fragment_seq_ivt(construction):
    template_tag, = parse_param(r'template=([fp]?\d+)', construction)

    t7_promoter = 'TAATACGACTCACTATA'
    seq = str(get_seq(template_tag))
    i = seq.find(t7_promoter) + len(t7_promoter)
    return DnaSeq(seq[i:]).transcribe()

def get_oligo_seq(tag):
    id = int(str(tag).strip('o'))
    oligo_db = read_oligo_db()

    # Ignore nonstandard nucleotides; they'd be too hard to deal with...

    # This regular expression works because `re.sub()` only substitutes the 
    # left-most occurrence of any overlapping patterns.  The non-greedy * is 
    # necessary to avoid eliding everything between the first and last 
    # nonstandard nucleotide.
    raw_seq = oligo_db.at[id,'Sequence']
    seq = re.sub(r'/.*?/', '', raw_seq)

    return DnaSeq(seq)

def is_circular(tag):
    return dispatch_to_tag(
            tag,
            p=lambda id: True,
            f=lambda id: False,
            o=lambda id: False,
    )

def is_double_stranded(tag):
    return dispatch_to_tag(
            tag,
            p=lambda id: True,
            f=is_fragment_double_stranded,
            o=lambda id: False,
    )

def is_fragment_double_stranded(tag):
    id = int(str(tag).strip('f'))
    fragment_db = read_fragment_db()
    construction = fragment_db.at[id,'Construction']
    method = construction.split(':')[0]
    return method != 'IVT'


def dispatch_to_tag(tag, **kwargs):
    if m := re.match(r'([pfo])(\d+)', tag):
        type_code, id = m.group(1), int(m.group(2))
        return kwargs[type_code](id)

    else:
        raise ValueError(f"unknown tag '{tag}'")

def parse_param(pattern, construction):
    if m := re.search(pattern, construction):
        return m.groups()
    else:
        raise ValueError(f"expected {pattern!r}, found {construction!r}")
@lru_cache(maxsize=None)
def read_fragment_db():
    df = pd.read_excel(work.fragment_db)
    df = df.set_index(
            df['Name'].apply(lambda x: int(x.strip('f'))),
    )
    return df

@lru_cache(maxsize=None)
def read_oligo_db():
    df = pd.read_excel(work.oligo_db)
    df = df.set_index(df.index + 2)
    return df


