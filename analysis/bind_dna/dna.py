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


def get_mw(tag):
    tag = str(tag)

    # Plasmid
    if m := re.match(r'p(\d+)', tag):
        seq = get_plasmid_seq(tag)
        double_stranded = True

    # Fragment
    elif m := re.match(r'f(\d+)', tag):
        seq = get_fragment_seq(tag)
        double_stranded = True

    # Oligo
    elif m := re.match(r'o(\d+)', tag):
        seq = get_oligo_seq(tag)
        double_stranded = False

    else:
        raise ValueError(f"unknown tag '{tag}'")

    return molecular_weight(seq, double_stranded=double_stranded)

def get_plasmid_path(id):
    id = int(str(id).strip('p'))
    return work.plasmid_dir / f'{id:>03}.dna'

def get_plasmid_seq(id):
    dna = snap.parse(get_plasmid_path(id))
    return DnaSeq(dna.sequence)

def get_fragment_seq(id):
    id = int(str(id).strip('f'))
    fragment_db = read_fragment_db()
    construction = fragment_db.at[id,'Construction']
    method = construction.split(':')[0]
    seq_from_construction = {
            'PCR': get_fragment_seq_pcr,
            'RE':  get_fragment_seq_digest,
    }
    seq_getter = seq_from_construction[method]
    return seq_getter(construction)

def get_fragment_seq_pcr(construction):
    template_id, = parse_param(r'template=(p?\d+)', construction)
    primer_ids = parse_param(r'primers=(\d+),(\d+)', construction)

    seq = get_plasmid_seq(template_id)
    primers = p = [get_oligo_seq(x) for x in primer_ids]
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
    template_id, = parse_param(r'template=(p?\d+)', construction)
    enzyme_name, = parse_param(r'enzyme=(\w+)', construction)

    seq = get_plasmid_seq(template_id)
    enzyme = AllEnzymes.get(enzyme_name)
    sites = enzyme.search(seq)
    site = -1 + one(
            sites,
            ValueError(f"{enzyme_name!r} does not cut {template_id!r}."),
            ValueError(f"{enzyme_name!r} cuts {template_id!r} {len(sites)} times."),
    )
    return seq[site:] + seq[:site]

def get_oligo_seq(id):
    id = int(str(id).strip('o'))
    oligo_db = read_oligo_db()

    # Ignore nonstandard nucleotides; they'd be too hard to deal with...

    # This regular expression works because `re.sub()` only substitutes the 
    # left-most occurrence of any overlapping patterns.  The non-greedy * is 
    # necessary to avoid eliding everything between the first and last 
    # nonstandard nucleotide.
    raw_seq = oligo_db.at[id,'Sequence']
    seq = re.sub(r'/.*?/', '', raw_seq)

    return DnaSeq(seq)


def parse_param(pattern, construction):
    if m := re.search(pattern, construction):
        return m.groups()
    else:
        raise ValueError(f"expected {pattern!r}, found {construction!r}")
