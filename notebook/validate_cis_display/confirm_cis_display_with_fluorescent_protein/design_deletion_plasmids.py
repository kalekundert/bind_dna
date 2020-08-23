#!/usr/bin/env python3

import sys
import math
import regex
import autoprop
import random
import numpy as np
import pandas as pd
import pytest
import functools
from more_itertools import one
from klab.cloning import mutagenesis
from pathlib import Path

random.seed(0)
np.random.seed(0)
next_id = 140

@autoprop
class Construct:

    def __init__(self, name, *segments):
        self.name = name
        self._segments = list(segments)

    def __repr__(self):
        return f'Construct(name={self.name!r}, segments={[x.name for x in self._segments]!r})'

    def __getitem__(self, key):
        return one(x for x in self._segments if x.key == key)

    @classmethod
    def from_seqs(cls, name, **seqs):
        self = cls(name)
        self._segments = [
                Segment(k, v)
                for k, v in seqs.items()
        ]
        return self

    @classmethod
    def from_prototype(cls, name, prototype, delete=[], replace={}, keep=None):
        self = cls(name)
        self._segments = prototype.segments[:]
        self.delete_segments(*delete)
        self.replace_segments(**replace)
        if keep is not None: self.keep_segments(*keep)
        return self

    def get_seq(self, start=None, end=None):
        return ''.join(x.seq for x in self.get_segments(start, end))

    def get_annotated_seq(self, start=None, end=None):
        return ''.join(x.annotated_seq for x in self.get_segments(start, end))

    def get_segments(self, start=None, end=None):

        class selector:

            def __init__(self):
                    self.found_start = (start == None)
                    self.found_end = False

            def __call__(self, segment):
                if not self.found_start:
                    self.found_start = (segment.key == start)
                    self.found_end = (end and segment.key == end)
                    return self.found_start

                elif not self.found_end:
                    self.found_end = (end and segment.key == end)
                    return True

                else:
                    return False

        return [x for x in filter(selector(), self._segments)]

    def have_matching_segment(self, segment):
        try:
            hit = self[segment.key]
            return hit.seq == segment.seq
        except ValueError:
            return False

    def keep_segments(self, *keepers):
        self._check_segments_exist(keepers)
        self._segments = [
                x for x in self.segments
                if x.key in keepers
        ]

    def delete_segments(self, *deletions):
        self._check_segments_exist(deletions)
        self._segments = [
                x for x in self.segments
                if x.key not in deletions
        ]

    def replace_segments(self, **replacement_seqs):
        self._check_segments_exist(replacement_seqs.keys())
        for i, segment in enumerate(self.segments):
            if segment.key in replacement_seqs:
                self._segments[i] = Segment(
                    segment.key,
                    replacement_seqs[segment.key],
                )

    def _check_segments_exist(self, names):
        unknown = set(names) - set(x.name for x in self.segments)
        if unknown:
            unknown_str = ','.join(repr(x) for x in names if x in unknown)
            raise ValueError(f"segment(s) not found: {unknown_str}")


@autoprop
class Segment:

    def __init__(self, key, seq):
        self._key = key
        self._name = None
        self._seq = seq

    def __repr__(self):
        seq = self.seq if len(self.seq) < 10 else self.seq[:10] + '…'
        if self.name == self.key:
            return f"Segment(key={self.key!r}, seq={seq!r})"
        else:
            return f"Segment(key={self.key!r}, name={self.name!r}, seq={seq!r})"

    def get_key(self):
        return self._key

    def get_name(self):
        return self._name or self.key

    def set_name(self, name):
        self._name = name

    def get_seq(self):
        return self._seq.strip('/')

    def set_seq(self, seq):
        self._seq = seq

    def get_annotated_seq(self):
        return self._seq

# `sbp primers list 3 -rj -v`
SR022 = 'gggttgtctcctctgatagc' # CTCC
SR091 = 'gtactcagagattgccggag' # AGAT
SR151 = 'ctaggggatggtccaatacg' # ATGG

def design_constructs():
    linkers = dict(
            # Flexible `GGGGSAAP` linker from the [Odegrip2004] supplement.  This 
            # is the linker I've been using so far, but it contains a G->C mutation 
            # (relative to the primers that were used to clone it) that converts 
            # the final amino acid from Ala to Pro.
            odegrip2004_pro='gggggaggaggatcagcggcccca',

            # Flexible `GGGGSAAA` linker implied by [Odegrip2004] primer lists, see 
            # above.
            odegrip2004_ala='gggggaggaggatcagcggccgca',

            # Flexible `SGSETPGTSESATPES` ("XTEN") linker from [Guilinger2014].
            guilinger2014='agtggcagcgaaacccctggtacatcggaatctgcgactccggaaagt',

            # Flexible `GSAGSAAGSGEF` linker from [Waldo1999].
            waldo1999='ggatccgctggctccgctgctggttctggcgaattc',

            # Rigid `LAEAAAKEAAAKEAAAKEAAAKAAA` linker from [Arai2001].
            arai2001='ctggcagaagcagcggcgaaagaggcggcggctaaggaagccgcggcaaaggaagcagccgctaaagccgctgca',
    )

    prototype = Construct.from_seqs(
            'FULL',
            fwd=SR022,
            buffer='',
            t7='taatacgactcactataggg',
            utr='cttaagta',
            rbs='taaggaggaaaaaatatg',
            strep='gctagctggagccacccgcagttcgaaaaaggcgcc',
            mwasabi='gtgagcaaaggcgaagaaaccaccatgggcgtgattaaaccggatatgaaaattaaactgaaaatggaaggcaacgtgaacggccatgcgtttgtgattgaaggcgaaggcgaaggcaaaccgtatgatggcaccaacaccattaatctggaagtgaaagaaggcgcgccgctgccgtttagctatgatattctgaccaccgcgtttagctatggcaaccgtgcgtttaccaaatatccggatgatattccgaactattttaaacagagctttccggaaggctatagctgggaacgtaccatgacctttgaagataaaggcattgtgaaagtgaaaagcgatattagcatggaagaagatagctttatttatgaaattcatctgaaaggcgagaactttccgccgaacggcccggtgatgcagaaagaaaccaccggctgggatgcgagcaccgaacgtatgtatgtgcgtgatggcgtgctgaaaggtgatgtgaaaatgaaactgctgctggaaggcggcggccatcatcgtgtggattttaaaaccatttatcgtgcgaaaaaagcggtgaaactgccggattatcattttgtggatcatcgtattgaaattctgaaccatgataaagattataacaaagtgaccgtgtatgaaattgcggtggcgcgtaacagcaccgatggcatggatgaactgtataaa',
            mwasabi_stop='',
            linker=linkers['odegrip2004_pro'],
            repa_n='actgatcttcaccaaacgtattaccggcaggtaaagaacccgaatccggtgttcactccccgtgaaggtgccggaacgctgaagttctgcgaaaaactgatggaaaaggcggtgggcttcacctcccgttttgatttcgccattcatgtggcgcatgcccgttcccgtggtctgcgtcggcgcatgccaccggtgctgcgtcgacgggctattgatgcgctgctgcaggggctgtgtttccactatgacccgctggccaaccgcgtccagtgttccatcaccacactggccattgagtgcggactggcgacagagtccggtgcaggaaaactctccatcacccgtgccacccgggccctgacgttcctgtcagagctgggactgattacctaccagacggaatatgacccgcttatcgggtgctacattccgaccgacatcacgttcacactggctctgtttgctgcccttgatgtgtctgaggatgcagtggcagctgcgcgccgcagtcgtgttgaatgggaaaacaaacagcgcaaaaagcaggggctggatacactgggtatggatgagctgatagcgaaagcttggcgttttgtgcgtgagcgtttccgcagttaccagacagagcttcagtcccgtggaataaaacgtgcccgtgcgcgtcgtgatgcgaacagagaacgtcaggacatcgtcaccctagtgaaacggcagctgacgcgtgaaatctcggaaggacgcttcactgctaatggtgaggcggtaaaacgcgaagtggagcgtcgtgtgaaggag',
            repa_rho='cgcatgattctgtcacgtaaccgcaattacagccggctggccacagcttctccctga',
            cis='aagtgatctcctcagaataatccggcctgcgccggaggcatccgcacgcctgaagcccgccggtgcacaaaaaaacagcgtcgcatgcaaaaaacaatctcatcatccaccttctggagcatccgattccccctgtttttaatacaaaatacgcctcagcgacggggaattt',
            ori='tgcttatccacatttaactgcaagggacttccccataaggttacaaccgttcatgtcataaagcgccagccgccagtcttacagggtgcaatgtatcttttaaacacctgtttatatctcctttaaactacttaattacattcatttaaaaagaaaacctattcactgcctgtcctgtggacagacag',
            rev=rc(SR151),
    )

    constructs = [
            Construct.from_prototype('p108', prototype,
                keep={'fwd', 'buffer', 'rev'},
            ),

            # Test is the DNA is shifted without any protein expression.
            Construct.from_prototype('p109', prototype,
                keep={'fwd', 'buffer', 't7', 'rev'},
            ),
            Construct.from_prototype('p110', prototype,
                keep={'fwd', 'buffer', 't7', 'ori', 'rev'},
            ),
            Construct.from_prototype('p111', prototype,
                keep={'fwd', 'buffer', 't7', 'repa_rho', 'cis', 'rev'},
            ),
            Construct.from_prototype('p112', prototype,
                keep={'fwd', 'buffer', 't7', 'repa_rho', 'cis', 'ori', 'rev'},
            ),

            # Test if the DNA is shifted with GFP expression.
            Construct.from_prototype('p113', prototype,
                delete={'linker', 'repa_n', 'repa_rho', 'cis', 'ori'},
                replace={'mwasabi_stop': 'tga/'},
            ),
            Construct.from_prototype('p114', prototype,
                delete={'linker', 'repa_n', 'repa_rho', 'cis'},
                replace={'mwasabi_stop': 'tga/'},
            ),
            Construct.from_prototype('p115', prototype,
                delete={'linker', 'repa_n', 'ori'},
                replace={'mwasabi_stop': 'tga/'},
            ),
            Construct.from_prototype('p116', prototype,
                delete={'linker', 'repa_n'},
                replace={'mwasabi_stop': 'tga/'},
            ),
            Construct.from_prototype('p117', prototype,
                delete={'linker', 'repa_n', 'ori'},
            ),
            Construct.from_prototype('p118', prototype,
                delete={'linker', 'repa_n'},
            ),

            # Test if the DNA is shifted by repA (not fused to anything).
            Construct.from_prototype('p119', prototype, 
                delete={'strep', 'mwasabi', 'mwasabi_stop', 'linker', 'cis', 'ori'},
            ),
            Construct.from_prototype('p120', prototype, 
                delete={'strep', 'mwasabi', 'mwasabi_stop', 'linker', 'ori'},
            ),
            Construct.from_prototype('p121', prototype, 
                delete={'strep', 'mwasabi', 'mwasabi_stop', 'linker'},
            ),

            # Test if the DNA is shifted by mWasabi-repA.
            Construct.from_prototype('p122', prototype, 
                delete={'cis', 'ori'},
            ),
            Construct.from_prototype('p123', prototype, 
                delete={'ori'},
            ),
            Construct.from_prototype('p124', prototype)
    ]

    from itertools import count
    p = count(125)

    for name, seq in linkers.items():
        if seq == prototype['linker'].seq:
            continue
        linker_constructs = [
                Construct.from_prototype(f'p{next(p)}', prototype, 
                    replace={'linker': seq},
                    delete={'cis', 'ori'},
                ),
                Construct.from_prototype(f'p{next(p)}', prototype, 
                    replace={'linker': seq},
                    delete={'ori'},
                ),
                Construct.from_prototype(f'p{next(p)}', prototype, 
                    replace={'linker': seq},
                ),
        ]

        # Give the linker segment a more meaningful name.
        for construct in linker_constructs:
            construct['linker'].name = name

        constructs += linker_constructs

    return prototype, constructs

def select_constructs(*names):
    return [x for x in constructs if x.name in names]

def design_buffer_seq(constructs, force=False):
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    # Make sure all the constructs are referencing the same buffer.
    buffer = one({x['buffer'] for x in constructs})
    assert buffer.seq == ''

    # Check if the buffer has been designed already.
    cache_path = Path('buffer.fa')
    if cache_path.exists() and not force:
        record = one(SeqIO.parse(cache_path, "fasta"))
        buffer.seq = str(record.seq)
        return int(record.id), buffer

    # Design the buffer sequence.
    construct_lens = [len(x.seq) for x in constructs]
    target_len = round_100(max(construct_lens) + 100)
    buffer_len = round_100(target_len - min(construct_lens) + 100)
    buffer.seq = optimize_buffer_seq(buffer_len)

    # Cache the result.
    cache_record = SeqRecord(
            Seq(buffer.seq),
            id=str(target_len),
            name='',
            description='',
    )
    SeqIO.write([cache_record], cache_path, "fasta")
    return target_len, buffer

def design_buffer_gblock(buffer):
    bsai = 'aGGTCTCa'
    return bsai + SR022[7:] + buffer.seq[:-8] + rc(bsai)

def design_buffer_deletions(constructs, target_len):
    primers = {}
    target_tm = calc_tm(SR022)
    target_q5_tm = 66

    for construct in constructs:
        del_len = len(construct.seq) - target_len
        possible_seqs = [
                construct['buffer'].seq[del_len:del_len+x]
                for x in range(10, 31)
        ]
        seq, tm = pick_primer_with_best_tm(possible_seqs, target_tm)
        q5_tm = int(round((tm - target_tm) + target_q5_tm))
        name = f'0_BUF_{construct.name}_TM{q5_tm}'

        primers[name] = seq

    return mutagenesis.consolidate_duplicate_primers(primers)

def design_segment_deletions(prototype, constructs):
    from functools import partial

    primers = {}
    next_id = 1

    def formatter(mismatch, construct, kind, tm, dir):
        nonlocal next_id

        def get_segment_name(seg):
            try:
                return construct[seg.key].name
            except ValueError:
                return seg.name

        segment = mismatch[{'R': 0, 'F': -1}[dir]]
        tag = get_segment_name(segment)
        name = f"{next_id}_{kind}_{tag}_TM{tm:.0f}_{dir}"
        next_id += 1
        return name

    for construct in constructs:
        if construct is prototype: continue

        for match_5, mismatch, match_3 in iter_mismatches(prototype, construct):
            kwargs = {
                    'start': match_5[0].key,
                    'end': match_3[-1].key,
            }
            pd = mutagenesis.PrimerDesigner()
            pd.min_overlap = 15
            pd.tm = 62
            pd.formatter = partial(formatter, mismatch, construct)
            pd.construct = construct.get_annotated_seq(**kwargs)

            # This is hacky, but we treat the backbone differently depending on 
            # if we're cloning a deletion or a linker.  For the former, we want 
            # all the primers aligned on the segment boundaries to minimize the 
            # number of primers we need.  We achieve this by replacing the 
            # mismatched sequence with 'X'.  For the latter, we want to be 
            # smarter in terms of making the smallest mutation possible.  We 
            # achieve this by using the full backbone sequence, so the primer 
            # designer can pick the shortest primers.

            if [x.name for x in mismatch] == ['linker']:
                pd.backbone = prototype.get_seq(**kwargs)

                # Another hack: the linker is really GC-rich, so we need to 
                # raise the target Tm a bit to get primers for the Odegrip2004 
                # (Ala) linker.
                if construct['linker'].name == 'odegrip2004_ala':
                    pd.tm = 63

            else:
                pd.backbone = \
                        prototype.get_seq(match_5[0].key, match_5[-1].key) + \
                        'X' + \
                        prototype.get_seq(match_3[0].key, match_3[-1].key)


            primers.update(pd.design_primers())

    return mutagenesis.consolidate_duplicate_primers(primers)

def pick_primer_with_best_tm(primers, target_tm):
    primer_tms = [
            (x, calc_tm(x))
            for x in primers
    ]
    primer_tms.sort(key=lambda x: abs(x[1] - target_tm))
    return primer_tms[0]

def iter_mismatches(prototype, construct):
    from itertools import groupby
    from more_itertools import windowed

    if prototype.seq == construct.seq:
        return

    # Convert the groups into lists so that windowing doesn't cause groups to 
    # get lost (as the iterator is advanced).
    groups = [
            (k, list(g))
            for k, g in groupby(
                prototype.segments,
                key=construct.have_matching_segment,
            )
    ]
    for window in windowed(groups, 3, step=2):
        pattern, mismatches = zip(*window)
        assert pattern == (True, False, True)
        yield mismatches

def fix_primer_names(primers, id):
    renamed = {}

    def resolve_dups(name):
        dups = name.split(', ')
        dup_parts = [
                x.split('_')[1:]
                for x in dups
        ]

        if len(dup_parts) == 1:
            return dup_parts[0]

        if dup_parts[0][0] == 'BUF':
            return keep_all_dups(dup_parts)
        else:
            return keep_unique_dups(dup_parts)

    def keep_all_dups(dup_parts):
        from itertools import groupby
        assert set(len(x) for x in dup_parts) == {3}

        parts = []
        for x in zip(*dup_parts):
            parts += [k for k, g in groupby(x)]
        return parts

    def keep_unique_dups(dup_parts):
        parts = []
        for a in dup_parts[0]:
            if all((a in x) for x in dup_parts[1:]):
                parts.append(a)

        if parts[0] not in ('MUT', 'DEL'):
            parts = ['DEL', *parts]

        return parts

    for name, seq in primers.items():
        name = name.upper()
        parts = resolve_dups(name)
        parts = [str(id)] + parts; id += 1
        renamed['_'.join(parts)] = seq

    return renamed

def report_buffer_gblock(buffer):
    from stepwise import tabulate

    dups = find_duplicate_kmers(buffer.seq)
    gblock = design_buffer_gblock(buffer)

    print("Duplicate k-mers:")
    print(dups)
    print()
    print(f">Buffer gBlock (len={len(gblock)} bp)")
    print(gblock)
    print()

def find_duplicate_kmers(seq):
    from itertools import count
    from collections import Counter

    dup_counts = {}

    for k in count(target_kmer(seq)):

        counter = Counter()
        for kmer in iter_kmers(seq, k):
            counter[kmer] += 1

        dup_counts[k] = len([x for x in counter.values() if x > 1])

        if dup_counts[k] == 0:
            break

    df = pd.DataFrame(dup_counts.items(), columns=['kmer', 'num_duplicates'])
    return df

def target_kmer(seq):
    from math import log, ceil
    return ceil(log(len(seq), 4))

@functools.lru_cache
def golden_gate_sites():
    from Bio.Restriction import BsaI, BbsI, BsmBI, BtgZI, SapI

    sites = []
    for enz in [BsaI, BbsI, BsmBI, BtgZI, SapI]:
        site = enz.site.lower()
        sites += [site, rc(site)]

    return sites

def calc_tm(seq):
    from Bio.SeqUtils.MeltingTemp import Tm_NN
    return Tm_NN(seq)

def round_100(x):
    return int(100 * math.ceil(x / 100))

def calc_gc(seq):
    from Bio.SeqUtils import GC
    return GC(seq) / 100

def rc(seq):
    from Bio.Seq import Seq
    from Bio.Alphabet import generic_dna
    return str(Seq(seq).reverse_complement())


def optimize_buffer_seq(n, scores=None):
    from tqdm import tqdm

    seq = ''.join(random.choices('atcg', k=n - len(SR091)))
    seq = mutate_golden_gate_sites(seq)

    # The `mutate_overused_kmers()` move is not included because it has a 
    # negative effect on performance.  I think the problem is that are lots of 
    # overused k-mers, and mutating all of them is just too big of a change.

    move_weights, moves = zip(*[
            (1, mutate_region),
    ])
    score_terms = [
            score_gc,
            score_repeats,
            score_golden_gate,
            score_inertness,
    ]

    def calc_score(seq):
        seq = SR022 + seq + rc(SR091)

        score = 0
        veto = False

        for score_term in score_terms:
            s, v = score_term(seq)
            score += s
            veto = veto or v

        return score, veto

    score, veto = calc_score(seq)
    best_seq, best_score = seq, score
    

    len0 = len(seq)
    N = len0 * 10
    T = 1 / (3 * len0)
    i = 0
    num_moves = 0
    num_accepted = 0

    with tqdm(total=N) as progress:
        while i < N:
            # Change the sequence somehow.
            move = random.choices(moves, weights=move_weights, k=1)[0]
            proposed_seq = move(seq)

            # Decide whether or not to keep the change.
            proposed_score, proposed_veto = calc_score(proposed_seq)
            p_accept = math.exp(-(proposed_score - score) / T) \
                    if proposed_score > score else 1
            is_accepted = (p_accept >= random.random())

            #debug(i, seq, proposed_seq, score, proposed_score, veto, proposed_veto, p_accept, is_accepted)

            if is_accepted:
                seq = proposed_seq
                score = proposed_score
                veto = proposed_veto
                num_accepted += 1
            num_moves += 1

            if scores is not None:
                scores.append(score)

            if veto:
                continue

            i += 1
            progress.update(1)

            if score < best_score:
                best_score = score
                best_seq = seq

            if score == 0:
                break

    print(f"Acceptance rate  (w/ veto): {num_accepted}/{num_moves}={100*num_accepted/num_moves:.2f}%")
    print(f"Acceptance rate (w/o veto): {num_accepted}/{i}={100*num_accepted/i:.2f}%")

    # Make sure we didn't accidentally change the length of the sequence.
    assert len(seq) == len0

    return (best_seq + rc(SR091)).lower()

def mutate_region(seq):
    k_max = target_kmer(seq)
    k = int(random.triangular(1, k_max, 1))
    i = random.randrange(len(seq) - k)
    return seq[:i] + ''.join(random.choices('atcg', k=k)) + seq[i+k:]

def mutate_golden_gate_sites(seq):
    return mutate_sites(seq, golden_gate_sites())

def mutate_overused_kmers(seq):
    from collections import Counter

    k = target_kmer(seq)
    counter = Counter()

    for kmer in iter_kmers(seq, k):
        counter[kmer] += 1

    overused = [
            kmer
            for kmer, counts in counter.items()
            if counts > 1
    ]
    return mutate_sites(seq, overused)

def mutate_sites(seq, sites):
    hits = [
            (i, i + len(site))
            for site in sites
            if (i := seq.find(site)) > 0
    ]
    for i,j in hits:
        shuffled_hit = ''.join(random.choices('atcg', k=j-i))
        seq = seq[:i] + shuffled_hit + seq[j:]

    return seq

def score_gc(seq, window_size=20):
    n_tot = 0
    n_bad = 0

    for kmer in iter_kmers(seq, window_size):
        n_bad += not (0.4 <= calc_gc(kmer) <= 0.6)
        n_tot += 1

    return n_bad / n_tot, False

def score_repeats(seq):
    from collections import Counter

    k = target_kmer(seq)
    counts = Counter()

    for kmer in iter_kmers(seq, k):
        counts[kmer] += 1

    n_bad = sum(x - 1 for x in counts.values() if x > 1)
    n_tot = sum(x for x in counts.values())
    return n_bad / n_tot, False

def score_golden_gate(seq):
    import re

    n_bad = 0
    for site in golden_gate_sites():
        n_bad += len(re.findall(site, seq))

    return n_bad, bool(n_bad)

def score_inertness(seq):
    features = [
            # T7 promoter
            'taatacgactcactata',

            # RBS
            'aggagg',
            'aggagg.....atg',
            'aggagg......atg',
            'aggagg.......atg',
            'aggagg........atg',

            # Rho
            'gccggaggcatccgcacgc',

            # DnaA binding site
            'ttatccaca',

            # repA binding motif
            'ca.ttaa.tg',   # [Giraldo1992]

            # AT 9-mer (conserved motif)
            'tttaaa',       # [Masai1987]

            # Weak terminators
            'tgaagcccgccggtgcacaaaaaaa',
            'ggagcatccgattccccctgttttt',
            'aaaatacgcctcagcgacggggaattt',

            # CIS
            'aagtgatctcctcagaataatccggcctgcgccggaggcatccgcacgcctgaagcccgccggtgcacaaaaaaacagcgtcgcatgcaaaaaacaatctcatcatccaccttctggagcatccgattccccctgtttttaatacaaaatacgcctcagcgacggggaattt',
            
            # oriR
            'tgcttatccacatttaactgcaagggacttccccataaggttacaaccgttcatgtcataaagcgccagccgccagtcttacagggtgcaatgtatcttttaaacacctgtttatatctcctttaaactacttaattacattcatttaaaaagaaaacctattcactgcctgtcctgtggacagacag'
    ]

    n_bad = 0
    for feature in iter_strandedness(features):
        max_errors = len(feature) // 10
        pattern = '(?:%s){s<=%.0f}' % (feature, max_errors)
        if regex.search(pattern, seq):
            n_bad += 1

    return n_bad, bool(n_bad)

def iter_kmers(seq, k):
    from more_itertools import windowed
    for kmer in windowed(seq, k):
        yield ''.join(kmer)

def iter_strandedness(seqs):
    for seq in seqs:
        yield seq
        yield rc(seq)


if __name__ == '__main__':
    prototype, constructs = design_constructs()
    #constructs = select_constructs('p109')
    #debug(constructs)

    target_len, buffer = design_buffer_seq(constructs, '-f' in sys.argv)
    buffer_primers = design_buffer_deletions(constructs, target_len)
    segment_primers = design_segment_deletions(prototype, constructs)
    primers = fix_primer_names({**buffer_primers, **segment_primers}, next_id)

    report_buffer_gblock(buffer)
    mutagenesis.report_primers_to_table(primers)


@pytest.mark.parametrize(
    'start,end,expected', [
            ('a', 'a', 'a'),
            ('a', 'b', 'ab'),
            ('a', 'c', 'abc'),
            ('a', 'd', 'abcd'),
            ('a', 'e', 'abcde'),
            
            ('b', 'b', 'b'),
            ('b', 'c', 'bc'),
            ('b', 'd', 'bcd'),
            ('b', 'e', 'bcde'),
            
            ('c', 'c', 'c'),
            ('c', 'd', 'cd'),
            ('c', 'e', 'cde'),
            
            ('d', 'd', 'd'),
            ('d', 'e', 'de'),

            ('e', 'e', 'e'),
])
def test_construct_get_seq(start, end, expected):
    p = Construct.from_seqs('P', a='a', b='b', c='c', d='d', e='e')
    assert p.seq == 'abcde'
    assert p.get_seq(start, end) == expected

@pytest.mark.parametrize(
    'kwargs,expected', [(
            {'delete': 'b'},
            [('a', 'b', 'cde')],
        ), (
            {'delete': 'bc'},
            [('a', 'bc', 'de')],
        ), (
            {'delete': 'bd'},
            [('a', 'b', 'c'), ('c', 'd', 'e')],
        ),
])
def test_iter_mismatches(kwargs, expected):
    prototype = Construct.from_seqs('P', a='a', b='b', c='c', d='d', e='e')
    construct = Construct.from_prototype('C', prototype, **kwargs)

    actual = [
            tuple(
                ''.join(x.seq for x in segs)
                for segs in mismatch
            )
            for mismatch in iter_mismatches(prototype, construct)
    ]
    assert actual == expected

@pytest.mark.parametrize(
    'seq,window,score', [
        ('aaaa', 2, 3/3),
        ('aaag', 2, 2/3),
        ('aagg', 2, 2/3),
        ('aaga', 2, 1/3),
        ('gaga', 2, 0/3),
])
def test_score_gc(seq, window, score):
    assert score_gc(seq, window) == (pytest.approx(score), False)

@pytest.mark.parametrize(
    'seq,score', [
        ('a', 0),

        ('aa', 1/2),
        ('at', 0/2),

        ('aaa', 2/3),
        ('aat', 1/3),
        ('att', 1/3),
        ('atc', 0/3),

        ('aaaa', 3/4),
        ('aaat', 2/4),
        ('aatt', 2/4),
        ('atcc', 1/4),
        ('atcg', 0/4),
])
def test_score_repeats(seq, score):
    assert score_repeats(seq) == (pytest.approx(score), False)

@pytest.mark.parametrize(
    'seq,score', [
        ('aaaaaa', 0),

        # Make sure SR091 doesn't cause an infinite loop.
        (SR091, 0),
        (rc(SR091), 0),

        # BsaI
        ('ggtctc', 1),
        ('gagacc', 1),

        # BbsI
        ('gaagac', 1),
        ('gtcttc', 1),

        # BsmBI
        ('cgtctc', 1),
        ('gagacg', 1),

        # BtgZI
        ('gcgatg', 1),
        ('catcgc', 1),

        # SapI
        ('gctcttc', 1),
        ('gaagagc', 1),
])
def test_score_golden_gate(seq, score):
    assert score_golden_gate(seq) == (score, bool(score))

@pytest.mark.parametrize(
    'seq,score', [
        ('', 0),

        # Make sure SR091 doesn't cause an infinite loop.
        (SR091, 0),
        (rc(SR091), 0),

        # 10% mismatches, rounded down
        ('taatacgactcactata', 1),
        ('Xaatacgactcactata', 1),
        ('XXatacgactcactata', 0),

        # Reverse complement recognized, too.
        ('tatagtgagtcgtatta', 1),
        ('Xatagtgagtcgtatta', 1),
        ('XXtagtgagtcgtatta', 0),

    ]
)
def test_score_inertness(seq, score):
    assert score_inertness(seq) == (score, bool(score))
