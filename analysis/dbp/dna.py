#!/usr/bin/env python3

import re
import autoprop
import pandas as pd
import autosnapgene as snap
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqUtils import molecular_weight, MeltingTemp
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
param_pattern = r'(?P<key>\w+)=((?P<value>[^"]\S*)|(?P<value_quoted>".*?"))(\s|$)'

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

def get_cols(tag):
    return dispatch_to_tag(
            tag,
            p=get_plasmid_cols,
            f=get_fragment_cols,
            o=get_oligo_cols,
    ).dropna().to_dict()

def get_conc_str(tag):
    try:
        return get_cols(tag)['Conc']
    except KeyError:
        raise ValueError(f"no concentration specified for {tag!r}")

def get_conc_nM(tag):
    conc_str = get_conc_str(tag)
    return parse_nanomolar(conc_str, get_mw(tag))

def get_conc_ng_uL(tag):
    conc_str = get_conc_str(tag)
    return parse_ng_uL(conc_str, get_mw(tag))

def get_protocol(tag):
    protocol_str = get_cols(tag)['Construction']
    return parse_protocol(protocol_str)


def get_plasmid_cols(tag):
    id = int(str(tag).strip('p'))
    plasmid_db = read_plasmid_db()
    return plasmid_db.loc[id]

def get_plasmid_path(tag):
    id = int(str(tag).strip('p'))
    return work.plasmid_dir / f'{id:>03}.dna'

def get_plasmid_seq(tag):
    dna = snap.parse(get_plasmid_path(tag))
    return DnaSeq(dna.sequence)

def get_fragment_cols(tag):
    id = int(str(tag).strip('f'))
    fragment_db = read_fragment_db()
    return fragment_db.loc[id]

def get_fragment_seq(tag):
    protocol = get_protocol(tag)
    return protocol.product_seq

def get_oligo_cols(tag):
    id = int(str(tag).strip('o'))
    oligo_db = read_oligo_db()
    return oligo_db.loc[id]

def get_oligo_seq(tag):
    raw_seq = get_oligo_cols(tag)['Sequence']
    return seq_from_str(raw_seq)

def get_oligo_tm(tag):
    """
    Return the melting temperature of the given oligo.
    """
    oligo_cols = get_oligo_cols(tag)
    name = oligo_cols['Name']
    seq = oligo_cols['Sequence']

    if m := re.search(r'_TM(\d+)', name):
        return float(m.group(1))

    else:
        seq = seq_from_str(seq)
        return MeltingTemp.Tm_Wallace(seq)

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
    protocol = get_protocol(tag)
    return protocol.method != 'IVT'


def dispatch_to_tag(tag, **kwargs):
    if m := re.match(r'([pfo])\d+', tag):
        type_code = m.group(1)
        return kwargs[type_code](tag)

    else:
        raise ValueError(f"unknown tag '{tag}'")

def parse_protocol(protocol_str):
    method, params_str = protocol_str.split(':', 1)
    params = parse_params(params_str)
    protocol_cls = {
            'PCR': PcrProtocol,
            'IPCR': InversePcrProtocol,
            'IVT': IvtProtocol,
            'RE':  DigestProtocol,
            'GG':  GoldenGateProtocol,
            'IDT': IdtProtocol,
    }
    return protocol_cls.get(method, Protocol)(method, params)

def parse_params(params_str):
    params = {}

    for m in re.finditer(param_pattern, params_str, re.VERBOSE):
        key = m.group('key')
        value = m.group('value') or m.group('value_quoted')[1:-1]

        if key in params:
            raise KeyError("dupplicate key {key!r} in {params_str!r}")

        params[key] = value

    return params

def parse_seconds(time_str):
    time_units = {
            's':        1,
            'sec':      1,
            'second':   1,
            'seconds':  1,
            'm':        60,
            'min':      60,
            'minute':   60,
            'minutes':  60,
            'h':        60*60,
            'hr':       60*60,
            'hour':     60*60,
            'hours':    60*60,
    }
    time_pattern_1 = fr'(?P<time>\d+)\s*(?P<unit>{"|".join(time_units)})'
    time_pattern_2 = fr'(?P<min>\d+)m(?P<sec>\d+)'

    if m := re.fullmatch(time_pattern_1, time_str):
        return int(m.group('time')) * time_units[m.group('unit')]
    if m := re.fullmatch(time_pattern_2, time_str):
        return 60 * int(m.group('min')) + int(m.group('sec'))

    raise ValueError(f"can't interpret {time_str!r} as a time")

def parse_celsius(temp_str):
    temp_pattern = fr'(?P<temp>\d+)\s*°?C'

    if m := re.fullmatch(temp_pattern, temp_str):
        return int(m.group('temp'))

    raise ValueError(f"can't interpret {temp_str!r} as a temperature")

def parse_microliters(vol_str):
    vol_pattern = fr'(?P<vol>\d+)\s*[µu]L'

    if m := re.fullmatch(vol_pattern, vol_str):
        return int(m.group('vol'))

    raise ValueError(f"can't interpret {vol_str!r} as a volume")

def parse_nanomolar(conc_str, mw):
    conc_pattern = r'(?P<conc>\d+)\s?(?P<unit>[nuµ]M|ng/[uµ]L)'
    unit_conversion = {
            'nM': 1,
            'uM': 1e3,
            'µM': 1e3,
            'ng/uL': 1e6 / mw,
            'ng/µL': 1e6 / mw,
    }

    if m := re.match(conc_pattern, conc_str):
        return float(m.group('conc')) * unit_conversion[m.group('unit')]
    else:
        raise ValueError(f"can't interpret {conc_str!r} as a concentration")

def parse_ng_uL(conc_str, mw):
    conc_pattern = r'(?P<conc>\d+)\s?(?P<unit>[nuµ]M|ng/[uµ]L)'
    unit_conversion = {
            'ng/uL': 1,
            'ng/µL': 1,
    }

    if m := re.match(conc_pattern, conc_str):
        return float(m.group('conc')) * unit_conversion[m.group('unit')]
    else:
        raise ValueError(f"can't interpret {conc_str!r} as a concentration")


@lru_cache(maxsize=None)
def read_plasmid_db():
    df = pd.read_excel(work.plasmid_db)
    df = df.set_index(df.index + 2)
    return df

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

def seq_from_str(raw_seq):
    # Ignore nonstandard nucleotides; they'd be too hard to deal with...

    # This regular expression works because `re.sub()` only substitutes the 
    # left-most occurrence of any overlapping patterns.  The non-greedy * is 
    # necessary to avoid eliding everything between the first and last 
    # nonstandard nucleotide.
    seq = re.sub(r'/.*?/', '', raw_seq)

    return DnaSeq(seq)

@autoprop
class Protocol:
    no_default = object()

    def __init__(self, method, params):
        self.method = method
        self.params = params

    def parse_param(self, key, pattern, default=no_default):
        if key not in self.params:
            if default is not Protocol.no_default:
                return default
            else:
                raise KeyError(f"{self.method} protocol missing required {key!r} parameter")

        if m := re.fullmatch(pattern, p := self.params[key]):
            return m.groups()
        else:
            raise ValueError(f"expected {pattern!r}, found {p!r}")

    def get_product_seq(self):
        raise NotImplementedError(f"{self.method}: support for product sequence not yet implemented.")

@autoprop
class PcrProtocol(Protocol):

    def get_template_tag(self):
        return one(self.parse_param('template', r'([fp]?\d+)'))

    def get_template_seq(self):
        return get_seq(self.template_tag)

    def get_primer_tags(self):
        return self.parse_param('primers', r'([o]?\d+),([o]?\d+)')

    def get_primer_seqs(self):
        return [get_oligo_seq(x) for x in self.primer_tags]

    def get_product_seq(self):
        seq = self.template_seq
        primers = p = self.primer_seqs
        primer_pairs = [
                (p[0], p[1].reverse_complement()),
                (p[1], p[0].reverse_complement()),
        ]

        for fwd, rev in primer_pairs:
            # Assume perfect complementarity in the last 15 bases.  This is a 
            # bit of a hack...
            primer_ends = fwd[-15:], rev[:15]
            i, j = sorted(seq.find(x) for x in primer_ends)
            if i > 0 and j > 0:
                break
        else:
            raise ValueError(f"{self.primer_tags[0]!r} and {self.primer_tags[1]!r} not found in {self.template_tag!r}")

        return fwd[:-15] + seq[i:j] + rev

    def get_product_len(self):
        return len(self.product_seq)

    def get_annealing_temp_celsius(self):
        if 'Ta' in self.params:
            return parse_celsius(self.params['Ta'])
        else:
            tms = [get_oligo_tm(x) for x in self.primer_tags]
            return min(tms) + 1

    def get_extension_time_seconds(self):
        if 'tx' in self.params:
            return parse_seconds(self.params['tx'])

        time_sec = 30 * self.product_len / 1000
        if time_sec <= 10: return 10
        if time_sec <= 15: return 15
        return (1 + time_sec // 30) * 30

    def get_scale(self):
        if 'scale' in self.params:
            return parse_microliters(self.params['scale'])

@autoprop
class InversePcrProtocol(PcrProtocol):
    pass

@autoprop
class DigestProtocol(Protocol):

    def get_template_tag(self):
        return one(self.parse_param('template', r'(p?\d+)'))

    def get_template_seq(self):
        return get_plasmid_seq(self.template_tag)

    def get_enzyme_name(self):
        return one(self.parse_param('enzyme', r'([\w\d-]+)'))

    def get_product_seq(self):
        seq = self.template_seq
        enzyme = AllEnzymes.get(self.enzyme_name)
        sites = enzyme.search(seq)
        site = -1 + one(
                sites,
                ValueError(f"{self.enzyme_name!r} does not cut {self.template_tag!r}."),
                ValueError(f"{self.enzyme_name!r} cuts {self.template_tag!r} {len(sites)} times."),
        )
        return seq[site:] + seq[:site]


@autoprop
class IvtProtocol(Protocol):

    def get_template_tag(self):
        return one(self.parse_param('template', r'([fp]?\d+)'))

    def get_template_seq(self):
        return get_seq(self.template_tag)

    def get_product_seq(self):
        t7_promoter = 'TAATACGACTCACTATA'
        seq = str(self.template_seq)
        i = seq.find(t7_promoter) + len(t7_promoter)
        return DnaSeq(seq[i:]).transcribe()

@autoprop
class GoldenGateProtocol(Protocol):

    def get_backbone_tag(self):
        return one(self.parse_param('bb', r'([fp]\d+)'))

    def get_insert_tags(self):
        tag = r'[fp]\d+'
        tags = one(self.parse_param('ins', fr'({tag}(?:,{tag})*)'))
        return tags.split(',')

    def get_enzyme_name(self):
        if 'enzyme' in self.params:
            return one(self.parse_param('enzyme', r'(\w+)'))
        else:
            return 'BsaI'


@autoprop
class IdtProtocol(Protocol):

    def get_product_seq(self):
        return self.params['seq']


