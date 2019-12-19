#!/usr/bin/env python3

"""\
Calculate how to dilute genes amplified by PCR.

Usage:
    dilute_amplicons.py <ng_uL> [options]

Arguments:
    <ng_uL>
        A file listing the concentration of each species (in ng/μL).  This can 
        either be a hand-written TOML file, or a TSV file exported by the 
        nanodrop.  In the former case, the "Sample Name" column must correspond 
        to the construct names known by this script (e.g. 11, 11 - ORI, etc.).

Options:
    -c --conc <µM>  [default: 5]
        The concentration achieve after dilution.  The default is 1 µM (1 
        pmol/µL).

    -v --volume <uL>  [default: 20]
        The volume of diluted amplicon to make.

    -V --stock-volume <uL>
        The volume of concentrated DNA that you want to use for the dilution.  
"""

import docopt
import pandas as pd
import toml
import re
import autosnapgene as snap
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqUtils import molecular_weight
from pathlib import Path

args = docopt.docopt(__doc__)
ng_uL_path = Path(args['<ng_uL>'])
final_uL = eval(args['--volume'])
final_uM = eval(args['--conc'])

if ng_uL_path.suffix == '.toml':
    df = pd.DataFrame(
            toml.load(args['<ng_uL>']).items(),
            columns=['amplicon', 'stock_ng_uL'],
    )

elif ng_uL_path.suffix == '.tsv':
    nanodrop = pd.read_csv(ng_uL_path, sep='\t')
    df = pd.DataFrame()
    df['amplicon'] = nanodrop['Sample Name']
    df['stock_ng_uL'] = nanodrop['Nucleic Acid(ng/uL)']
    
else:
    print(f"Expected '*.toml' or '*.tsv', not: '{ng_uL_path}'")
    raise SystemExit

PROJECT_DIR = Path(__file__).parents[1]
PLASMID_DIR = PROJECT_DIR / 'sequences' / 'plasmids'

DnaSeq = lambda x: Seq(x.upper(), generic_dna)

def get_mw_kda(name):
    name = str(name)

    utr = DnaSeq('GGGCTTAAGTATAAGGAGGAAAAAAT')
    y_tag_xmni = DnaSeq('GGCTCAAGGGCGGGGGGCGGCGGGGAAAA')

    dna = snap.parse(PLASMID_DIR / f'{int(name):03d}.dna')
    seq = DnaSeq(dna.sequence.upper())

    i = seq.find(utr)
    j = seq.find(y_tag_xmni) + len(y_tag_xmni)

    rna = seq[i:j].transcribe()
    return molecular_weight(rna) / 1000

try:
    df['mw_kda'] = df['amplicon'].apply(get_mw_kda)
except KeyError as err:
    print(f"Unknown amplicon: {err}")
    raise SystemExit

df['stock_uM'] = df['stock_ng_uL'] / df['mw_kda']
df['stock_uL'] = final_uL * final_uM / df['stock_uM']
df['water_uL'] = final_uL - df['stock_uL']

pd.set_option('display.precision', 2)
print(f"target volume: {final_uL:.2f} µL")
print(f"target conc:   {final_uM:.2f} µM (pmol/µL)")
print()
print(df)
