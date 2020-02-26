#!/usr/bin/env python3

"""\
Calculate how to dilute genes amplified by PCR.

Usage:
    dilute <tsv> [options]
    dilute <tag> <ng_uL> [options]

Arguments:
    <tsv>
        A TSV file exported by the nanodrop listing the concentration of each 
        species (in ng/μL).  The "Sample Name" column must contains "tags" 
        referring to either plasmids (p01), fragments (f01), or oligos (o01).

    <tag>
        An individual tag, as described above.

    <ng_uL>
        An individual ng/µL measurement.

Options:
    -c --conc <nM>
        The concentration achieve after dilution.  The default is read from the 
        "Conc" column of the database containing the given species, or 75 nM if 
        no such concentration is specified.

        Note: 1 pmol/µL is 1000 nM.
        Note: 75 nM corresponds to 125 ng/μL of the DHFR control template.  

    -v --volume <uL>  [default: 20]
        The volume of diluted amplicon to make.

    -V --stock-volume <uL>
        The volume of concentrated DNA that you want to use for the dilution.  
"""

import docopt
import pandas as pd
import bind_dna as dbp
import stepwise

eval_or_none = lambda x: x if x is None else eval(x)

args = docopt.docopt(__doc__)
final_uL = eval(args['--volume'])
stock_uL = eval_or_none(args['--stock-volume'])

def get_type(tag):
    try:
        return dbp.get_cols(tag)['Type']
    except:
        return 'DNA'

def get_final_conc(tag):
    if args['--conc']:
        return eval(args['--conc']) 
    else:
        return dbp.get_conc(tag)

if tsv := args['<tsv>']:
    nanodrop = pd.read_csv(tsv, sep='\t')
    df = pd.DataFrame()
    df['amplicon'] = nanodrop['Sample Name']
    df['stock_ng_uL'] = nanodrop['Nucleic Acid(ng/uL)']
    
else:
    df = pd.DataFrame([{
        'amplicon':     args['<tag>'],
        'stock_ng_uL':  float(args['<ng_uL>']),
    }])

df['type'] = df['amplicon'].apply(get_type)
df['mw_da'] = df['amplicon'].apply(dbp.get_mw)
df['stock_nM'] = 1e6 * df['stock_ng_uL'] / df['mw_da']
df['final_nM'] = df['amplicon'].apply(get_final_conc)

if stock_uL is not None:
    df['stock_uL'] = stock_uL
    df['water_uL'] = stock_uL * (df['stock_nM'] / df['final_nM'] - 1)
else:
    df['stock_uL'] = final_uL * df['final_nM'] / df['stock_nM']
    df['water_uL'] = final_uL - df['stock_uL']

pd.set_option('display.precision', 2)

types = '/'.join(df['type'])
if len(x := set(df['final_nM'].dropna())) == 1:
    final_nM = f' to {x.pop():.2f} nM'
else:
    final_nm = ''

df = df[[
        'amplicon',
        'mw_da',
        'stock_ng_uL',
        'stock_nM',
        'final_nM',
        'stock_uL',
        'water_uL',
]]

protocol = stepwise.Protocol()
protocol += f"""\
Dilute the purified {types}{final_nM}:

{df}
"""
print(protocol)
