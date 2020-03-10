#!/usr/bin/env python3

"""\
Calculate how to dilute genes amplified by PCR.

Usage:
    dilute <tags_or_tsv> [<ng_uL>] [options]

Arguments:
    <tags_or_tsv>
        A comma-separated list of tags referring to either plasmids (p01), 
        fragments (f01), or oligos (o01).

        ~or~

        The path to a TSV file exported by the nanodrop listing the 
        concentration of each species (in ng/μL).  The "Sample Name" column 
        must contains "tags" matching the format described above.

    <ng_uL>
        A comma-separated list of ng/µL measurements.  Not used if a TSV path 
        is specified (in that case concentrations are read from the TSV).  The 
        default is read from the "Conc" column of the database containing the 
        given species.  Note that this is also the default for the target 
        concentration, so you must specify one or the other.

Options:
    -c --conc <nM>
        The concentration achieve after dilution.  The default is read from the 
        "Conc" column of the database containing the given species, or 75 nM if 
        no such concentration is specified.

        Note: 1 pmol/µL is 1000 nM.
        Note: 75 nM corresponds to 125 ng/μL of the DHFR control template.  

    -v --volume <uL>  [default: 20]
        The volume of diluted DNA/RNA to make.

    -V --stock-volume <uL>
        The volume of concentrated DNA/RNA that you want to use for the 
        dilution.
"""

import docopt
import pandas as pd
import bind_dna as dbp
import stepwise
from pathlib import Path

eval_or_none = lambda x: x if x is None else eval(x)

args = docopt.docopt(__doc__)
tags = args['<tags_or_tsv>']
final_uL = eval(args['--volume'])
stock_uL = eval_or_none(args['--stock-volume'])

def get_type(tag):
    try:
        return dbp.get_cols(tag)['Type']
    except:
        return 'DNA'

def get_final_nM(tag):
    if args['--conc']:
        return eval(args['--conc']) 
    else:
        return dbp.get_conc_nM(tag)

if Path(tags).exists():
    nanodrop = pd.read_csv(tags, sep='\t')
    df = pd.DataFrame()
    df['tag'] = nanodrop['Sample Name']
    df['stock_ng_uL'] = nanodrop['Nucleic Acid(ng/uL)']
    
else:
    cli_tags = tags.split(',')
    cli_concs = \
            [float(x) for x in args['<ng_uL>'].split(',')] \
            if args['<ng_uL>'] else \
            [dbp.get_conc_ng_uL(x) for x in cli_tags]

    df = pd.DataFrame({
        'tag': cli_tags,
        'stock_ng_uL': cli_concs,
    })

df['type'] = df['tag'].apply(get_type)
df['mw_da'] = df['tag'].apply(dbp.get_mw)
df['stock_nM'] = 1e6 * df['stock_ng_uL'] / df['mw_da']
df['final_nM'] = df['tag'].apply(get_final_nM)

if stock_uL is not None:
    df['stock_uL'] = stock_uL
    df['water_uL'] = stock_uL * (df['stock_nM'] / df['final_nM'] - 1)
else:
    df['stock_uL'] = final_uL * df['final_nM'] / df['stock_nM']
    df['water_uL'] = final_uL - df['stock_uL']

pd.set_option('display.precision', 2)

types = '/'.join(set(df['type']))
if len(x := set(df['final_nM'].dropna())) == 1:
    final_nM = f' to {x.pop():.2f} nM'
else:
    import sys
    print(x, file=sys.stderr)
    final_nM = ''

protocol = stepwise.Protocol()
protocol += f"""\
Dilute the purified {types}{final_nM} [1]:

{df[['tag','stock_uL','water_uL']]}
"""

protocol.footnotes[1] = f"""\
Concentrations:

{df[['tag', 'mw_da','stock_ng_uL','stock_nM','final_nM']]}
"""

print(protocol)
