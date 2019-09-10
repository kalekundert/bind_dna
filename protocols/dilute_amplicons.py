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
    -c --conc <nM>  [default: 75]
        The concentration achieve after dilution.  The default is 75 nM, which 
        corresponds to 125 ng/μL of the DHFR control template.

    -v --volume <uL>  [default: 20]
        The volume of diluted amplicon to make.
"""

import docopt
import pandas as pd
import toml
from pathlib import Path

args = docopt.docopt(__doc__)
ng_uL_path = Path(args['<ng_uL>'])
final_uL = eval(args['--volume'])
final_nM = eval(args['--conc'])

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

# Molecular weights include primer overhangs (unlike the first time I did this 
# calculation).
mw_da = {
        'Zif':                  217898,
        '11':                   981161,
        '11 - ORI':             760127,

        '18':                   42058.08,
        '18 + Cy5':             42058.08 + 739,

        '23':                   1341517.95,
        '23 - ORI':             1341517.95 - 222428.32,
        '24':                   1306315.78,
        '24 - ORI':             1306315.78 - 222428.32,
        '25':                   1339657.88,
        '25 - ORI':             1339657.88 - 222428.32,
        '26':                   1332262.68,
        '26 - ORI':             1332262.68 - 222428.32,
        '27':                   1275392.36,
        '27 + Cy5':             1275392.36 + 739,
        '38 + Cy5':             1314324.47 + 739,
        '38 - repA + Cy5':       550553.03 + 739,
        '39 + Cy5':             1267973.64 + 739,
        '40 + Cy5':             1308759.94 + 739,
        '41 + Cy5':             1280950.03 + 739,
        '41 - repA + Cy5':       517178.59 + 739,
        '42 + Cy5':             1277239.69 + 739,

        'pKBK034': 1.33e6,
        'pKBK035': 1.78e6,
        'pKBK027': 2.38e6,
        'pKBK038': 2.42e6,
        'pKBK039': 2.37e6,
        'pKBK040': 2.41e6,
}

try:
    df['mw_da'] = df['amplicon'].apply(lambda x: mw_da[str(x)])
except KeyError as err:
    print(f"Unknown amplicon: {err}")
    raise SystemExit

df['stock_nM'] = 1e6 * df['stock_ng_uL'] / df['mw_da']
df['stock_uL'] = final_uL * final_nM / df['stock_nM']
df['water_uL'] = final_uL - df['stock_uL']

pd.set_option('display.precision', 2)
print(f"target volume: {final_uL:.2f} uL")
print(f"target conc:   {final_nM:.2f} nM")
print()
print(df)
