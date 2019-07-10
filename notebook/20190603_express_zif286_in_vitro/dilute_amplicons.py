#!/usr/bin/env python3

"""\
Calculate how to dilute genes amplified by PCR.

Usage:
    dilute_amplicons.py <ng_uL> [options]

Arguments:
    <ng_uL>
        A TOML file listing the concentration of each species (in ng/μL).

Options:
    -c --conc <nM>  [default: 75]
        The concentration achieve after dilution.  The default is 75 nM, which 
        corresponds to 125 ng/μL of the DHFR control template.

    -v --volume <uL>  [default: 20]
        The volume of diluted amplicon to make.
"""

import docopt
import toml

args = docopt.docopt(__doc__)
ng_uL = toml.load(args['<ng_uL>'])
final_uL = eval(args['--volume'])
final_nM = eval(args['--conc'])

# Molecular weights include primer overhangs (unlike the first time I did this 
# calculation).
mw_da = {
        '11':       981_161,
        '11 - ORI': 760_127,
}
stock_nM = {
        k: 1e6 * ng_uL[k] / mw_da[k]
        for k in ng_uL
}
stock_uL = {
        k: final_uL * final_nM / stock_nM[k]
        for k in ng_uL
}

for k in ng_uL:
    print(f"""\
{k}
  water: {final_uL - stock_uL[k]:.2f} μL
  DNA:   {stock_uL[k]:.2f} μL ({ng_uL[k]:.1f} ng/μL, {stock_nM[k]:.1f} nM)
""")



