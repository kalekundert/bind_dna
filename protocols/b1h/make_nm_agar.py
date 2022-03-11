#!/usr/bin/env python3

"""
Prepare permissive and selective plates for B1H assays.

Usage:
    make_nm_agar.py <n> [-P] [-S] [-s <concs>] [-a <antibiotics>]

Options:
    -P --no-permissive
        Don't make permissive plates (e.g. +His, +Ura).

    -S --no-selective
        Don't make selective plates (e.g. −His, −Ura, +3AT).

    -s --3at-concs <mM,...>             [default: 10]
        What concentrations of 3-AT to include in the selective plates.  
        Specify multiple concentrations separated by commas.

    -a --antibiotics <names,...>        [default: carb]
        What antibiotics to include in the plates.  Must be either 'carb', 
        'kan', or 'carb,kan'.
"""

import docopt
import stepwise

from stepwise import pl, ul, table
from stepwise_mol_bio import round_up_to_1_sig_fig
from fractions import Fraction

args = docopt.docopt(__doc__)

concs = [float(x) for x in args['--3at-concs'].split(',')]
antibiotics = set(args['--antibiotics'].split(','))
n_conds = (
        (not args['--no-permissive']) +
        (not args['--no-selective']) * len(concs)
)
n_plates = int(args['<n>'])
n_extra_plates = 0.6
plate_volume_mL = 25

p = stepwise.Protocol()

def nm_agar_recipe(n_plates, n_conds):
    # This function is a bit hacky because master mixes don't support multiple 
    # units.  See #57.

    volume_mL = plate_volume_mL * (n_plates + n_extra_plates) * n_conds * 1.2
    volume_mL = round(volume_mL, -1)
    k = volume_mL / 500

    pre_autoclave = table(
            header=['Reagent', 'Final', 'Quantity'],
            rows=[
                    ['M9 salts', '1x', f'{5.25*k:.2f} g'],
                    ['−His/−Ura DO supplement [2]', '', f'{375*k:.2f} mg'],
                    ['water', '', f'to {437*k:.2f} mL']
            ],
            align='<<>>',
    )

    antibiotic_rows = []
    if 'kan' in antibiotics:
        antibiotic_rows += ['kanamycin', '50 mg/mL', '25 µg/mL, 0.5x', f'{250*k:.2f} µL'],
    if 'carb' in antibiotics:
        antibiotic_rows += ['carbenicillin', '100 mg/mL', '100 µg/mL, 1x', f'{500*k:.2f} µL'],

    post_autoclave = table(
            header=['Reagent', 'Stock', 'Final', 'Volume'],
            rows=[
                    ['M9-agar', '', '', f'{437*k:.2f} mL'],
                    ['glucose', '200 mg/mL', '4 mg/mL', f'{10*k:.2f} mL'],
                    ['thiamine', '200 mg/mL', '4 mg/mL', f'{10*k:.2f} mL'],
                    ['MgSO₄', '100 mM', '1 mM', f'{5*k:.2f} mL'],
                    ['ZnSO₄', '100 mM', '10 µM', f'{50*k:.2f} µL'],
                    ['CaCl₂', '1 M', '100 µM', f'{50*k:.2f} µL'],
                    *antibiotic_rows,
                    ['IPTG', '1 M', '10 µM', f'{5*k:.2f} µL'],
            ],
            align='<>>>',
    )

    p = stepwise.Protocol()

    p += pl(
            f"Prepare ≈{458*k:.0f} mL 12/11x NM-agar [1]:",
            pre_autoclave,
            ul(
                f"Add {7.5*k:.2g} g agar",
                "Autoclave at 121°C for 15 min.",
                "Cool to ≈55°C [3].",
            ),
            post_autoclave,
    )

    p.footnotes[1] = "Noyes et al. (2008) DOI: 10.1093/nar/gkn048"
    p.footnotes[2] = "This supplement isn't exactly what [Noyes2008] calls for: it's missing 9 of the amino acids."
    p.footnotes[3] = pl(
            "3-AT is heat labile and will be destroyed if added to media hotter that 55°C.  Store plates containing 3-AT at 4°C for up to 2 months.",
            "https://tinyurl.com/3xz8euwt",
    )

    return p

def permissive_recipe(n_plates):
    permissive = stepwise.MasterMix()
    permissive.volume = 225, 'mL'
    permissive.show_concs = True
    permissive['NM-agar'].stock_conc = Fraction(12,11), 'x'
    permissive['NM-agar'].volume = 206.25, 'mL'
    permissive['histidine'].stock_conc = 100, 'mM'
    permissive['histidine'].volume = 13.5, 'mL'
    permissive['uracil'].stock_conc = 20, 'mM'
    permissive['uracil'].volume = 2.25, 'mL'
    permissive['water'].solvent = True
    permissive.hold_ratios.volume = \
            plate_volume_mL * (n_plates + n_extra_plates), 'mL'
    return permissive

def selective_recipe(n_plates, _3at_conc):
    selective = stepwise.MasterMix()
    selective.volume = 225, 'mL'
    selective.show_concs = True
    selective['NM-agar'].stock_conc = Fraction(12,11), 'x'
    selective['NM-agar'].volume = 206.25, 'mL'
    selective['3-AT'].stock_conc = 1000, 'mM'
    selective['3-AT'].hold_stock_conc.conc = _3at_conc, 'mM'
    selective['water'].solvent = True
    selective.hold_ratios.volume = \
            plate_volume_mL * (n_plates + n_extra_plates), 'mL'
    return selective

p += nm_agar_recipe(n_plates, n_conds)

if not args['--no-permissive']:
    permissive = permissive_recipe(n_plates)
    p += pl(
            f"Prepare {permissive.volume} permissive media:",
            permissive
    )

if not args['--no-selective']:
    for conc in concs:
        selective = selective_recipe(n_plates, conc)
        p += pl(
                f"Prepare {selective.volume} selective media with {conc:g} mM 3-AT:",
                selective
        )

p += "Pour 25 mL plates."

p.print()

