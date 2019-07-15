#!/usr/bin/env python3

"""\
Plot data from experiments to determine optimal qPCR parameters for the primers 
I'll be using in my preliminary experiments.

Usage:
    optimize_qpcr.py
"""

import docopt
import bio96
import pandas as pd
import matplotlib.pyplot as plt

from pathlib import Path
from color_me import ucsf

def plot_ta_gradient(toml_path):
    toml_path = Path(toml_path)
    svg_path = toml_path.stem + '.svg'

    df = load_cq(toml_path)
    cq = df.groupby('ta')['cq'].agg(['mean', 'std'])
    cq = cq.reset_index(level='ta')

    fig, ax = plt.subplots()

    ax.plot(df['ta'], df['cq'], '+', color=ucsf.blue[0])
    ax.plot(cq['ta'], cq['mean'], color=ucsf.blue[0])

    ax.set_xlabel("Ta (°C)")
    ax.set_ylabel("Cq")

    fig.tight_layout()
    fig.savefig(svg_path)

def plot_primer_conc(toml_path):
    toml_path = Path(toml_path)
    svg_path = toml_path.stem + '.svg'

    df = load_cq(toml_path)
    i = ['skpp_202_conc', 'qpcr_60_conc']
    cq = df.groupby(i)['cq'].agg(['mean', 'std'])
    cq = cq.reset_index(level=i)
    img = cq.pivot(
            index='skpp_202_conc',
            columns='qpcr_60_conc',
            values='mean',
    )

    fig, ax = plt.subplots()

    mat = ax.matshow(img)
    ax.set_xticklabels([''] + list(img.columns))
    ax.set_yticklabels([''] + list(img.index))
    ax.set_xlabel("[QPCR_60_REV] (nM)")
    ax.set_ylabel("[skpp-202-R] (nM)")

    mat.set_clim(18.6, 19.6)
    bar = fig.colorbar(mat)

    fig.tight_layout()
    fig.savefig(svg_path)

def load_cq(toml_path):

    def parse_csv(toml_path):
        df = pd.read_csv(toml_path / 'Quantification Cq Results.csv')
        df = df.rename({'Well': 'well0', 'Cq': 'cq'}, axis='columns')
        return df[['well0', 'cq']]

    return bio96.load(
            toml_path,
            parse_csv,
            merge_cols={'well0': 'well0'},
            path_guess='{0.stem}/',
    )

        
if __name__ == '__main__':
    args = docopt.docopt(__doc__)

    plot_ta_gradient('20190610_optimize_ta.toml')
    plot_primer_conc('20190610_optimize_primer_conc.toml')

