#!/usr/bin/env python3

"""\
Usage:
    plot_kinetic <toml>
"""

import wellmap
import matplotlib.pyplot as plt
from dbp.plate_reader import BiotekExperiment
from more_itertools import one, unique_everseen as uniq
from pathlib import Path

def load_fam_fluorescence(p):
    expt = BiotekExperiment(p)
    return expt.kinetic['450,521']

def plot_kinetic(df, x, y, *, x_label='', y_label='', labels=[]):
    from matplotlib.offsetbox import AnchoredText

    num_rows = df.row_i.max() - df.row_i.min() + 1
    num_cols = df.col_j.max() - df.col_j.min() + 1

    fig, axes = plt.subplots(
            num_rows, num_cols,
            sharex=True,
            sharey=True,
            figsize=(12,12),
            tight_layout=True,
    )

    for (i, j), g in df.groupby(['row_i', 'col_j']):
        ax = axes[i,j]

        ax.plot(g[x], g[y])

        label_values = {k: one(uniq(g[k])) for k in labels}

        title = '\n'.join(
                v.format(label_values[k])
                for k, v in labels.items()
        )
        title_box = AnchoredText(
                title,
                frameon=True,
                loc='upper left',
        )

        ax.add_artist(title_box)

    for ax in axes[:, 0]:
        ax.set_ylabel(y_label)
    for ax in axes[-1, :]:
        ax.set_xlabel(x_label)
                
    return fig, axes

if __name__ == '__main__':
    import docopt
    args = docopt.docopt(__doc__)
    toml_path = Path(args['<toml>'])

    df = wellmap.load(
            toml_path=toml_path,
            data_loader=load_fam_fluorescence,
            merge_cols=True,
            path_guess='{0.stem}.xlsx',
    )

    fig, axes = plot_kinetic(
            df, x='minutes', y='read',
            x_label='time [min]',
            y_label='RFU',
            labels={
                'title': '{}',
                'rnase_U_mL': 'RNase H: {} U/mL',
            },
    )

    plt.savefig(toml_path.with_suffix('.svg'))
    plt.show()
