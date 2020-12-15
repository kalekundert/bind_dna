#!/usr/bin/env python3

"""\
Usage:
    compare_slopes <toml>
"""

import wellmap
import docopt
import numpy as np
import matplotlib.pyplot as plt
from dbp.plate_reader import BiotekExperiment
from scipy.stats import linregress
from color_me.ucsf import iter_colors
from more_itertools import one, unique_everseen as uniq
from pathlib import Path

def load_plate_reader(p):
    expt = BiotekExperiment(p)
    return expt.kinetic['450,521']

def calc_linear_fits(df):

    fits = {}

    for well, g in df.groupby(['well']):
        m, b, r, p, err = linregress(g['minutes'], g['read'])

        fits[well] = (m, b, r, p, err)

    return fits

def plot_linear_fits(df, fits):
    fig, axes = plt.subplots(1, 2, sharex=True, sharey=True)

    axes[0].set_title('Controls')
    axes[1].set_title('Experiments')

    for ax in axes[:0]:
        ax.set_ylabel('RFU')
    for ax in axes[-1:]:
        ax.set_xlabel('time [min]')

    for control, gi in df.groupby(['control']):
        ax = axes[0 if control else 1]

        for color, (well, gj) in iter_colors(gi.groupby(['well'])):
            rnase_U_mL = one(uniq(gj['rnase_U_mL']))
            title = one(uniq(gj['title']))
            label = f'{title}' if control else f'{rnase_U_mL} U/mL'

            plot_linear_fit(ax, gj, fits[well], color=color, label=label)

        ax.legend(loc='best')

    return fig, axes

def plot_linear_fit(ax, df, fit, *, color, label):
    x, y = df['minutes'], df['read']

    ax.plot(
            x, y,
            marker='+',
            linestyle='none',
            color=color,
    )

    m, b, r, p, err = fit
    x_fit = np.linspace(x.min(), x.max())
    y_fit = m * x_fit + b

    ax.plot(
            x_fit, y_fit,
            linestyle=':',
            color=color,
            label=label,
    )

def plot_slopes(df, fits):
    fig, ax = plt.subplots()
    ax.set_ylabel('slope [RFU/min]')

    x = 0
    x_ticks = []
    x_labels = []

    groups = df.groupby(['control', 'rnase_U_mL', 'well'])
    for color, ((control, _, well), g) in iter_colors(sorted(groups)):
        m, b, r, p, err = fits[well]

        rnase_U_mL = one(uniq(g['rnase_U_mL']))
        title = one(uniq(g['title']))
        label = f'{title}' if control else f'{rnase_U_mL} U/mL'

        if title == 'FAM':
            continue

        ax.plot(
                [x,x],
                [0,m],
                color=color,
                linewidth=6,
        )
        #ax.plot(
        #        [x,x],
        #        [0,m],
        #        marker='_',
        #        color='black',
        #)
        ax.errorbar(
                [x],
                [m],
                yerr=err,
                color=color,
                capsize=3,
        )
        x_ticks.append(x)
        x_labels.append(label)
        x += 1

    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_labels, rotation='vertical')
    ax.set_xlim(min(x_ticks) - 0.5, max(x_ticks) + 0.5)
    ax.set_ylim(0, ax.get_ylim()[1])

    fig.tight_layout()

    return fig, ax


if __name__ == '__main__':
    args = docopt.docopt(__doc__)
    toml_path = Path(args['<toml>'])
    df = wellmap.load(
            toml_path=toml_path,
            data_loader=load_plate_reader,
            merge_cols=True,
            path_guess='{0.stem}.xlsx',
    )
    fits = calc_linear_fits(df)

    plot_linear_fits(df, fits)
    plt.savefig(toml_path.stem + '_fits.svg')

    plot_slopes(df, fits)
    plt.savefig(toml_path.stem + '_slopes.svg')
    plt.show()
