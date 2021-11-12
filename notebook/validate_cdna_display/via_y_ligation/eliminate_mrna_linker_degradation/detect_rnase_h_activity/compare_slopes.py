#!/usr/bin/env python3

"""\
Compare the slopes of fluorescent signals between wells in an RNase H assay.

Usage:
    compare_slopes <toml> [-I] [-y <min,max>]

Arguments:
    <toml>
        A wellmap file describing the layout of the experiment.  The layout 
        should be associated with data from a kinetic run measuring 450/521 
        fluorescence exported from a BioTek plate reader.  The following 
        parameters should be specified:

        - "control": The name of the control, for wells that are just testing 
          that the assay worked.  This fields should be set to '' for all 
          experimental wells.

        - "fit_start_min": Exclude time points earlier than the given value (in 
          minutes) from the linear fit.

        - "fit_stop_min": Exclude time points later than the given value (in 
          minutes) from the linear fit.

        - "style.label": A format string that will be used to create a label 
          for each well.

        - "style.color_by": A list of the field names that should be considered 
          when choosing colors.

        - "style.marker_by": A list of the field names that should be 
          considered when choosing markers.

Options:
    -y --y-lim <min,max>
        The minimum and maximum values to display on the y-axis.

    -I --no-interactive
        Don't show an interactive GUI.
"""

import wellmap
import docopt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
from dbp.plate_reader import BiotekExperiment
from scipy.stats import linregress
from color_me.ucsf import iter_colors, dark_grey
from matplotlib.lines import Line2D
from more_itertools import one, unique_everseen as uniq
from itertools import cycle
from pathlib import Path

class Style:

    @classmethod
    def from_extras(cls, df, extras):
        return cls(df, **extras)


    def __init__(
            self, df,
            label='{well}',
            color_by=['well'],
            marker_by=(),
            sort_by=(),
    ):
        self.label = label
        self.color_by = color_by
        self.marker_by = marker_by
        self.sort_by = sort_by

        def pick_by_group(cols, iter_values):
            if not cols:
                return {}

            keys = df.groupby(cols).groups.keys()
            return {
                    self._normalize_key(k): v
                    for k, v in iter_values(keys)
            }

        self.colors = pick_by_group(
                self.color_by, 
                lambda ks: (
                    (k, color)
                    for color, k in iter_colors(ks)
                )
        )
        self.markers = pick_by_group(
                self.marker_by, 
                lambda ks: (
                    (k, (n, 2, 0))
                    for n, k in enumerate(ks, start=3)
                )
        )

    def pick_color(self, g):
        try:
            key = self._key_from_group(g, self.color_by)
            return self.colors[key]

        except KeyError:
            return dark_grey[0]

    def pick_marker(self, g):
        try:
            key = self._key_from_group(g, self.marker_by)
            return self.markers[key]

        except KeyError:
            return '+'

    def format_label(self, g):
        params = {
                col: ','.join(map(str, uniq(g[col])))
                for col in g
        }

        if not isinstance(self.label, list):
            labels = [self.label]
        else:
            labels = self.labels

        for label in labels:
            try:
                return label.format(**params)
            except KeyError:
                continue

    def _key_from_group(self, g, cols):
        try:
            key = tuple(
                    one(uniq(g[col]))
                    for col in cols
            )
            return self._normalize_key(key)

        except ValueError:
            raise KeyError(cols)

    def _normalize_key(self, k):
        if not isinstance(k, tuple):
            k = (k,)

        return tuple(
                None if isinstance(x, float) and math.isnan(x) else x
                for x in k
        )

def load_plate_reader(p):
    expt = BiotekExperiment(p)
    return expt.kinetic['450,521']

def get_fit_indices(df):
    return (df['minutes'] >= df['fit_start_min']) & \
           (df['minutes'] <= df['fit_stop_min'])

def calc_linear_fits(df):
    fits = {}

    for well, g in df.groupby(['well']):
        i = get_fit_indices(df)
        m, b, r, p, err = linregress(g['minutes'][i], g['read'][i])

        fits[well] = (m, b, r, p, err)

    return fits

def record_fits(df, fits, style, path):
    table = []
    for (well, control), g in df.groupby(['well', 'control']):
        fit = fits[well]
        row = dict(
                well=well,
                label=control or style.format_label(g),
                slope=fit[0],
                intercept=fit[1],
                r2=fit[2]**2,
        )
        table.append(row)

    table = pd.DataFrame(table)
    table.to_csv(path)

def plot_linear_fits(df, fits, style):
    fig, axes = plt.subplots(
            1, 2,
            sharex=True,
            sharey=True,
            figsize=[8.0, 4.8],
            constrained_layout=True,
    )

    axes[0].set_title('Controls')
    axes[1].set_title('Experiments')
    axes[0].set_ylabel('RFU')
    axes[0].set_xlabel('time [min]')
    axes[1].set_xlabel('time [min]')

    artists = {True: [], False: []}

    for control, gi in df.groupby(['control']):
        ax = axes[0 if control else 1]

        for well, gj in gi.groupby(['well']):
            artists[bool(control)] += plot_linear_fit(
                    ax, gj, fits[well],
                    color=style.pick_color(gj),
                    marker=style.pick_marker(gj),
                    label=control or style.format_label(gj),
            )

    blank = Line2D([], [], linestyle='none')
    axes[1].legend(
            handles=artists[True] + [blank] + artists[False],
            bbox_to_anchor=(1.05, 1.00),
            loc='upper left',
            borderaxespad=0,
    )

    return fig, axes

def plot_linear_fit(ax, df, fit, *, color, marker, label):
    i = get_fit_indices(df)
    x, y = df['minutes'], df['read']

    ax.plot(
            x[i], y[i],
            marker=marker,
            linestyle='none',
            color=color[0],
            mfc='none',
    )
    ax.plot(
            x[~i], y[~i],
            marker=marker,
            linestyle='none',
            color=color[2],
            mfc='none',
    )

    m, b, r, p, err = fit
    x_fit = np.linspace(x.min(), x.max())
    y_fit = m * x_fit + b

    ax.plot(
            x_fit, y_fit,
            linestyle=':',
            color=color[0],
    )

    artist = Line2D(
            [], [], 
            marker=marker,
            linestyle=':',
            color=color[0],
            label=label,
    )
    return [artist]

def plot_slopes(df, fits, style, ylim):
    fig, ax = plt.subplots(constrained_layout=True)
    ax.set_ylabel('slope [RFU/min]')

    x = 0
    x_ticks = []
    x_labels = []

    groups = df.groupby(['control', *style.sort_by, 'well'])
    for (control, *_, well), g in sorted(groups):
        m, b, r, p, err = fits[well]
        color = style.pick_color(g)[0]
        label = control or style.format_label(g)

        ax.plot(
                [x,x],
                [0,m],
                color=color,
                linewidth=6,
        )
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

    if ylim:
        ax.set_ylim(*ylim)
    else:
        ax.set_ylim(0, ax.get_ylim()[1])

    return fig, ax


if __name__ == '__main__':
    args = docopt.docopt(__doc__)
    toml_path = Path(args['<toml>'])
    df, extras = wellmap.load(
            toml_path=toml_path,
            data_loader=load_plate_reader,
            merge_cols=True,
            path_guess='{0.stem}.xlsx',
            extras='style',
    )
    if 'fit_start_min' not in df:
        df['fit_start_min'] = 0
    if 'fit_stop_min' not in df:
        df['fit_stop_min'] = max(df['minutes'])

    df['control'] = df['control'].fillna('')
    style = Style.from_extras(df, extras)
    fits = calc_linear_fits(df)

    if args['--y-lim']:
        ylim = map(float, args['--y-lim'].split(','))
    else:
        ylim = None

    record_fits(df, fits, style, toml_path.stem + '_fits.csv')

    plot_linear_fits(df, fits, style)
    plt.savefig(toml_path.stem + '_fits.svg')

    plot_slopes(df, fits, style, ylim)
    plt.savefig(toml_path.stem + '_slopes.svg')

    if not args['--no-interactive']:
        plt.show()

