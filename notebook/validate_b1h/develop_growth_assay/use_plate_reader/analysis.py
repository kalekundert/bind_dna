#!/usr/bin/env python3

"""
Usage:
    analysis.py od <path> [<wells>] [-o <path>]
    analysis.py rate <path> [-o <path>]
    analysis.py cmp <path> [-o <path>]

Options:
    -o --output <path>
        Save the generated figure to the given path.  The file type will be 
        inferred from the file extension.
"""

import wellmap
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import warnings

from kbkbio.plate_reader import BiotekExperiment
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel
from sklearn.exceptions import ConvergenceWarning
from scipy.optimize import minimize_scalar
from numpy import inf
from joblib import Memory
from natsort import natsorted, natsort_key
from color_me import ucsf
from docopt import docopt
from functools import cached_property

memory = Memory('cache')

class GrowthModel:
    """
    Fit a growth curve as a Gaussian process.

    This approach allows the maximal growth rate to be calculated, in a way 
    that (i) uses all of the available data, (ii) is not too sensitive to 
    noise, and (iii) can easily accommodate replicates.

    The growth rate can in fact be calculated analytically, but doing so 
    requires manually calculating the derivative of the kernel function, so 
    it's not trivial to change the kernel.  I've chosen to use the RBF kernel, 
    which has the correlation between any two points decay exponentially with 
    respect to the distance between them.  This kernel is notable for being 
    infinitely differentiable and therefore very smooth.  Growth curves are 
    also very smooth, so this makes sense and seems to work well.
    """

    def __init__(self, time_min, od600):
        if isinstance(time_min, pd.Series):
            time_min = time_min.values
        if isinstance(od600, pd.Series):
            od600 = od600.values

        x_train = time_min.reshape(-1, 1)
        y_train = od600.reshape(-1, 1)

        rbf = RBF(length_scale=1e2, length_scale_bounds=(1e1, 1e4))
        noise = WhiteKernel(noise_level_bounds=(1e-32, 1))

        self.gpr = GaussianProcessRegressor(rbf + noise)
        
        with warnings.catch_warnings():
            # Very flat growth curves issue a convergence warning because the 
            # correlation time seems like it should be longer than 1e4.  But 
            # allowing very long (e.g. 1e16) correlation times doesn't actually 
            # improve the fit for these curves, and having such a big search 
            # space can cause poor fits for less flat curves.  Fundamentally, 
            # the 1e1-1e4 range makes sense for this kind of data, so we don't 
            # need to acknowledge warnings about it.
            warnings.filterwarnings(
                    action='ignore',
                    message='The optimal value found for dimension 0 of parameter k1__length_scale is close to the specified upper bound',
                    category=ConvergenceWarning,
            )
            self.gpr.fit(x_train, y_train)

    def predict_od600(self, time_min):
        x = time_min.reshape(-1, 1)
        return self.gpr.predict(x)

    def predict_rate(self, time_min):
        # https://stats.stackexchange.com/questions/373446/computing-gradients-via-gaussian-process-regression
        l = self.gpr.kernel_.k1.length_scale
        a = self.gpr.alpha_

        x = np.asarray(time_min).reshape(-1, 1)
        x_train = self.gpr.X_train_.reshape(1, -1)
        dx = x_train - x

        f = np.exp(-(dx**2) / (2*l**2))
        df = f * dx / l**2

        return df @ a

    def calc_max_rate(self, fast=True):
        # Dual-annealing is the most consistent, but takes 32s compared to 17s 
        # for brute force.  In theory, I like brute force more, though, because 
        # its faster and not stochastic.
        x_train = self.gpr.X_train_.ravel()
        y_train = self.gpr.y_train_.ravel()

        # We're using a local minimizer, so we need to make sure our initial 
        # guess is close to the global minimum.  This approach only works 
        # bewouldn't work if our functtakes advantage of the fact that we know 
        # we have a smooth function, 
        #dy_train = self.predict_rate(x_train)
        #i = np.argmax(dy_train)
        #bounds = (
        #        x_train[max(i - 1, 0)],
        #        x_train[min(i + 1, len(x_train) - 1)],
        #)

        #res = minimize_scalar(
        #        lambda x: -self.predict_rate(x),
        #        #bracket=(x_train[i], x_train[i+1]),
        #        bounds=bounds,
        #        method='bounded',
        #)
        #res.fun *= -1

        from scipy.optimize import brute, dual_annealing, optimize, OptimizeResult

        def objective(x):
            return -self.predict_rate(x)

        bounds = [(x_train.min(), x_train.max())]

        # x_min = brute(
        #         objective,
        #         ranges=bounds,
        #         finish=optimize,
        #         Ns=100,
        # )

        # res = OptimizeResult()
        # res.x = x_min
        # res.fun = self.predict_rate(x_min)
        # return res

        res = dual_annealing(
                 objective,
                 bounds=bounds,
        )

        res.fun *= -1

        # 0.001 is the smallest change that can be recorded by the plate 
        # reader.  Enforcing this as the minimum rate avoids problematic 
        # scenarios like negative rates, or rates smaller than could be 
        # measured.
        min_rate = 0.001 / (x_train.max() - x_train.min())

        if min_rate > res.fun.item():
            res.fun = np.array([[min_rate]])
            res.x = np.array([[-1]])

        return res

def load(path):
    return wellmap.load(
            '20220323_b1h_s4_s5_s16_s22.toml',
            data_loader=load_plate_reader,
            merge_cols=True,
            path_guess='{0.stem}.xlsx',
    )

def load_plate_reader(path):
    expt = BiotekExperiment(path)
    return expt.kinetic[600]

@memory.cache
def calc_max_growth_rates(df):
    cols = [
            'well', 'row_i', 'col_j',
            'system', 'strain', 'target', 'selection', 'dilution',
    ]
    gb = df.groupby(cols, as_index=False)
    return gb.apply(calc_max_growth_rate)

def calc_max_growth_rate(df):
    model = GrowthModel(df['time_min'].values, df['read'].values)
    res = model.calc_max_rate()
    return pd.Series({
        'max_rate': res.fun.item(),
        'max_rate_time_min': res.x.item(),
    })

def plot_growth_curves(df, wells):
    sele = df['well'].isin(wells)
    df = df[sele]

    # Figure out the dimensions of the plot

    i_min, i_max = min(df['row_i']), max(df['row_i'])
    j_min, j_max = min(df['col_j']), max(df['col_j'])

    w = i_max - i_min + 1
    h = j_max - j_min + 1

    fig, axes = plt.subplots(
            w, h,
            squeeze=False,
            sharex=True,
            sharey=True,
    )

    for (well, i, j), g in df.groupby(['well', 'row_i', 'col_j']):
        plot_growth_curve(g, well, ax=axes[i - i_min, j - j_min])

    for i, ax1 in enumerate(axes[:,-1]):
        row = wellmap.row_from_i(i + i_min)
        ax2 = ax1.twinx()
        ax2.set_ylabel(row, rotation='horizontal', ha='left')
        ax2.set_yticks([])

    for j, ax in enumerate(axes[0,:]):
        col = wellmap.col_from_j(j + j_min)
        ax.set_title(col)

    for ax in axes[:,0]:
        ax.set_ylabel('OD600')

    for ax in axes[-1,:]:
        ax.set_xlabel('time (min)')

    return fig

def plot_growth_curve(df, well, *, ax):
    # Plot the raw data:

    df = df.query('well == @well')
    if df.empty:
        raise ValueError(f"no measurements found for well {well}")

    x, y = df['time_min'], df['read']

    ax.plot(x, y, '+', color=ucsf.dark_grey[0])

    # Plot the fitted growth curve:
    
    model = GrowthModel(x, y)

    x_fit = np.linspace(0, max(x), 10 * len(x))
    y_fit = model.predict_od600(x_fit)

    ax.plot(x_fit, y_fit, color=ucsf.blue[0])
    ax.set_xlim(-60, max(x) + 60)

    ax2 = ax.twinx()
    ax2.plot(
            x_fit,
            model.predict_rate(x_fit),
            color=ucsf.red[1],
            linestyle='--',
    )

    res = model.calc_max_rate()

    ax2.axvline(
            res.x.item(),
            color=ucsf.light_grey[0],
            linestyle=':',
    )
    ax2.set_yticks([])

    return ax

def plot_growth_rates(df, *, ax):
    img = gb.pivot('row_i', 'col_j', 'rate')

    ax = sns.heatmap(img, annot=True, fmt='.2e', cmap='viridis')
    #ax = sns.heatmap(img, annot=True, fmt='.0f', cmap='viridis')
    #plt.imshow(img)
    plt.show()

def plot_growth_rates_vs_dilutions(df):
    gb = df.groupby(['system'])
    fig, axes = plt.subplots(
            len(gb), 1,
            sharex=True,
            sharey=True,
            figsize=(6, 4*len(gb)),
    )

    for (system, g), ax in zip(natsorted(gb), axes):
        ax.set_title(system)
        plot_growth_rate_vs_dilution(g, ax=ax)

    for ax in axes:
        ax.set_ylabel('OD600/min')

    axes[-1].set_xlabel('initial OD600')
    fig.tight_layout()

def plot_growth_rate_vs_dilution(df, *, ax):
    colors = {
            'TGG': ucsf.red[0],
            'AAA': ucsf.blue[0],
    }
    line_styles = {
            'His+Ura': '--',
            '10 mM 3-AT': '-',
    }

    for (selection, target), g in df.groupby(['selection', 'target']):
        initial_od = 0.1 * 10**(1-g['dilution'])
        ax.loglog(
                initial_od, g['max_rate'],
                label=f'{target}, {selection}',
                color=colors[target],
                linestyle=line_styles[selection],
        )

        i = (g['max_rate_time_min'] == -1)
        ax.loglog(
                initial_od[i], g['max_rate'][i],
                linestyle='none',
                marker='o',
                markerfacecolor='none',
                markeredgecolor=colors[target],
        )

        ax.set_xlim(min(initial_od), max(initial_od))
        ax.legend(loc='best')

def plot_growth_comparisons(df):
    fig, ax = plt.subplots(figsize=(3, 4))
    xticks = []
    xtick_labels = []

    colors = {
            's4/s5': ucsf.dark_grey[0],
            's16/s22': ucsf.blue[0],
    }

    def iter_groups(df):

        def by_comparisons(item):
            (system, selection), g = item
            return (
                    natsort_key(system),
                    ['His+Ura', '10 mM 3-AT'].index(selection),
            )

        cols = ['system', 'selection']
        groups = sorted(df.groupby(cols), key=by_comparisons)

        for i, (key, g) in enumerate(groups):
            yield i, *key, g

    for i, system, selection, g in iter_groups(df):
        tgg_max_rate = g.query('target == "TGG"')['max_rate'].max()
        aaa_max_rate = g.query('target == "AAA"')['max_rate'].max()
        fold_change = tgg_max_rate / aaa_max_rate

        ax.plot(
                [i, i],
                [0, fold_change],
                color=colors[system],
                linewidth=5,
                solid_capstyle='butt',
        )

        xticks.append(i)
        xtick_labels.append(f'{system}, {selection}')

    ax.set_xticks(xticks, labels=xtick_labels, rotation='vertical')
    ax.set_xlim(-0.5, i+0.5)
    ax.set_ylabel('fold change in max growth rate\n[TGG/AAA]')

    fig.tight_layout()
    return fig

def wells_from_cli(df, wells):
    if not wells:
        return set(df['well'])
    else:
        return {
                wellmap.well_from_ij(i, j)
                for i, j in wellmap.iter_well_indices(args['<wells>'])
        }



if __name__ == '__main__':
    args = docopt(__doc__)
    df = load(args['<path>'])

    if args['od']:
        wells = wells_from_cli(df, args['<wells>'])
        plot_growth_curves(df, wells)

    if args['rate']:
        df = calc_max_growth_rates(df)
        plot_growth_rates_vs_dilutions(df)

    if args['cmp']:
        df = calc_max_growth_rates(df)
        plot_growth_comparisons(df)

    if path := args['--output']:
        plt.savefig(path)
    else:
        plt.show()
