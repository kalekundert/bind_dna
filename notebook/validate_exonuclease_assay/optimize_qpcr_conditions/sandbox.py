#!/usr/bin/env python3

from bind_dna import qpcr
from matplotlib import pyplot as plt

import numpy as np

pcr = qpcr.load_pcr_data('20190611_check_efficiency.toml')
melt = qpcr.load_melt_data('20190611_check_efficiency.toml')
cq = qpcr.load_cq_data('20190611_check_efficiency.toml')

def plot_pcr(ax, df, style=lambda k: {}):
    colors 
    for (template_pg, well), df in df.groupby(['template_pg', 'well']):
        ax.semilogy(df['cycle'], df['rfu'])

def plot_melt(ax, df, style=lambda k: {}):
    for key, g in df.groupby('well'):
        ax.plot(g['temperature'], g['drfu_dtemp'], **style(key))

    x = df['temperature']
    ax.set_xlim(min(x), max(x))

def plot_efficiency(ax, df):
    from scipy.stats import linregress

    x = df['template_pg']
    y = df['cq']
    x_log = np.log10(x)

    ax.plot(x, y, '+')
    ax.set_xscale('log')

    m, b, r, p, err = linregress(x_log, y)
    x_fit = np.linspace(*ax.get_xlim())
    y_fit = np.polyval((m, b), np.log10(x_fit))

    fits = {
            'm': m,
            'b': b,
            'r2': r**2,
            'factor': 10**(-1/m),
            'efficiency': 100 * (10**(-1/m) - 1),
    }

    pprint(fits)

    ax.plot(x_fit, y_fit)
    ax.set_ylabel('Cq')
    ax.set_xlabel('DNA (pg)')
    ax.set_ylim(0, 40)
    ax.grid(True)


    



fig, ax = plt.subplots(1, 3)

plot_pcr(ax[0], pcr)
plot_melt(ax[1], melt)
plot_efficiency(ax[2], cq)

plt.show()

#print(pcr)
#print(melt)
#print(cq)

