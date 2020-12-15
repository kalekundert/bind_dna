#!/usr/bin/env python3

import wellmap
import numpy as np
import matplotlib.pyplot as plt
import color_me
import itertools

from dbp.plate_reader import BiotekExperiment
from scipy.stats import linregress

def load_plate_reader(p):
    expt = BiotekExperiment(p)
    return expt.kinetic['450,521']

df = wellmap.load(
        toml_path='20201211_meaure_o42_o213.toml',
        data_loader=load_plate_reader,
        merge_cols=True,
)

fig, ax = plt.subplots()
x_fit = np.linspace(min(df['minutes']), max(df['minutes']))

iter_colors = lambda it: zip(itertools.cycle(color_me.ucsf.cycle), it)

for color, ((beacon, dnazyme), g) in iter_colors(df.groupby(['beacon', 'dnazyme'])):
    x, y = g['minutes'], g['read']
    m, b, r, p, err = linregress(x, y)
    y_fit = m * x_fit + b

    ax.plot(
            x, y,
            label=f'beacon={beacon}\ndnazyme={dnazyme}\nm={m:.2f} RFU/min\nR²={r**2:.2f}',
            marker='+',
            linestyle='none',
            color=color,
    )
    ax.plot(
            x_fit, y_fit,
            linestyle=':',
            color=color,
    )

ax.legend(bbox_to_anchor=(1, 1), loc='upper left')
ax.set_xlabel('time [min]')
ax.set_ylabel('RFU')

plt.tight_layout()
plt.savefig('20201211_measure_o42_o213.svg')
plt.show()
