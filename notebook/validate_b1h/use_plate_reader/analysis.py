#!/usr/bin/env python3

import wellmap
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

from kbkbio.plate_reader import BiotekExperiment
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel
from scipy.optimize import minimize_scalar
from numpy import inf
from functools import cached_property

def load_plate_reader(path):
    expt = BiotekExperiment(path)
    return expt.kinetic[600]

df = wellmap.load(
        '20220323_b1h_s4_s5_s16_s22.toml',
        data_loader=load_plate_reader,
        merge_cols=True,
        path_guess='{0.stem}.xlsx',
)
#df = df.query('well == "D3"')

# grid = sns.relplot(
#         data=df,
#         x='time_min',
#         y='read',
#         row='row',
#         col='col',
# )

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
        x_train = time_min.reshape(-1, 1)
        y_train = od600.reshape(-1, 1)

        rbf = RBF(length_scale=1e2, length_scale_bounds=(1e1, 1e16))
        noise = WhiteKernel(noise_level_bounds=(1e-32, 1))

        self.gpr = GaussianProcessRegressor(rbf + noise)
        self.gpr.fit(x_train, y_train)

    def predict_od600(self, time_min):
        x = time_min.reshape(-1, 1)
        return self.gpr.predict(x)

    def predict_rate(self, time_min):
        # https://stats.stackexchange.com/questions/373446/computing-gradients-via-gaussian-process-regression
        l = self.gpr.kernel_.k1.length_scale
        a = self.gpr.alpha_

        x = time_min.reshape(-1, 1)
        x_train = self.gpr.X_train_.reshape(1, -1)
        dx = x_train - x

        f = np.exp(-(dx**2) / (2*l**2))
        df = f * dx / l**2

        return df @ a

    @cached_property
    def max_rate(self):
        x_train = self.gpr.X_train_
        res = minimize_scalar(
                lambda x: -self.predict_rate(x),
                bounds=(min(x_train), max(x_train)),
                method='bounded',
        )
        res.fun *= -1
        return res

def calc_growth_rate(df):
    model = GrowthModel(df['time_min'].values, df['read'].values)
    res = model.max_rate
    return pd.Series({
        'rate': res.fun.item(),
        'rate_x': res.x.item(),
    })


gb = df.groupby(['well', 'row_i', 'col_j'], as_index=False)

gb = gb.apply(calc_growth_rate)

img = gb.pivot('row_i', 'col_j', 'rate')
print(img)

ax = sns.heatmap(img, annot=True, fmt='.2e', cmap='viridis')
#ax = sns.heatmap(img, annot=True, fmt='.0f', cmap='viridis')
#plt.imshow(img)
plt.show()




# fig, ax = plt.subplots(2, 1, sharex=True)

# ax[0].plot(df['time_min'], df['read'])

# ax[1].semilogy(df['time_min'], df['read'])


# df = df.query('read <= 0.20')
# p = fit_growth_model(df['time_min'], df['read'])
# debug(p)

# t_fit = np.linspace(0, 800)

# ax[0].plot(t_fit, growth_model(t_fit, *p))
# ax[1].semilogy(t_fit, growth_model(t_fit, *p))


# grid.ax.set_ylim(0, 1)


# Calculate growth rates:
# - Discard data points above a certain OD.
#   - Justification is that OD becomes non-linear at high values, and cells 
#     stop growing eventually.  The principled cutoff value probably depends on 
#     the instrucment and the strain, but it's probably good enough just to 
#     pick an arbitrary point.  Maybe 10% of final value?  That wouldn't work 
#     for flat traces though.  Maybe just a fixed cutoff, e.g. 0.15?
#
# - Fit exponential curve
#   - two parameters: initial OD and doubling time.
#
# - Fit Gaussian process
#   - get maximal growth rate
#
# Questions:
# - Are fits good?
#   
#   - Superimpose fit over raw data.
#   - Lots of plots, just show them all?  
#
# - Are growth rates constant despite changes in initial OD?
#
#   - This should be the case, and it would be good, because it means I don't 
#     need to dilute carefully.
#
#   - Bar graph: growth rate vs dilution.
#
# - Does the one-plasmid system work as well as the two-plasmid one?
#
#   - 3 variables: system, target, media
#   - Just make bar graph of growth rates?

#plt.show()


