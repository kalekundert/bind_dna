#!/usr/bin/env python3

# The 3.0, 2.0, and 1.5 bands in the ladder are well-resolved and comparable in 
# intensity to my samples of interest.  I know how much mass is in these bands, 
# so I can use them to make a standard curve that relates intensity to the 
# quantity of DNA loaded on the gel.  In turn, I can use this to compare the 
# results from the gel to the results from the Nanodrop.

import numpy as np
import polars as pl
import matplotlib.pyplot as plt
from scipy.stats import linregress
from color_me import ucsf

xlsx_path = '20220901_compare_purification_f186.xlsx'

df_std = pl.read_excel(xlsx_path, 2)
df_unk = pl.read_excel(xlsx_path, 1)

# I loaded 25 ng of ladder (10 µL at 2.5 ng/µL):
load_mass_ug = 25 / 1000

df_std = df_std.select([
    load_mass_ug * pl.col("Band mass (ng per 1 µg load)").alias('mass_ng'),
    pl.col('Intensity (px)').alias('intensity_px'),
])
df_unk = df_unk.select([
    pl.col('Construct').alias('construct'),
    pl.col('Purification').alias('purification'),
    pl.col('Volume (µL)').alias('volume_uL'),
    pl.col('Intensity (px)').alias('intensity_px'),
])

def fit_std_curve(df_std):
    mass_ng = df_std.select('mass_ng').to_numpy().flatten()
    intensity_px = df_std.select('intensity_px').to_numpy().flatten()

    fit = linregress(intensity_px, mass_ng)

    def predict(x):
        return x * fit.slope + fit.intercept

    return fit, predict

fit, predict = fit_std_curve(df_std)

df_unk = df_unk.with_columns([
    pl.col('intensity_px').map(predict).alias('mass_ng'),
])

debug(df_unk)

x_fit = np.linspace(0, 5000)

plt.plot(
        x_fit,
        predict(x_fit),
        label=f'R²={fit.rvalue**2:.6f}',
        linestyle='--',
        color=ucsf.dark_grey[0],
)
plt.plot(
        df_std['intensity_px'].view(),
        df_std['mass_ng'].view(),
        label='standards',
        linestyle='none',
        marker='+',
        markeredgecolor=ucsf.blue[0],
        markerfacecolor='none',
)
plt.plot(
        df_unk['intensity_px'].view(),
        df_unk['mass_ng'].view(),
        label='unknowns',
        linestyle='none',
        marker='o',
        markeredgecolor=ucsf.blue[0],
        markerfacecolor='none',
)
plt.legend(loc='best')
plt.xlabel('Intensity (px)')
plt.ylabel('Mass (ng)')
plt.savefig('fit_load_mass.svg')
plt.show()

