#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from dbp import densiometry
from dbp.densiometry import BandFitParams
from scipy.stats import linregress
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from dataclasses import dataclass
from itertools import count
from more_itertools import flatten
from math import inf
from pathlib import Path
from joblib import Memory
from copy import copy
from numbers import Real

memory = Memory('.cache', verbose=0)

class LadderBand:

    def __init__(self, len_bp, mass_ng_per_ug, center_guess, **kwargs):
        self.len_bp = len_bp
        self.mass_ng_per_ug = mass_ng_per_ug
        self.ignore = kwargs.pop('ignore', False)

        # Force every band to fit close to where it's supposed to.
        kwargs.setdefault('center_min', center_guess - 45)
        kwargs.setdefault('center_max', center_guess + 45)
        self.p0 = BandFitParams(center_guess, **kwargs)

class Lane:

    def fit(self, plot):
        self.x, self.y = densiometry.curve_from_image(plot)
        self.ps = cached_fit_bands(self.x, self.y, self.p0)

class LadderLane(Lane):

    def __init__(self, i, bands, mass_ug, center_offset=0):
        self.i = i
        self.bands = bands
        self.mass_ug = mass_ug
        self.p0 = [copy(x.p0) for x in self.bands]

        for p0 in self.p0:
            p0._center_guess += center_offset

class UnknownLane(Lane):

    def __init__(self, i, label, bands, load_uL):
        self.i = i
        self.label = label
        self.bands = bands
        self.load_uL = load_uL
        self.p0 = self.bands

    def predict_conc(self, plot, std):
        self.fit(plot)

        self.pxs = np.array([densiometry.gaussian_area(*p) for p in self.ps])
        self.ngs = std.predict(self.pxs)

        self.px = self.pxs.sum()
        self.ng = self.ngs.sum()
        self.ng_uL = self.ng / self.load_uL

neb_1kb_plus = [
        # Ignore these bands because they are too closely spaced to be 
        # reliable.
        LadderBand(10000,  40, 314, ignore=True),
        LadderBand( 8000,  40, 331, ignore=True),
        LadderBand( 6000,  48, 355, ignore=True),
        LadderBand( 5000,  40, 373, ignore=True),
        LadderBand( 4000,  32, 399, ignore=True),

        LadderBand( 3000, 120, 439),
        LadderBand( 2000,  40, 509),
        LadderBand( 1500,  57, 579),
        LadderBand( 1200,  45, 641),

        LadderBand( 1000, 122, 686),
        LadderBand(  900,  34, 718),
        LadderBand(  800,  31, 748),
        LadderBand(  700,  27, 785),
        LadderBand(  600,  23, 824),

        LadderBand(  500, 124, 871),
        LadderBand(  400,  49, 944),
        # Off the gel, but partially visible in some lanes.
        LadderBand(  400,  49, 967, ignore=True),
]
p2 = [
        BandFitParams(237, center_min=200, center_max=250),
        BandFitParams(290, center_min=260, center_max=350),
        BandFitParams(430, center_min=400, center_max=470),
]
ladders = [
        LadderLane(2, neb_1kb_plus, 1, ),
        LadderLane(3, neb_1kb_plus, 1/3, center_offset=-10),
        LadderLane(4, neb_1kb_plus, 1/9, center_offset=-20),
]
unknowns = [
        UnknownLane(5, 'K', p2, 1),
        UnknownLane(6, 'K', p2, 1/3),
        UnknownLane(7, 'K', p2, 1/9),
        UnknownLane(8, 'T', p2, 1),
        UnknownLane(9, 'T', p2, 1/3),
        UnknownLane(10, 'T', p2, 1/9),
]

img = densiometry.load_image('20220917_compare_comp_cells_v2_plots.tif')
divs = densiometry.find_divisions(img)
plots = dict(densiometry.iter_plots(count(1), img, divs))

linear_range_cutoff_px = 8000

def calc_standard_curve(plots, ladders):
    std_px = []
    std_ng = []

    for ladder in ladders:
        ladder.fit(plots[ladder.i])

        for band, p in zip(ladder.bands, ladder.ps):
            if band.ignore:
                continue

            px = densiometry.gaussian_area(*p)
            ng = ladder.mass_ug * band.mass_ng_per_ug

            if px > linear_range_cutoff_px:
                continue

            std_px.append(px)
            std_ng.append(ng)

    std_px = np.reshape(std_px, (-1, 1))
    std_ng = np.reshape(std_ng, (-1, 1))

    std = LinearRegression(fit_intercept=False, positive=True)
    std.fit(std_px, std_ng)

    r2 = r2_score(std_ng, std.predict(std_px))

    #res = linregress(std_px, std_ng)
    return StdCurve(std, r2), std_px, std_ng

@memory.cache
def cached_fit_bands(x, y, p0):
    return densiometry.fit_bands(x, y, p0)

def plot_lane(lane, **kwargs):
    x, y, ps = lane.x, lane.y, lane.ps

    plt.plot(x, y)
    plt.plot(x, densiometry.gaussian_sum(x, *flatten(ps)), **kwargs)

    heights = [p[0] for p in ps]
    centers = [p[1] for p in ps]
    plt.plot(centers, heights, '+')

class StdCurve:

    def __init__(self, model, r2):
        self.model = model
        self.r2 = r2

    def predict(self, xs):
        scalar = isinstance(xs, Real)
        xs = np.asarray(xs).reshape(-1, 1)
        ys = self.model.predict(xs)
        return ys.item() if scalar else ys


# Fits/regressions

std, px, ng = calc_standard_curve(plots, ladders)

unk_px = []
unk_ng = []
results = []

for unk in unknowns:
    unk.predict_conc(plots[unk.i], std)
    unk_px.extend(unk.pxs)
    unk_ng.extend(unk.ngs)

    row = {
            'Label': unk.label,
            'Load (µL)': unk.load_uL,
            'Signal (px)': unk.px,
            'Mass (ng)': unk.ng,
            'Conc (ng/µL)': unk.ng_uL,
    }
    results.append(row)

# Plotting

w, h = 6, 2 + 1 + 2
plt.figure(figsize=(10, 10))

# Standard curve

plt.subplot(h, w, (1,12))
plt.axline(
        (x := 0, std.predict(x)),
        (x := linear_range_cutoff_px, std.predict(x)),
        linestyle='--',
        label=f'R²={std.r2:.4f}',
)
plt.plot(px, ng, '+', label="Ladder bands")
plt.plot(unk_px, unk_ng, 'o', mfc='none', label="Unknown bands")
plt.axvline(linear_range_cutoff_px, linestyle=':')
plt.xlabel('px')
plt.ylabel('ng')
plt.legend(loc='best')

# Ladder fits

for i, ladder in enumerate(ladders):
    plt.subplot(h, w, (j := 13 + 2*i, j + 1))
    plt.title(f'ladder ({1000 * ladder.mass_ug:.1f} ng)')
    plot_lane(ladder)
    plt.ylim(0, 600)

# Unknown fits

for i, unk in enumerate(unknowns):
    plt.subplot(h, w, (j := 19 + 2*i, j + 1))
    plt.title(f'{unk.label} ({unk.load_uL:.2f} µL)')
    plot_lane(unk, label=f'{unk.px:.1f} px\n{unk.ng:.1f} ng')
    plt.ylim(0, 600)
    plt.legend(loc='best', fontsize='small')

plt.tight_layout()

df = pd.DataFrame(results)

py_path = Path(__file__)
svg_path = py_path.parent / f'{py_path.stem}_fits.svg'
csv_path = py_path.with_suffix('.csv')

plt.savefig(svg_path)
df.to_csv(csv_path, index=False)

plt.show()
