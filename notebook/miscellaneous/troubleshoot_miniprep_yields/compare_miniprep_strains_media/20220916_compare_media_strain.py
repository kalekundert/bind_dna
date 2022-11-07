#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from dbp import densiometry
from dbp.densiometry import BandFitParams
from scipy.stats import linregress
from dataclasses import dataclass
from itertools import count
from more_itertools import flatten
from math import inf
from pathlib import Path
from joblib import Memory

memory = Memory('.cache', verbose=0)

class LadderBand:

    def __init__(self, len_bp, mass_ng_per_ug, center_guess, **kwargs):
        self.len_bp = len_bp
        self.mass_ng_per_ug = mass_ng_per_ug
        self.ignore = kwargs.pop('ignore', False)

        # Force every band to fit close to where it's supposed to.
        kwargs.setdefault('center_min', center_guess - 20)
        kwargs.setdefault('center_max', center_guess + 20)
        self.p0 = BandFitParams(center_guess, **kwargs)

class Lane:

    def fit(self, plot):
        self.x, self.y = densiometry.curve_from_image(plot)
        self.ps = cached_fit_bands(self.x, self.y, self.p0)

class LadderLane(Lane):

    def __init__(self, i, bands, mass_ug):
        self.i = i
        self.bands = bands
        self.mass_ug = mass_ug
        self.p0 = [x.p0 for x in self.bands]

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
        self.ngs = std.intercept + std.slope * self.pxs

        self.px = self.pxs.sum()
        self.ng = self.ngs.sum()
        self.ng_uL = self.ng / self.load_uL

neb_1kb_plus = [
        # Ignore these bands because they are too closely spaced to be 
        # reliable.
        LadderBand(10000,  40, 331, ignore=True),
        LadderBand( 8000,  40, 351, ignore=True),
        LadderBand( 6000,  48, 370, ignore=True),
        LadderBand( 5000,  40, 387, ignore=True),
        LadderBand( 4000,  32, 409, ignore=True),

        LadderBand( 3000, 120, 450),
        LadderBand( 2000,  40, 527),
        LadderBand( 1500,  57, 598),
        LadderBand( 1200,  45, 665),

        LadderBand( 1000, 122, 721),
        LadderBand(  900,  34, 752),
        LadderBand(  800,  31, 791),
        LadderBand(  700,  27, 835),
        LadderBand(  600,  23, 883),

        LadderBand(  500, 124, 934),
        LadderBand(  400,  49, 1002),
        # Ignore this band because it's truncated by the end of the gel.
        LadderBand(  300,  37, 1069, ignore=True),
]
p2 = [
        BandFitParams(500),
        BandFitParams(330),
        BandFitParams(270),
]
p1124 = [
        BandFitParams(140),
        BandFitParams(300),
]
ladders = [
        LadderLane(2, neb_1kb_plus, 1),
        LadderLane(3, neb_1kb_plus, 1/3),
        LadderLane(4, neb_1kb_plus, 1/9),
]
unknowns = [
        UnknownLane(5, 'K2', p2, 1),
        UnknownLane(6, 'K2', p2, 0.1),
        UnknownLane(7, 'K4', p1124, 1),
        UnknownLane(8, 'K4', p1124, 0.1),
        UnknownLane(9, 'T2', p2, 1),
        UnknownLane(10, 'T2', p2, 0.1),
        UnknownLane(11, 'T4', p1124, 1),
        UnknownLane(12, 'T4', p1124, 0.1),
]

img = densiometry.load_image('20220916_compare_media_strains_laser_plots.tif')
divs = densiometry.find_divisions(img)
plots = dict(densiometry.iter_plots(count(1), img, divs))

linear_range_cutoff_px = 10000

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

    res = linregress(std_px, std_ng)
    return res, std_px, std_ng

@memory.cache
def cached_fit_bands(x, y, p0):
    return densiometry.fit_bands(x, y, p0)

def plot_lane(lane):
    x, y, ps = lane.x, lane.y, lane.ps

    plt.plot(x, y)
    plt.plot(x, densiometry.gaussian_sum(x, *flatten(ps)))

    heights = [p[0] for p in ps]
    centers = [p[1] for p in ps]
    plt.plot(centers, heights, '+')

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

w, h = 6, 2 + 1 + 4
plt.figure(figsize=(10, 10))

# Standard curve

plt.subplot(h, w, (1,12))
plt.plot(px, ng, '+')
plt.plot(unk_px, unk_ng, 'o', mfc='none')
plt.axline(
        (0, std.intercept),
        slope=std.slope,
        linestyle='--',
        label=f'R²={std.rvalue**2:.4f}',
)
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
    plt.subplot(h, w, (j := 19 + 3*i, j + 2))
    plt.title(f'{unk.label} ({unk.load_uL} µL)')
    plot_lane(unk)
    plt.ylim(0, 600)

plt.tight_layout()

df = pd.DataFrame(results)

py_path = Path(__file__)
svg_path = py_path.parent / f'{py_path.stem}_fits.svg'
csv_path = py_path.with_suffix('.csv')

plt.savefig(svg_path)
df.to_csv(csv_path, index=False)

plt.show()
