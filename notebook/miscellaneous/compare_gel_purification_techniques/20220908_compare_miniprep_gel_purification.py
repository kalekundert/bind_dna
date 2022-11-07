#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd

from dbp import densiometry
from dbp.densiometry import BandFitParams
from scipy.stats import linregress
from dataclasses import dataclass
from itertools import count
from more_itertools import flatten
from math import inf
from pathlib import Path

def load_image(path):
    img = plt.imread(path)
    return 1 - img / img.max()

class LadderBand:

    def __init__(self, len_bp, mass_ng_per_ug, center_guess, **kwargs):
        self.len_bp = len_bp
        self.mass_ng_per_ug = mass_ng_per_ug
        self.ignore = kwargs.pop('ignore', False)

        # Force every band to fit close to where it's supposed to.
        kwargs.setdefault('center_min', center_guess - 2)
        kwargs.setdefault('center_max', center_guess + 2)
        self.fit_params = BandFitParams(center_guess, **kwargs)

class LadderLane:

    def __init__(self, i, bands, mass_ug):
        self.i = i
        self.bands = bands
        self.mass_ug = mass_ug

neb_1kb_plus = [
        # Ignore these bands because they are too closely spaced to be 
        # reliable.
        LadderBand(10000,  40, 236, ignore=True),
        LadderBand( 8000,  40, 239, ignore=True),
        LadderBand( 6000,  48, 242, ignore=True),
        LadderBand( 5000,  40, 244, ignore=True),
        LadderBand( 4000,  32, 260, ignore=True),

        LadderBand( 3000, 120, 278),
        LadderBand( 2000,  40, 319),
        LadderBand( 1500,  57, 350),
        LadderBand( 1200,  45, 375),

        LadderBand( 1000, 122, 394),
        LadderBand(  900,  34, 406),
        LadderBand(  800,  31, 421),
        LadderBand(  700,  27, 435),
        LadderBand(  600,  23, 453),

        LadderBand(  500, 124, 474),
        LadderBand(  400,  49, 507),
        LadderBand(  300,  37, 539),
        LadderBand(  200,  32, 583),
        # Ignore this band because it is too close to the end of the gel, and 
        # is overwhelmed by background signal.
        LadderBand(  100,  61, 634, ignore=True),
]
ladders = [
        LadderLane(2, neb_1kb_plus, 0.05),
]
unknowns = {
        3: 'p240',
        4: 'p242A',
        5: 'p242B',
        6: 'p2',
        7: 'p236',

        8: 'crude',
        9: 'freeze-and-squeeze',
        10: 'electroelution',
        11: 'QIAEX II',
}

n = 7

img = load_image('20220908_compare_miniprep_gel_purification_densiometry.tif')
divs = densiometry.find_divisions(img)
plots = dict(densiometry.iter_plots(count(1), img, divs))

plt.figure(figsize=(10, 10))

plt.subplot(n, 2, 5)
plt.title('ladder')

def calc_standard_curve(plots, ladders):
    std_px = []
    std_ng = []

    def find_closest_band(center, bands):
        best_dx = inf

        for band in bands:
            if band.ignore:
                continue

            dx = abs(center - band.fit_params.center_guess)
            if dx < best_dx:
                best_band = band
                best_dx = dx

        return best_dx, best_band
            

    for ladder in ladders:
        x, y = densiometry.curve_from_image(plots[2])
        p0 = [x.fit_params for x in ladder.bands]
        ps = densiometry.fit_bands(x, y, p0)

        # Hacky, but whatever.
        plt.plot(x, y)
        plt.plot(x, densiometry.gaussian_sum(x, *flatten(ps)))

        heights = [p[0] for p in ps]
        centers = [p[1] for p in ps]
        plt.plot(centers, heights, '+')

        for band, p in zip(ladder.bands, ps):
            if band.ignore:
                continue

            px = densiometry.gaussian_area(*p)
            ng = ladder.mass_ug * band.mass_ng_per_ug

            std_px.append(px)
            std_ng.append(ng)

    res = linregress(std_px, std_ng)
    return res, std_px, std_ng

std, px, ng = calc_standard_curve(plots, ladders)

plt.subplot(n, 2, (1,4))
plt.plot(px, ng, '+')
plt.axline((0, std.intercept), slope=std.slope, linestyle='--')

results = []

for i, unk in unknowns.items():
    x, y = densiometry.curve_from_image(plots[i])
    bands = [BandFitParams(center_guess=300)]
    p = densiometry.fit_bands(x, y, bands)[0]

    px = densiometry.gaussian_area(*p)
    ng = std.slope * px + std.intercept

    results.append({
        'sample': unk,
        'px': px,
        'ng': ng,
    })

    plt.subplot(n, 2, (1,4))
    plt.plot([px], [ng], label=unk, marker='o', markerfacecolor='none')

    plt.subplot(n, 2, i+3)
    plt.plot(x, y)
    plt.plot(x, densiometry.gaussian(x, *p))
    plt.title(unk)
    plt.ylim(0, 400)


plt.subplot(n, 2, (1,4))
plt.legend(loc='best', fontsize='x-small')
plt.xlabel('px')
plt.ylabel('ng')
plt.xlim(left=0)
plt.ylim(bottom=0)

plt.tight_layout()

df = pd.DataFrame(results)

py_path = Path(__file__)
svg_path = py_path.parent / f'{py_path.stem}_fits.svg'
csv_path = py_path.with_suffix('.csv')

plt.savefig(svg_path)
df.to_csv(csv_path)

plt.show()
