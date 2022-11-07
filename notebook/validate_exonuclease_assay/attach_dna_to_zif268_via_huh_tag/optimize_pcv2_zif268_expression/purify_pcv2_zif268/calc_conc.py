#!/usr/bin/env python3

"""
Calculate protein concentration via SDS PAGE/gel densiometry.

Usage:
    calc_conc.py <conf> [-O | -G]

Options:
    -O --svg-only
        Output the plots to files, but do not launch the GUI.

    -G --gui-only
        Launch the GUI, but do not save the plot to files.
"""

import docopt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import nestedtext as nt
import pint

from dbp import densiometry
from scipy.stats import linregress
from matplotlib.colors import Normalize
from pathlib import Path
from statistics import mean, stdev
from color_me import ucsf
from more_itertools import one, flatten, all_equal
from pydantic import BaseModel, validator, root_validator
from typing import List, Dict, Tuple, Union, Any, Optional
from numpy.typing import ArrayLike
from dataclasses import dataclass

# I can't get pint to work with pandas.  Even with `pint_pandas` imported, I 
# still get `UnitStrippedWarning` whenever I do anything.

UNIT = pint.UnitRegistry()
UNIT.default_format = '~P'
UNIT.setup_matplotlib()
pint.set_application_registry(UNIT)

import warnings
warnings.filterwarnings('error', 'unit')
warnings.simplefilter('error', pint.UnitStrippedWarning)

class LaneConfig(BaseModel):

    @property
    def is_standard(self):
        return False

    @property
    def is_unknown(self):
        return False

class StandardLaneConfig(LaneConfig):
    band: str
    conc: Any

    @validator('conc', pre=True)
    def parse_conc(cls, value):
        return UNIT(value)

    @property
    def bands(self):
        return [self.band]
    @property
    def is_standard(self):
        return True


class UnknownLaneConfig(LaneConfig):
    bands: List[str]
    dilution: float

    @validator('bands', pre=True)
    def parse_bands(cls, value):
        import csv
        return next(csv.reader([value]))

    @property
    def is_unknown(self):
        return True

class AnonymousLaneConfig(LaneConfig):

    @root_validator(pre=True)
    def no_bands(cls, values):
        if 'conc' in values:
            raise ValueError(f"standard expected")
        if 'dilution' in values:
            raise ValueError(f"unknown expected")

class BandConfig(BaseModel):
    center: int
    range: Tuple[float, float]

    @validator('range', pre=True)
    def parse_range(cls, value):
        try:
            lo, hi = value.split('-')
        except:
            raise ValueError(f"cannot interpret {value!r} as a range, expected: '<low>-<high>'") from None

        return float(lo), float(hi)

class Config(BaseModel):
    image: Path
    lanes: List[Union[AnonymousLaneConfig, UnknownLaneConfig, StandardLaneConfig]]
    bands: Dict[str, BandConfig]



@dataclass
class Band:
    name: str
    conf: BandConfig
    p: Tuple[float, float, float, float]

    def y_fit(self, x):
        return densiometry.gaussian(x, *self.p)

    @property
    def area_fit(self):
        return densiometry.gaussian_area(*self.p)

@dataclass
class Lane:
    conf: LaneConfig
    bands: List[Band]
    x: ArrayLike
    y: ArrayLike

    @property
    def y_fit(self):
        ps = flatten(b.p for b in self.bands)
        return densiometry.gaussian_sum(self.x, *ps)

@dataclass
class StdCurve:
    m: pint.Quantity
    b: float
    r: float
    p: float
    m_err: float
    b_err: float
    min_area: float
    max_area: float
    conc_unit: pint.Unit

    @property
    def r2(self):
        return self.r**2

    def calc_conc(self, area):
        return (area - self.b) / self.m

    def calc_area(self, conc):
        return conc * self.m + self.b

def load_image(path):
    img = plt.imread(path)
    return 1 - img / img.max()

def fit_imagej_plots(conf, img):
    divs = densiometry.find_divisions(img)
    iter_plots = densiometry.iter_plots(conf.lanes, img, divs)

    lanes = []
    standards = []
    unknowns = {}

    for i, (lane_conf, plot) in enumerate(iter_plots):
        if not lane_conf.is_standard and not lane_conf.is_unknown:
            continue

        # Figure out what bands we are expecting.
        band_confs = [conf.bands[b] for b in lane_conf.bands]
        band_params = [
                densiometry.BandFitParams(
                    center_guess=b.center,
                    center_min=b.range[0] if b.range else None,
                    center_max=b.range[1] if b.range else None,
                )
                for b in band_confs
        ]

        # Fit each band.
        x, y = densiometry.curve_from_image(plot)
        ps = densiometry.fit_bands(x, y, band_params)

        lane = Lane(
                x=x,
                y=y,
                conf=lane_conf,
                bands=[
                    Band(
                        name=b,
                        conf=conf.bands[b],
                        p=p,
                    )
                    for b, p in zip(lane_conf.bands, ps)
                ],
        )
        lanes.append(lane)

    return lanes

def find_standards(lanes):
    rows = []
    conc_unit = min(
            lane.conf.conc.units
            for lane in lanes
            if lane.conf.is_standard
    )

    for lane in lanes:
        if not lane.conf.is_standard:
            continue

        # We already check that exactly one band is specified for each standard 
        # (when parsing the user input file).  We rely on that fact here, but 
        # we don't need to worry about making a nice error message.
        band = one(lane.bands)

        row = dict(
                conc=lane.conf.conc.m_as(conc_unit),
                area=band.area_fit,
                name=band.name,
        )
        rows.append(row)

    df = pd.DataFrame(rows)

    df.conc_unit = conc_unit

    if not all_equal(df['name']):
        raise ValueError("multiple standards not currently supported")

    return df

def find_unknowns(lanes):
    rows = []

    for lane in lanes:
        if not lane.conf.is_unknown:
            continue

        for band in lane.bands:
            row = dict(
                    dilution=lane.conf.dilution,
                    area=band.area_fit,
                    name=band.name,
            )
            rows.append(row)

    return pd.DataFrame(rows)

def fit_std_curve(standards):
    result = linregress(standards['conc'], standards['area'])
    conc_unit = standards.conc_unit

    return StdCurve(
            m=result.slope,
            b=result.intercept,
            r=result.rvalue,
            p=result.pvalue,
            m_err=result.stderr,
            b_err=result.intercept_stderr,
            min_area=min(standards['area']),
            max_area=max(standards['area']),
            conc_unit=conc_unit,
    )

def calc_unknown_concs(unknowns, std_curve):
    unknowns['raw_conc'] = std_curve.calc_conc(unknowns['area'])
    unknowns['corrected_conc'] = unknowns['raw_conc'] * unknowns['dilution']
    unknowns['within_std_curve'] = (
            (std_curve.min_area < unknowns['area']) &
            (unknowns['area'] < std_curve.max_area)
    )
    unknowns.conc_unit = std_curve.conc_unit

    i = unknowns['within_std_curve']
    concs = unknowns[i]\
            .groupby('name')\
            .agg({'corrected_conc': [np.mean, np.std]})\
            .rename(columns={'corrected_conc': 'conc'})

    concs.conc_unit = std_curve.conc_unit
    return concs

def plot_fits(lanes):
    n_rows = len(lanes)
    n_cols = 1

    fig, axes = plt.subplots(
            n_rows,
            n_cols,
            squeeze=True,
            figsize=(3 * n_cols, 1 * n_rows),
    )

    for ax, lane in zip(axes, lanes):
        ax.plot(lane.x, lane.y, **raw_data_style(lane))
        ax.plot(lane.x, lane.y_fit, **fit_style(lane))
        ax.legend(loc='best', fontsize='xx-small')

    fig.tight_layout()
    return fig

def plot_results(standards, unknowns, std_curve, results):
    fig, axes = plt.subplots(1, 2, squeeze=True)

    plot_std_curve(axes[0], standards, std_curve)
    plot_concs(axes[1], unknowns, results)

    fig.tight_layout()

    return fig

def plot_std_curve(ax, standards, std_curve):
    # Plot the underlying data
    ax.plot(
            standards['conc'],
            standards['area'],
            marker='+',
            color=ucsf.dark_grey[0],
            linestyle='none',
    )

    # Plot the linear regression
    x_fit = np.linspace(*ax.get_xlim())
    y_fit = std_curve.calc_area(x_fit)
    ax.plot(
            x_fit, y_fit,
            linestyle='solid',
            color=ucsf.blue[0],
            label=f'R²={std_curve.r2:g}',
    )
    ax.legend(loc='best')

    # Decorate axes
    ax.set_title('Standard Curve')
    ax.set_xlabel(f"concentration [{std_curve.conc_unit:~P}]")
    ax.set_ylabel(f"signal [px]")
    ax.set_ylim(0, ax.get_ylim()[1])

def plot_concs(ax, unknowns, results):
    x = 0
    ticks = []
    tick_labels = []
    conc_unit = results.conc_unit

    for name, g in unknowns.groupby('name'):

        # Plot the mean
        dist = results.loc[name]['conc']
        ax.plot(
                [x, x],
                [0, dist['mean']],
                color=ucsf.blue[0],
                linewidth=7,
        )

        # Plot the individual data points
        i = g['within_std_curve']

        y1 = g['corrected_conc'][i]
        y2 = g['corrected_conc'][~i]

        x1 = [x] * len(y1)
        x2 = [x] * len(y2)

        ax.plot(
                x1, y1,
                marker='+',
                markerfacecolor='none',
                markeredgecolor=ucsf.dark_grey[0],
                linestyle='none',
        )
        ax.plot(
                x2, y2,
                marker='o',
                markerfacecolor='none',
                markeredgecolor=ucsf.dark_grey[0],
                linestyle='none',
        )

        ticks.append(x)
        tick_labels.append(f"{name}\n{dist['mean']:.2f} (±{dist['std']:.2f} sd) {conc_unit:~P}")

        x += 1

    ax.set_title('Unknowns')
    ax.set_xticks(ticks)
    ax.set_xticklabels(tick_labels)
    ax.set_ylabel(f"concentration [{conc_unit:~P}]")


def raw_data_style(lane):
    return dict(
            color=ucsf.dark_grey[0],
    )

def fit_style(lane):
    if lane.conf.is_unknown:
        label=f"{','.join(lane.conf.bands)}: {lane.conf.dilution}x"
    else:
        label=f"{lane.conf.band}: {lane.conf.conc}"

    return dict(
            color=ucsf.blue[0],
            label=label,
    )


if __name__ == '__main__':
    args = docopt.docopt(__doc__)
    conf = Config.parse_obj(nt.load(args['<conf>']))

    img = load_image(conf.image)
    lanes = fit_imagej_plots(conf, img)

    fig1 = plot_fits(lanes)
    if not args['--gui-only']:
        fig1.savefig(conf.image.parent / (conf.image.stem + '_fits.svg'))

    standards = find_standards(lanes)
    unknowns = find_unknowns(lanes)
    std_curve = fit_std_curve(standards)
    results = calc_unknown_concs(unknowns, std_curve)

    fig2 = plot_results(standards, unknowns, std_curve, results)
    if not args['--gui-only']:
        fig2.savefig(conf.image.parent / (conf.image.stem + '_results.svg'))

    if not args['--svg-only']:
        plt.show()

