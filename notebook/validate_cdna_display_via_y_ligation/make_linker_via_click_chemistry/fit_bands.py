#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import more_itertools as mi

from math import *
from scipy.optimize import curve_fit
from pathlib import Path

img_path = Path('20200804_click_linker_658_plots.png')
labels = [
        dict(
            name='16h',
            bands=2,
        ),
        dict(
            name='4h',
            bands=2,
        ),
        dict(
            name='1h',
            bands=2,
        ),
        dict(
            name='15m',
            bands=2,
        ),
        dict(
            name='o126 only',
            bands=1,
        ),
]

def find_divisions(img):
    """
    ImageJ produces a single images where the plots for each lane are stacked 
    on top of each other.  The first step is to figure out where these plots 
    start and stop, so we can break each one out and deal with it individually.  
    The basic approach to doing this is to look for horizontal black lines.

    Some complications:

    - There's no divider between the first plot and a bit of header text at the 
      top of the image.  I empirically worked out that the text takes up 16 
      pixels, and account for that explicitly.

    - There are two lines of black pixels at the bottom of image, for no 
      apparent reason.
    """

    is_divider = np.all(img, axis=1)
    divisions = np.arange(len(is_divider))[is_divider]
    divisions = sorted([16, *divisions])[1:-1]
    return divisions

def iter_plots(labels, img, divs):
    for label, (i, j) in mi.zip_equal(labels, mi.pairwise(divs)):
        yield label, img[i+1:j,1:-1]

def path_from_image(img):
    x = np.arange(img.shape[1])
    y = np.argmax(img, axis=0)
    ii = y > 0
    return x[ii], y[ii]

def fit_bands(x, y, n=2):
    """
    Fit Gaussian curves to each band in the data.

    The following assumptions are made:

    - The number of bands is given, and is either 1 or 2.
    - All of the Gaussians have the same width.
    """

    # The parameters are formatted like so:
    # (a1, b1, a2, b2, ..., aN, bN, c)

    def iter_params(p):
        for pi in mi.chunked(p[:-1], 2):
            yield (*pi, p[-1])

    def f(x, *p):
        return sum((gaussian(x, *pi) for pi in iter_params(p)))

    # The initial guess has to be pretty good, or the optimize doesn't seem to 
    # find a good solution.  This is fine as long as I'm checking the results 
    # by hand, but it might be more robust/generalizable to use a global 
    # optimization routine instead of a non-linear curve fit.
    if n == 1:
        p0 = (300, 550, 10)
    elif n == 2:
        p0 = (300, 550, 50, 600, 10)

    bounds = (
            n*(  0, min(x)) + (-20,),
            n*(500, max(x)) + ( 20,),
    )

    p, cov = curve_fit(
            f, x, y,
            p0=p0,
            bounds=bounds,
    )
    return tuple(iter_params(p))

def gaussian(x, a, b, c):
    return a * np.exp(-(x - b)**2 / (2 * c**2))

def gaussian_area(a, b, c):
    # https://en.wikipedia.org/wiki/Gaussian_function
    return a * c * sqrt(2 * pi)

img = 1 - plt.imread(img_path)
divs = find_divisions(img)

n_rows = len(divs) - 1
n_cols = 1
ax = None

for i, (label, plot) in enumerate(iter_plots(labels, img, divs), 1):
    x, y = path_from_image(plot)
    params = fit_bands(x, y, n=label['bands'])

    if not ax:
        ax = plt.subplot(n_rows, n_cols, i, label=label['name'])
    else:
        plt.subplot(n_rows, n_cols, i, label=label['name'], sharex=ax)

    plt.plot(x, y, label=label['name'])

    y_fits = []
    a_fits = []

    for p in params:
        y_fit = gaussian(x, *p);    y_fits.append(y_fit)
        a_fit = gaussian_area(*p);  a_fits.append(a_fit)
    for y_fit, a_fit in zip(y_fits, a_fits):
        plt.plot(x, y_fit, label=f'A={a_fit:.1f} px²; {100 * a_fit / sum(a_fits):.2f}%')
    if label['bands'] > 1:
        plt.plot(x, sum(y_fits))

    plt.legend(loc='upper left')

plt.gcf().set_size_inches(6 * n_cols, 2 * n_rows)
plt.tight_layout()
plt.savefig(img_path.with_suffix('.svg'))
plt.show()
