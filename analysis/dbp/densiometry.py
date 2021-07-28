#!/usr/bin/env python3

import numpy as np
import autoprop

from scipy.optimize import curve_fit
from more_itertools import zip_equal, pairwise, chunked, flatten
from math import sqrt, pi

@autoprop
class BandFitParams:

    def __init__(self, center_guess, *,
            center_min=None,
            center_max=None,

            width_guess=10,
            width_min=1,
            width_max=None,

            height_guess=None,
            height_min=0,
            height_max=None,

            baseline_guess=0,
            baseline_min=0,
            baseline_max=None,
    ):
        self._center_guess = center_guess
        self._center_min = center_min
        self._center_max = center_max

        self._width_guess = width_guess
        self._width_min = width_min
        self._width_max = width_max

        self._height_guess = height_guess
        self._height_min = height_min
        self._height_max = height_max

        self._baseline_guess = baseline_guess
        self._baseline_min = baseline_min
        self._baseline_max = baseline_max

    def get_center_guess(self):
        return self._center_guess

    def get_center_min(self, x):
        return self._center_min or min(x)

    def get_center_max(self, x):
        return self._center_max or max(x)

    def get_width_guess(self):
        return self._width_guess

    def get_width_min(self):
        return self._width_min

    def get_width_max(self):
        return self._width_max or 5 * self.width_guess

    def get_height_guess(self, x, y):
        return self._height_guess or max(
            y[int(self.get_center_min(x)):int(self.get_center_max(x))]
        )

    def get_height_min(self):
        return self._height_min

    def get_height_max(self, y):
        return self._height_max or 3 * max(y)

    def get_baseline_guess(self):
        return self._baseline_guess

    def get_baseline_min(self):
        return self._baseline_min

    def get_baseline_max(self, y):
        return self._baseline_max or max(y)

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
    for label, (i, j) in zip_equal(labels, pairwise(divs)):
        yield label, img[i+1:j,1:-1]

def path_from_image(img):
    x = np.arange(img.shape[1])
    y = np.argmax(img, axis=0)
    ii = y > 0
    return x[ii], y[ii]

def fit_bands(x, y, bands):
    """
    Fit Gaussian curves to each band in the data.
    """

    p0 = list(flatten(
            (
                b.get_height_guess(x, y),
                b.center_guess,
                b.width_guess,
                b.baseline_guess,
            )
            for b in bands
    ))
    bounds_min = list(flatten(
            (
                b.height_min,
                b.get_center_min(x),
                b.width_min,
                b.baseline_min,
            )
            for b in bands
    ))
    bounds_max = list(flatten(
            (
                b.get_height_max(y),
                b.get_center_max(x),
                b.width_max,
                b.get_baseline_max(y),
            )
            for b in bands
    ))
    p, cov = curve_fit(
            gaussian_sum,
            x, y,
            p0=p0,
            bounds=(bounds_min, bounds_max),
    )
    return list(chunked(p, 4))

def gaussian(x, a, b, c, d):
    return a * np.exp(-(x - b)**2 / (2 * c**2)) + d

def gaussian_sum(x, *p):
    return sum((gaussian(x, *p_i) for p_i in chunked(p, 4)))

def gaussian_area(a, b, c, d):
    # https://en.wikipedia.org/wiki/Gaussian_function
    return a * c * sqrt(2 * pi)
