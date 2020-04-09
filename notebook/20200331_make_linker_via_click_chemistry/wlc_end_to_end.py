#!/usr/bin/env python3

import numpy as np
from math import exp
from scipy.optimize import root_scalar
from matplotlib.pyplot import *

# l: max length
# r2: mean squared end-to-end distance
# p: persistence length

monomer_lengths_A = {
        'spacer18': 6 * 3.5,
        'ssDNA': 6.3,
}
persistence_lengths_A = {
        'spacer18': 3.8,
        'ssDNA': 25.0,
}

# The arm needs to have at least this length.
path_length_A = 33.4 + 62.2 + 44.0 + 45.3
end_to_end_dist_A = 83.4

for monomer in monomer_lengths_A:
    r = end_to_end_dist_A; r2 = r**2
    P = persistence_lengths_A[monomer]

    # Really, we're in the L >> P domain, so I don't need to do all this fancy 
    # solving.
    #
    #   R**2 = 2 * P * L
    #   L = R**2 / 2*P

    x0 = r2/(2*P)   # Limit if L >> P
    x1 = r2         # Limit if L << P

    f_r2 = lambda L: 2*P*L * (1 - (P/L)*(1 - np.exp(-L/P)))
    f = lambda L: f_r2(L) - r2

    root = root_scalar(f, x0=x0, x1=x1)

    polymer_len = max(root.root, path_length_A)
    n_monomers = polymer_len / monomer_lengths_A[monomer]

    debug(monomer, x0, root, path_length_A, polymer_len, n_monomers)
    print()

