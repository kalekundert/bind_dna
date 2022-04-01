#!/usr/bin/env python3

import numpy as np

kd = np.array([
1.7e-10,
2.0e-8,
2.2e-9,
1.7e-8,
1.0e-10,
2.4e-9,
4.4e-9,
5.1e-9,
2.5e-10,
1.0e-10,
7.0e-8,
])

err = np.array([
0.07e-10,
0.78e-8,
0.30e-9,
0.08e-8,
0.01e-10,
0.19e-9,
0.42e-9,
0.47e-9,
0.91e-10,
0.23e-10,
1.21e-8,
])

kd_nM = kd / 1e-9
err_nM = err / 1e-9

for x in err_nM:
    print(f'{x:.2f}')

