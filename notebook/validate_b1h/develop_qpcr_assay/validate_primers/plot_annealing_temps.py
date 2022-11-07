#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_excel('20220620_optimize_qpcr_tm.xlsx')

# I know seaborn would probably be better for this...

def by_order(item):
    (gene, target, primers), group = item

    target_order = {"3’": 0, "mid": 1, "5’": 2}
    gene_order = {'GFP': 0, 'RFP': 1}

    return target_order[target], gene_order[gene]

groups = df.groupby(['Gene', 'Target', 'Primers'])

for (gene, target, primers), g in sorted(groups, key=by_order):
    plt.plot(g['Ta'], g['Intensity (px)'], label=f'{gene}, {target}, {primers.replace(",", "+")}')

plt.xlabel('$T_A$ (°C)')
plt.ylabel('product band (px)')
plt.legend()
plt.tight_layout()
plt.savefig('20220620_optimize_qpcr_tm_densiometry.svg')
plt.show()
