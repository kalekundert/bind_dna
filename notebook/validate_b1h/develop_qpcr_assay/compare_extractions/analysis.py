#!/usr/bin/env python3

# Analysis:
# - See raw Cq for −RT controls
# - Show Cq curves, grouped in some reasonable way.

import wellmap
import wellmap_qpcr
import matplotlib.pyplot as plt
import pandas as pd

from pathlib import Path

pd.set_option("display.max_columns", 100)
pd.set_option("display.max_rows", 100)
pd.set_option("display.width", 10000)

wellmap_path = Path('20220603_compare_rna_purification.toml')
df = wellmap.load(
        wellmap_path,
        data_loader=wellmap_qpcr.load_cq,
        merge_cols=True,
        path_guess='{0.stem}',
)

df = df\
        .groupby(['primers', 'target', 'rt', 'extraction'])\
        .apply(wellmap_qpcr.agg_cq)

print(df)

# df = wellmap_qpcr.calc_Δcq(
#         df_expt=df.query('primers == "sr1,sr2"').droplevel('primers'),
#         df_ref =df.query('primers == "sr3,sr5"').droplevel('primers'),
# )
df = wellmap_qpcr.calc_Δcq(
        df_expt=df.loc['sr1,sr2'],
        df_ref =df.loc['sr3,sr5'],
)

df = wellmap_qpcr.calc_ΔΔcq(
        df_expt=df.loc['2TGG'],
        df_ref =df.loc['2AAA'],
)

fig, ax = plt.subplots(figsize=(5, 8))
ticks = []
tick_labels = []

for i, (index, row) in enumerate(df.iterrows()):
    ax.plot(
            [i, i],
            [0, row['fold_change']],
            linewidth=5,
            solid_capstyle='butt',
    )

    rt, extraction = index
    label = f'{extraction}, {"−+"[rt]}RT'

    ticks.append(i)
    tick_labels.append(label)

ax.set_xlim(-0.5, i+0.5)
ax.set_xticks(ticks, tick_labels, rotation='vertical')
ax.set_ylabel('gene expression\n2TGG rel. to 2AAA')

plt.tight_layout()
plt.savefig(wellmap_path.with_suffix('.svg'))
plt.show()
