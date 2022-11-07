#!/usr/bin/env python3

import wellmap
import wellmap_qpcr
import matplotlib.pyplot as plt

df = wellmap.load(
        '20220507_try_b1h_qpcr_assay.toml',
        data_loader=wellmap_qpcr.load_cq,
        merge_cols=True,
        path_guess='{0.stem}',
)

df = df\
        .groupby(['primers', 'target', 'architecture'])\
        .apply(wellmap_qpcr.agg_cq)

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

fig, ax = plt.subplots()
ticks = []
tick_labels = []

for i, (arch, row) in enumerate(df.iterrows()):
    ax.plot(
            [i, i],
            [0, row['fold_change']],
            linewidth=5,
            solid_capstyle='butt',
    )

    ticks.append(i)
    tick_labels.append(arch)

ax.set_xlim(-0.5, i+0.5)
ax.set_xticks(ticks, tick_labels)
ax.set_xlabel('single-plasmid scaffold')

ax.set_ylabel('gene expression\n2TGG rel. to 2AAA')

plt.savefig('20220507_try_b1h_qpcr_assay.svg')
plt.show()

# Plots:
# - [x] Cq heatmap
#   - wellmap_qpcr
#
# - [x] gene expression
#   - bar plot
#   - ΔCt or ΔΔCt?
