#!/usr/bin/env python3

import wellmap
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

def load_cq(path):
    csv_dir = path.parent / path.stem
    csv_path = csv_dir / 'Quantification Cq Results.csv'
    return pd.read_csv(csv_path)[['Well', 'Cq']]

def format_label(slug):
    low, high = slug.split('_')
    return f'{low}–{high}°C'

path = Path('20200803_optimize_ta.toml')

df = wellmap.load(
        path,
        data_loader=load_cq,
        merge_cols={'well0': 'Well'},
)

print(df)

for label, df_plate in df.groupby('plate'):
    plt.plot(
            df_plate['temperature_C'],
            df_plate['Cq'],
            marker='+',
            linestyle='none',
            label=format_label(label),
    )

plt.legend(loc='best')
plt.savefig(path.with_suffix('.svg'))
plt.show()
