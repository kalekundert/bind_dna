# first line: 310
@memory.cache
def calc_max_growth_rates(df):
    cols = [
            'well', 'row_i', 'col_j',
            'sample', 'full_sample', 'group', 'reference', 'initial_od', 'outlier',
    ]
    gb = df.groupby(cols, as_index=False)
    return gb.apply(calc_max_growth_rate)
