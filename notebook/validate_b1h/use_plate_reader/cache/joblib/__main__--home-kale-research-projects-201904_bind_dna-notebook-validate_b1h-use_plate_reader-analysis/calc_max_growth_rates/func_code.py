# first line: 177
@memory.cache
def calc_max_growth_rates(df):
    cols = [
            'well', 'row_i', 'col_j',
            'system', 'strain', 'target', 'selection', 'dilution',
    ]
    gb = df.groupby(cols, as_index=False)
    return gb.apply(calc_max_growth_rate)
