# first line: 147
@memory.cache
def cached_fit_bands(x, y, p0):
    return densiometry.fit_bands(x, y, p0)
