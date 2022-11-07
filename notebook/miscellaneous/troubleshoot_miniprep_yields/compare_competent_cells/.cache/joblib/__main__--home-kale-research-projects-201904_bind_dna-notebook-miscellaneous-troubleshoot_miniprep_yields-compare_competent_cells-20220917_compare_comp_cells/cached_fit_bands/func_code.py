# first line: 137
@memory.cache
def cached_fit_bands(x, y, p0):
    return densiometry.fit_bands(x, y, p0)
