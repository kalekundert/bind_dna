# first line: 326
@memory.cache
def calc_growth_model(time_min, od600):
    """
    Thin wrapper for GrowthModel constructor that adds caching.
    """
    return GrowthModel(time_min, od600)
