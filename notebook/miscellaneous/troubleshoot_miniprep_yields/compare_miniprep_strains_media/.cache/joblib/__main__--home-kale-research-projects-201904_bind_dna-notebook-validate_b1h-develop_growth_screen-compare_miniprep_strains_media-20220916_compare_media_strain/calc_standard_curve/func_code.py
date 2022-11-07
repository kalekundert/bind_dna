# first line: 72
@memory.cache
def calc_standard_curve(plots, ladders):
    std_px = []
    std_ng = []

    def find_closest_band(center, bands):
        best_dx = inf

        for band in bands:
            if band.ignore:
                continue

            dx = abs(center - band.fit_params.center_guess)
            if dx < best_dx:
                best_band = band
                best_dx = dx

        return best_dx, best_band

    i = 1

    for ladder in ladders:
        ladder.x, ladder.y = densiometry.curve_from_image(plots[ladder.i])

        p0 = [x.fit_params for x in ladder.bands]
        ladder.ps = densiometry.fit_bands(ladder.x, ladder.y, p0)

        for band, p in zip(ladder.bands, ladder.ps):
            if band.ignore:
                continue

            px = densiometry.gaussian_area(*p)
            ng = ladder.mass_ug * band.mass_ng_per_ug

            std_px.append(px)
            std_ng.append(ng)

    res = linregress(std_px, std_ng)
    return res, std_px, std_ng
