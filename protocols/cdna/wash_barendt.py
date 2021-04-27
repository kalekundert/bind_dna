#!/usr/bin/env python3

import stepwise
import appcli
import autoprop

from stepwise import pl, ul

@autoprop
class WashBarendt(appcli.App):
    """\
Remove unligated linker by ultrafiltration.

Usage:
    wash_barendt [-v <µL>]

Options:
    -v --volume <µL>        [default: ${app.volume_uL}]
        The volume to dilute the purified mRNA to.  Note that this must be 
        greater than 15 µL, since that is the dead volume of the spin filter.

    -m --mwco <kDa>         [default: ${app.mwco_kDa}]
        The MWCO of the spin filter.  According to Millipore Sigma, this should 
        be 2-3x smaller than the molecular weight of the ligated product: 
        https://tinyurl.com/4ffxu8zb
"""
    __config__ = [
            appcli.DocoptConfig(),
    ]

    volume_uL = appcli.param(
            '--volume',
            cast=float,
            default=15,
    )
    mwco_kDa = appcli.param(
            '--mwco',
            cast=int,
            default=100,
    )

    def get_protocol(self):
        p = stepwise.Protocol()
        p += pl(
                "Remove unligated linker by ultrafiltration:",
                s := ul(
                    "Bring reaction to 500 µL with 8M urea.",
                    f"Load onto a {self.mwco_kDa} kDa MWCO spin-filter [1].",
                    "Spin 14000g, 15 min.",
                    "Wash with 500 µL 8M urea.",
                    "Wash with 500 µL nuclease-free water.",
                    "Wash with water again.",
                    "Wash with water again, but spin for 30 min.",
                    "Invert the filter into a clean tube and spin 1000g, 2 min to collect ligated product in a volume of ≈15 µL.",
                ),
        )
        if self.volume_uL > 15:
            s += f"Dilute to {self.volume_uL:g} µL with nuclease-free water"

        p.footnotes[1] = 'https://tinyurl.com/4ffxu8zb'
        p.footnotes[2] = 'Final urea concentration: ≈200 µM'

        return p

if __name__ == '__main__':
    app = WashBarendt.from_params()
    appcli.load(app)
    app.protocol.print()
