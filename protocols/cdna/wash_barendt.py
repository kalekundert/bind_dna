#!/usr/bin/env python3

import stepwise
import appcli
import autoprop

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
"""
    __config__ = [
            appcli.DocoptConfig(),
    ]

    volume_uL = appcli.param(
            '--volume',
            cast=float,
            default=15,
    )

    def get_protocol(self):
        s = stepwise.Step(
                "Remove unligated linker by ultrafiltration:",
                br='\n',
        )

        s += "Bring reaction to 500 µL with 8M urea."
        s += "Load onto a 100 kDa MWCO spin-filter [1]."
        s += "Spin 14000g, 15 min."
        s += "Wash with 500 µL 8M urea."
        s += "Wash with 500 µL nuclease-free water."
        s += "Wash with water again."
        s += "Wash with water again, but spin for 30 min [1]."
        s += """\
                Invert the filter into a clean tube and spin 1000g, 2 min to 
                collect ligated product in a volume of ≈15 µL."""

        if self.volume_uL > 15:
            s += f"Dilute to {self.volume_uL:g} µL with nuclease-free water"

        p = stepwise.Protocol()
        p += s

        p.footnotes[1] = 'Amicon UFC510024'
        p.footnotes[2] = 'Final urea concentration: ≈200 µM'

        return p

if __name__ == '__main__':
    app = WashBarendt.from_params()
    appcli.load(app)
    app.protocol.print()
