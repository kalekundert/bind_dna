#!/usr/bin/env python3

import stepwise
import appcli
import autoprop

from stepwise import Quantity, pl, ul
from stepwise_mol_bio import Main, DirectDilution, Gel, Stain
from freezerbox.stepwise.dilute import Dilute
from appcli import Key, DocoptConfig
from inform import plural

@autoprop.cache
class MeasureProteinConc(Main):
    """\
Measure protein concentration running a gel (with a standard curve), 
visualizing with SYPRO orange, and comparing band intensities.

Usage:
    measure_protein_conc <unknown>... [-s <steps>] [-S <steps>] [-v <µL>]
        [-x <fold>] [-d <dilute>] [-g <gel>]

Argument:
    <unknown>
        The name of one or more protein samples of unknown concentration.

Options:
    -s --unknown-steps <int>    [default: ${app.unknown_steps}]
        How many dilutions of the unknown sample to make.  You may want to 
        lower this, e.g., to fit more samples on a single gel.  This option 
        takes precedence over `--steps`, if both are specified.

    -S --steps <int>            [default: ${app.standard_steps}]
        How many dilutions of both the unknown and reference samples to make.  
        Use `--unknown-steps` to specify different numbers of steps for the 
        reference and unknown samples.

    -v --unknown-volume <µL>
        What volume of each unknown dilution to prepare.  This should be at 
        least 5 µL, in order to have enough to load on the gel.  20 µL is 
        recommended, to avoid having to pipet small volumes. 

    -x --unknown-dilution <x>  [default: ${app.unknown_dilution}]
        The greatest dilution of the unknown sample to run on the gel.  The 
        default is chosen to match the standard curve.

    -d --unknown-dilution-cmd <cmd>
        Use the given command to make dilutions of the unknown sample.  Any 
        occurrences of "{unk}" in the command will be replace by the 
        comma-separated names of the unknown proteins.  By default, direct 
        dilutions are made (since they are more accurate than serial 
        dilutions). 

    -g --gel-cmd <cmd>
        Use the given command to run the gel.  Any occurrences of "{n}" in the 
        command will be replaced by the combined number of standard and unknown 
        dilutions to run.  By default, 5 µL of each dilution are loaded on a 
        4-12% Bolt/MES SDS PAGE gel.

References:

1. Knight MI & Chambers PJ. Problems associated with determining protein 
   concentration. Mol Biotechnol 23, 19–28 (2003).
"""
    __config__ = [
            DocoptConfig,
    ]

    unknowns = appcli.param('<unknown>')
    unknown_steps = appcli.param(
            Key(DocoptConfig, '--unknown-steps', cast=int),
            Key(DocoptConfig, '--steps', cast=int),
            default=7,
    )
    standard_steps = appcli.param('--steps', cast=int, default=7)
    unknown_volume = appcli.param('--unknown-volume', cast=float, default=20)
    unknown_dilution = appcli.param('--unknown-dilution', cast=float, default=1/8)
    unknown_dilution_cmd = appcli.param('--unknown-dilution-cmd', default=None)
    gel_cmd = appcli.param('--gel-cmd', default=None)

    def get_protocol(self):
        p = stepwise.Protocol()
        p += self.standard_dilution_protocol
        p += self.unknown_dilution_protocol
        p += self.electrophoresis_protocol
        p += stepwise.load("stain sypro-orange")
        p += f"Quantify band intensities, then calculate concentrations by linear regression."
        return p

    def get_standard_dilution_protocol(self):
        p = stepwise.Protocol()
        p += pl(
                "Prepare 500 µL 80 µg/mL BSA:",
                ul(
                    "480 µL water",
                    "20 µL 2 mg/mL BSA standard (from a fresh ampule)",
                ),
        )

        dd = DirectDilution(
                volume=Quantity(20, 'µL'),
                steps=self.standard_steps,
        )
        dd.set_conc_high_low('80 µg/mL', '10')
        dd.material = "BSA standard"

        p += pl(
                "Prepare the standard curve as follows:",
                dd.dilution_table,
        )

        return p

    def get_unknown_dilution_protocol(self):
        unk = ', '.join(self.unknowns)

        if self.unknown_dilution_cmd:
            return stepwise.load(self.unknown_dilution_cmd.format(unk=unk))

        dd = DirectDilution(
                volume=Quantity(self.unknown_volume, 'µL'),
                steps=self.unknown_steps,
        )
        dd.set_conc_high_low('1x', self.unknown_dilution)
        dd.material = "unknown"

        p = stepwise.Protocol()
        p += pl(
                f"Prepare the unknown {plural(self.unknowns):sample/s} ({unk}) as follows:",
                dd.dilution_table,
        )
        return p

    def get_electrophoresis_protocol(self):
        n = self.standard_steps + len(self.unknowns) * self.unknown_steps

        if self.gel_cmd:
            return stepwise.load(self.gel_cmd.format(n=n))

        gel = Gel('bolt/mes')
        gel.sample_mix.num_reactions = n
        gel.sample_mix['sample'].name = "standards/samples"
        gel.sample_mix['sample'].volume = '5 µL'
        gel.sample_mix['sample'].stock_conc = None
        gel.stain = None

        return gel.protocol


if __name__ == '__main__':
    MeasureProteinConc.main()
