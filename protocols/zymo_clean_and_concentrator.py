#!/usr/bin/env python3

import stepwise
import appcli
import autoprop

from inform import plural
from stepwise import pl, ul

@autoprop
class ZymoCleanAndConcentrator(appcli.App):
    """\
Usage:
    zymo_clean_and_concentrator [-v <µL>] [-s <µL>] [-d <when>] [-n <int>]

Options:
    -v --elute-volume <µL>      [default: 25]
        The volume of water to elute the RNA from the column in.

    -s --sample-volume <µL>
        The volume of the input RNA sample, in µL.

    -d --dnase-treatment <when>
        Add a DNase I treatment step to the protocol:

        pre: Before loading the sample onto the column
        post: After loading the sample onto the column.

    -n --num-samples <int>
        The number of samples to process.  This is only relevant when doing 
        DNase I treatment, in which case it will be used to help setup a master 
        mix.
"""
    __config__ = [
            appcli.DocoptConfig(),
    ]

    elute_volume_uL = appcli.param(
            '--elute-volume',
            cast=float,
            default=None,
    )
    sample_volume_uL = appcli.param(
            '--sample-volume',
            cast=float,
            default=None,
    )
    dnase_treatment = appcli.param(
            '--dnase-treatment',
            default=None,
    )
    num_samples = appcli.param(
            '--num-samples',
            cast=int,
            default=None,
    )

    def get_protocol(self):
        p = stepwise.Protocol()

        if self.dnase_treatment == 'pre':
            mm = stepwise.MasterMix("""\
                    Reagent                   Stock    Volume  MM?
                    ======================  =======  ========  ===
                    nuclease-free water              to 50 µL   +
                    DNA digestion buffer        10x      5 µL   +
                    DNase I [1]              1 U/µL      5 µL   +
                    RNA sample               <10 µg     40 µL   -
            """)

            step = "Setup a DNase I reaction for each sample:"
            if self.num_samples:
                step = f"Setup {plural(self.num_samples):# DNase I reaction/s}:"
                mm.num_reactions = self.num_samples

            if self.sample_volume_uL:
                mm['RNA sample'].volume = self.sample_volume_uL, 'µL'

            p += pl(step, mm)
            p += "Incubate at room temperature for 15 min."

            p.footnotes[1] = "Reconstitute lyophilized DNase I (#E1009-A; 250U) with 275 µL nuclease-free water and mix by gentle inversion.  Store frozen aliquots."

        if self.dnase_treatment == 'post':
            mm = stepwise.MasterMix("""\
                    Reagent                   Stock    Volume  MM?
                    ======================  =======  ========  ===
                    DNA digestion buffer        10x     75 µL   +
                    DNase I [1]              1 U/µL      5 µL   +
            """)
            if self.num_samples:
                mm.num_reactions = self.num_samples

            p += pl("Prepare 80 µL DNase I reaction mix for each sample:", mm)

        s = ul()
        p += pl("Purify RNA using Zymo Clean & Concentrator spin columns:", s)

        if self.sample_volume_uL:
            if self.sample_volume_uL < 50:
                s += f"Add {50 - self.sample_volume_uL:g} nuclease-free water to each sample"

            v = max(50, self.sample_volume_uL)
            s += f"Add {2*v:g} µL RNA binding buffer; mix."
            s += f"Add {3*v:g} µL >95% ethanol; mix."

        else:
            s += "If necessary, bring each sample to 50 µL."
            s += "Add 2 volumes RNA binding buffer; mix."
            s += "Add 3 volumes >95% ethanol; mix."

        s += f"Load sample on spin column."
        s += f"Spin 1 min, 16000g; discard flow-through."
        
        if self.dnase_treatment == 'post':
            s += f"Add 400 µL RNA wash buffer."
            s += f"Spin 30s, 16000g; discard flow-through."
            s += f"Add 80 µL DNase I reaction mix."
            s += f"Incubate at room temperature for 15 min."

        s += f"Add 400 µL RNA prep buffer."
        s += f"Spin 30s, 16000g; discard flow-through."
        s += f"Add 700 µL RNA wash buffer."
        s += f"Spin 30s, 16000g; discard flow-through."
        s += f"Add 400 µL RNA wash buffer."
        s += f"Spin 60s, 16000g; discard flow-through."
        s += f"Add {self.elute_volume_uL:g} µL nuclease-free water."
        s += f"Spin 30s, 16000g."

        return p


if __name__ == '__main__':
    app = ZymoCleanAndConcentrator.from_params()
    appcli.load(app)
    app.protocol.print()
