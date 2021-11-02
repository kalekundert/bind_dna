#!/usr/bin/env python3

import stepwise
import appcli
import autoprop

from stepwise import Quantity, Q, pl, ul
from stepwise_mol_bio import Main, Lyophilize
from appcli import DocoptConfig
from inform import plural

@autoprop.cache
class SplitLigation(Main):
    """\
Ligate a puromycin linker to an mRNA transcript via split ligation.

This protocol is based on [Doshi2014].

Usage:
    split_ligation_doshi <mrna> [<linker>] [<splint>] [--debug]

Options:
    -d --debug
        Save aliquots after each step, to ensure that the reaction worked as 
        expected, and include steps describing how to setup −ligation, 
        −digestion, and −purification controls.
"""

    __config__ = [
            DocoptConfig,
    ]

    mrna = appcli.param('<mrna>')
    splint = appcli.param('<splint>', default='o277')
    linker = appcli.param('<linker>', default='o278')
    debug = appcli.param('--debug', default=False)

    def get_protocol(self):
        p = stepwise.Protocol()
        if self.debug:
            p += self.control_protocol
        p += self.annealing_protocol
        p += self.ligation_protocol
        p += self.digestion_protocol
        p += self.purification_protocol
        p += self.concentration_protocol
        return p

    def get_control_protocol(self):
        p = stepwise.Protocol()

        # These concentrations are not exactly the same as in the reactions.  
        # The actual concentrations change after each step, so it wouldn't be 
        # easy to match them exactly.  Instead, I chose round numbers that are 
        # just slightly less dilute than the actual concentrations.

        # I'm not concerned that these volumes are too small to pipet 
        # accurately.  I'm not going to be comparing the intensities of these 
        # bands to anything.

        p += pl(
                "Prepare a 500 nM mRNA-only control:",
                ul(
                    f"4.75 µL nuclease-free water",
                    f"0.25 µL 10 µM {self.mrna}",
                ),
        )
        p += pl(
                "Prepare a 2 µM linker-only control:",
                ul(
                    f"4.90 µL nuclease-free water",
                    f"0.10 µL 100 µM {self.linker}",
                ),
        )
        p += pl(
                "Prepare a 2 µM splint-only control:",
                ul(
                    f"4.90 µL nuclease-free water",
                    f"0.10 µL 100 µM {self.splint}",
                ),
        )
        return p

    def get_annealing_reaction(self):
        # reaction volume:
        # - [Doshi2014] calls for a total reaction volume of 60 µL, but that 
        #   includes the 10 µL of reagents that will be added after annealing.  
        #   So the volume at this step is only 50 µL.
        #
        # mRNA volume:
        # - [Doshi2014] calls for 5 µg mRNA in a 60 µL reaction.
        # - Nb mRNA used by [Doshi2014] has MW=143365.24
        # - That works out to: 581 nM (final conc)
        # - 10 µM seems to be a reasonable stock concentration, based on my 
        #   previous IVT reactions with comparable-length templates.
        # 
        # oligo volumes:
        # - [Doshi2014] calls for 2.67 µM final concentrations of both oligos.
        rxn = stepwise.MasterMix("""\
                Reagent               Stock      Volume  MM?
                ===================  ======  ==========  ===
                nuclease-free water          to 50.0 µL   +
                ssRNA ligase buffer     10x      6.0 µL   +
                mRNA                  10 µM      3.5 µL
                splint               100 µM      1.6 µL   +
                linker               100 µM      1.6 µL
        """)
        rxn['mRNA'].name = f"{self.mrna} (mRNA)"
        rxn['splint'].name = f"{self.splint} (splint)"
        rxn['linker'].name = f"{self.linker} (puro)"

        # control aliquot volumes:
        # - I just need enough to run a gel, which is 3.2 µL (see notes for 
        #   `urea/ligate/doshi2014` gel preset).
        # - Round to 5 µL to be safe.

        if self.debug:
            rxn.hold_ratios.volume += 10, 'µL'

        return rxn

    def get_annealing_protocol(self):
        p = stepwise.Protocol()
        rxn = self.annealing_reaction

        p += pl(
                f"Setup {plural(rxn.num_reactions):# annealing reaction/s}{p.add_footnotes('[Doshi2014] DOI:10.1038/srep06760')}:",
                rxn,
        )

        # This annealing protocol is not specified by [Doshi2014].  I got this 
        # protocol from Erkin.  I imagine that a lot of protocols would work, 
        # though.
        p += pl(
                "Incubate as follows:",
                ul(
                    "82°C for 2 min",
                    "Cool to 25°C over 5 min",
                ),
        )
        return p

    def get_ligation_protocol(self):
        p = stepwise.Protocol()
        k = 1

        if self.debug:
            k = (self.annealing_reaction.volume - Q('5 µL')) / Q('50 µL')
            p += pl(
                    "Remove a 5 µL −ligase aliquot.  Store at −20°C.",
            )

        # enzyme volume:
        # - [Doshi2014] specifies this volume of enzyme, but doesn't specify 
        #   the stock concentration (or the catalog number).  I'm assuming that 
        #   they used NEB T4 ssRNA ligase (M0204) since most of their other 
        #   reagents are from NEB.

        p += pl(
                "Add the following to the reaction:",
                ul(
                    f"{6*k:.2g} µL 10 mM ATP",
                    f"{1*k:.2g} µL RNase inhibitor (murine)",
                    f"{3*k:.2g} µL 10 U/µL ssRNA ligase (NEB M0204)",
                ),
        )
        p += "Incubate at room temperature overnight."
        return p

    def get_digestion_reaction(self):
        return stepwise.MasterMix("""\
                Reagent                Stock     Volume  MM?
                ====================  ======  =========  ===
                nuclease-free water           to 240 µL   +
                λ exonuclease buffer     10x      24 µL   +
                λ exonuclease         5 U/µL       5 µL   +
                ligation reaction                 60 µL
        """)

    def get_digestion_protocol(self):
        p = stepwise.Protocol()
        rxn = self.digestion_reaction

        if self.debug:
            p += pl(
                    "Remove a 5 µL −digestion aliquot.  Store at −20°C.",
            )

        p += pl(
                f"Setup {plural(rxn.num_reactions):# digestion reaction/s}:",
                rxn,
        )
        p += "Incubate at 37°C for 45 min."
        p += pl(
                "Denature the exonuclease:",
                ul(
                    "Add 60 µL 1M tris-HCl, pH 7.5",
                    "Incubate at 65°C for 5 min.",
                ),
        )
        return p

    def get_purification_protocol(self):
        p = stepwise.Protocol()

        if self.debug:
            # Aliquot volume:
            # - The digestion reaction is 4x more dilute than the annealing/ 
            #   ligation reaction.
            # - So I would need a 20 µL aliquot to get the same amount of 
            #   material.
            # - I think 10 µL is the most I can load on a gel without a 
            #   concentration step, so that what I went for.
            p += pl(
                    "Remove a 10 µL −purification aliquot.  Store at −20°C.",
            )

        # Binding buffer:
        # - The NEB protocol calls for the RNA to be in "lysis/binding buffer":
        #   - 100 mM tris-HCl, pH=7.5
        #   - 500 mM LiCl
        #   - 0.5% LDS
        #   - 1 mM EDTA
        #   - 5 mM DTT
        #
        # - I'm not sure if it's really necessary to use this buffer, but the 
        #   safe thing is to assume that it is.
        #   
        #   - The LDS could help prevent non-specific binding.  It could also 
        #     just be there to help lyse cells, though.
        #
        #   - LiCl is used (at higher concentrations, e.g. 2.5-4M) to 
        #     specifically precipitate RNA and not DNA.  Still, 500 mM is a 
        #     lot, and it might help the RNA get close to the bead and/or 
        #     anneal with the oligo-dT.
        #
        #     https://www.thermofisher.com/us/en/home/references/ambion-tech-support/rna-isolation/general-articles/the-use-of-licl-precipitation-for-rna-purification.html
        #
        # - The ligation/digestion reaction already has >200 mM tris, so I 
        #   don't need any more of that.
        #
        # - Perhaps I could make a 10x mix of the remaining chemicals, and add 
        #   it to the reaction.  I have all the necessary chemicals on hand.

        p += pl(
                "Prepare 100 µL 10x lysis/binding buffer (−tris):",
                ul(
                    "20.5 µL water",
                    "50 µL 10M LiCl",
                    "25 µL 20% (w/v) LDS",
                    "2 µL 500 mM EDTA",
                    "2.5 µL 1 M DTT",
                ),
        )

        # Wash buffer volume:
        # - NEB specifies volumes of wash buffer to use for 50 µL and 100 µL 
        #   beads.  [Doshi2014] calls for 70 µL beads, so I interpolated.

        f = 'https://tinyurl.com/3mmtjk5h'
        p += pl(
                f"Purify using magnetic oligo-dT(25 beads){p.add_footnotes(f)}:",
                ul(
                    "Equilibrate 70 µL beads in 200 µL lysis/binding buffer.",
                    "Add 33 µL 10x lysis/binding buffer (−tris) to the digestion reaction.",
                    "Remove buffer from beads and immediately add digestion reaction.",
                    "Agitate at room temperature for 10 min.",
                    "Discard supernatant.",
                ),
                ul(
                    pl(
                        "Repeat 2x:",
                        ul(
                            "Add 350 µL wash buffer 1",
                            "Agitate at room temperature for 1 min.",
                            "Discard supernatant.",
                        ),
                        br='\n',
                    )
                ),
                ul(
                    pl(
                        "Repeat 2x:",
                        ul(
                            "Add 350 µL wash buffer 2",
                            "Agitate at room temperature for 1 min.",
                            "Discard supernatant.",
                        ),
                        br='\n',
                    ),
                ),
                ul(
                    "Add 350 µL low salt buffer",
                    "Agitate at room temperature for 1 min.",
                    "Discard supernatant.",
                ),
                ul(
                    pl(
                        "Repeat 2x:",
                        ul(
                            "Add 30 µL nuclease-free water.",
                            "Vortex gently to resuspend beads.",
                            "Agitate at 50°C for 2 min.",
                            "Recover supernatant.",
                        ),
                        br='\n',
                    ),
                ),
        )
        return p

    def get_concentration_protocol(self):
        # [Doshi2014] concentrates by precipitation, but doing so adds salt to 
        # the reaction.  Maybe that's not a bad thing, though...
        lyo = Lyophilize(volume=Quantity(6.5, 'µL'))
        return lyo.protocol


if __name__ == '__main__':
    SplitLigation.main()

    



                    



# Anneal:
# - 581 nM mRNA (paper calls for 5 µg; mRNA MW=143365.24; volume=60 µL)
# - 2.67 µM (final) splint
# - 2.67 µM (final) puromycin linker
#   - ≈5x excess of splint/linker is expensive, but won't cause problems 
#     because it'll be destroyed in the exonuclease step.
# - 1x ssRNA ligase buffer
# - Water to 60 µL
#
# Annealing protocol unspecified, but I assume that 95°C for 2 min is 
# reasonable.
#
# Ligate:
# - 1 mM ATP
# - 1 µL RNase inhibitor (NEB)
# - 3 µL ssRNA ligase
#
# Incubate at room temperature overnight.
#
# Digest:
# - 5 µL lambda exonuclease (NEB)
# - 1x lambda exonuclease buffer
# - water to 240 µL
#
# Incubate at 37°C for 45 min.
# Add 60 µL 1M tris, pH 7.5
# Incubate at 65°C for 5 min.
#
# Bead purification
# - Follow manufacturer's protocol (NEB, OligodT₂₅ beads)
# - Elute twice in 30 µL 10 mM Tris-HCl, pH 7.5
#   - Might do nuclease free water instead.
# 
# Concentrate:
# - Protocol calls for precipitation.
# - Might do lyophilization instead.
#
