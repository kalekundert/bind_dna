#!/usr/bin/env python3


import stepwise
import appcli
import autoprop

from stepwise import pl, ul, dl, table
from stepwise_mol_bio import Main, BindableReagent, bind
from freezerbox import ReagentConfig
from appcli import Key, DocoptConfig
from itertools import groupby

@autoprop
class PlateAssay(Main):
    """\
    Usage:
        plate_assay.py <strains>...

    Arguments:
        <strains>
            The names of the strains to test.  These strains must be present in 
            the FreezerBox database.
    """
    __config__ = [
            DocoptConfig,
    ]

    class Strain(BindableReagent):
        antibiotics = appcli.param(ReagentConfig)

    strains = appcli.param(
            Key(DocoptConfig, '<strains>', cast=lambda tags: [
                PlateAssay.Strain(x) for x in tags
            ]),
            get=bind
    )

    def __init__(self, strains):
        self.strains = strains

    def get_protocol(self):
        # Reference:
        # - [Noyes2008]
        # - Expected results:
        #   - Supplementary Figure 5
        # - Protocol:
        #   - Supplementary Methods: "Zif268 activity on mutant finger 1 binding sites"

        p = stepwise.Protocol()

        # Histidine concentration:
        # - [Noyes2008] calls for 0.1% His.
        # - I assume this is 0.1% (w/v) = 1 mg/mL
        # - MW = 155.15 g/mol; 1 mg/mL = 6.44 mM
        # - My stock histidine is 100 mM, so it's more convenient to use a molar 
        #   concentration for this (e.g. 6 mM).  But it'd probably be smartest to just 
        #   copy [Noyes2008] exactly...
        #
        # Uracil concentration:
        # - [Noyes2008] consistently uses 200 µM uracil for permissive conditions.
        # - I could also just use my −His (+Ura) supplement.
        #   - Not sure what the resulting uracil concentration would be.
        #   - This approach is fine if I'm only making permissive plates, but if I'm 
        #     making both plates at once, it'll be better to add uracil.
        #
        # 3-AT concentrations:
        # - From previous experience, the AAA-target strain (s5) cannot survive with 
        #   even 0 mM 3-AT.  So there isn't much point titrating low concentrations of 
        #   3-AT (originally I did 1, 2, 4 mM).
        # - [Noyes2008] uses 5 and 10 mM 3-AT for selections.  Figure S5 furthermore 
        #   shows that the TGG-target strain can survive 50 mM 3-AT.
        # - My goal is just to show that the strains and plasmids perform as expected, 
        #   and eventually to evaluate my single-plasmid constructs.  For the latter, 
        #   it might be useful to do a titration to see how much more/less stringent 
        #   the single plasmid is.  But I think it's reasonable to use 10 mM 3-AT to 
        #   just ask, "do these plasmids/strains mostly work".

        by_antibiotics = lambda x: x.antibiotics
        strains_by_antibiotics = [
                (k, list(g))
                for k, g in groupby(
                    sorted(self.strains, key=by_antibiotics),
                    key=by_antibiotics,
                )
        ]

        def format_media(antibiotics):
            return '+'.join(['LB', *antibiotics])

        def format_strains(strains):
            return ','.join(x.tag for x in strains)

        f = "[Noyes2008] uses 2xYT instead of LB.  I don't think this detail matters, although the cells may grow faster in richer media."
        p += pl(f"Grow overnight cultures{p.add_footnotes(f)}:",
                dl(*(
                    (format_media(antibiotics), format_strains(strains))
                    for antibiotics, strains in strains_by_antibiotics
                )),
                ul(
                    "Use fresh colonies (i.e. either transformed or restreaked within 1 week).",
                    "Incubate at 37°C overnight with shaking."
                ),
        )

        # Acclimate to NM:
        # - [Noyes2008] has a 1h growth in NM + uracil to "acclimate the cells to NM 
        #   media".
        # - This step doesn't make total sense to me:
        #   - Why +Ura but not +His?
        #   - Why not do the whole day culture in rich NM media?  Maybe they grow 
        #     slowly?  I'm tempted to do this, but the OD for "log phase" will be 
        #     different in NM media.
        # - I'll switch the day culture step to NM once I measure the OD for log phase.  
        #   Until then, I'll just keep this the same as I did the first time.
        p += pl(
                "Grow day cultures:",
                ul(
                    f"Inoculate 3 mL {'/'.join(format_media(x) for x, _ in strains_by_antibiotics)} with 30 µL saturated overnight culture.",
                    "Incubate at 37°C with shaking until OD≈0.4 (5-6h)",
                ),
        )

        # Number of washes:
        # - [Noyes2008] calls for 4 washes.
        # - I only did 2 washes the first time I did this protocol, and the selection 
        #   was still very strict.
        # - Probably best to err on the side of caution, though.
        p += pl(
                "Wash cells with minimal media:",
                ul(
                    pl(
                        "Repeat 4x:",
                        ul(
                            "Spin 3500g, 3 min, 4°C.",
                            "Resuspend in 1 mL NM.",
                        ),
                        br='\n',
                    ),
                    "Dilute to OD=0.1",
                ),
        )

        # Serial dilution:
        # - [Noyes2008]_ just plated ≈500 cfus for each condition.
        # - I plated 6 10x dilutions for each construct.
        # - I think that the serial dilution is a better approach, because it makes it 
        #   easier to quantify/interpret large changes in fitness.
        p += stepwise.load('serial 90µL 0.1OD / 10 6 -m cells -d NM')

        # Antibiotics:
        # - Once I'm testing my single plasmid, I should plate everything on Carb 
        #   plates, and then also plate the two-plasmid system on Carb+Kan plates.
        p += pl(
                "Plate 5 µL of each dilution in triplicate on permissive/selective plates with the appropriate antibiotics.",
                ul(
                    "Don't draw a grid on the plate directly.  It makes it hard to count colonies.  Instead, place a guide underneath the plate and spot based on that.",
                )
        )

        return p


if __name__ == '__main__':
    PlateAssay.main()
