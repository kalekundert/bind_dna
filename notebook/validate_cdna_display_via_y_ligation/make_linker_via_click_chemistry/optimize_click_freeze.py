#!/usr/bin/env python3

import stepwise
from stepwise import pl, ul

p = stepwise.Protocol()

p += stepwise.load('cond optimize_click_freeze_cond.xlsx')

rxn = stepwise.MasterMix("""\
        Reagent   Stock  Volume  MM?
        =======  ======  ======  ===
        PBS          2x    2 µL   -
        o126     400 µM    1 µL   +
        o233     400 µM    1 µL   +
""")
rxn.num_reactions = 2
rxn.extra_percent = 0
rxn.hold_ratios.volume = '3 µL'

p += pl(
    f"Setup {rxn.num_reactions} click reactions:",
    rxn,
)
p += "Divide each reaction in three."

p += pl(
    "Incubate the reactions as follows:",
    ul(
        "−20°C overnight",
        "25°C for 4h, −20°C overnight",
        "25°C overnight",
    ),
)

p += stepwise.load(f"gel urea {rxn.num_reactions + 2} -c 1730 -S")
p.insert_footnotes(
        pl(
            "The product (o236) has MW = 17285 g/mol, so:",
            "100 µM = 1730 ng/µL",
            br='\n',
        ),
)
p += stepwise.load("stain sybr-green-ii/page")
p += stepwise.load("laser blue red")

p.print()
