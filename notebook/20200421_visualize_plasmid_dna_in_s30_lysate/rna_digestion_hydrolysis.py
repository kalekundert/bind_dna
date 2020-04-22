#!/usr/bin/env python3

from stepwise import Protocol, MasterMix, load

def hydroxide_hydrolysis():
    rxn = MasterMix("""\
            Reagent         Stock    Volume
            =============  ======  ========
            water                  to 10 µL  
            IVTT reaction            2.5 µL
            EDTA           500 mM   1.25 µL
            NaOH              1 M   1.67 µL
    """)

    return f"""\
Setup an NaOH hydrolysis reaction:

{rxn}

- Add EDTA before NaOH to avoid precipitating Mg²⁺.
- Incubate at 95°C for 10 min.
"""

def carbonate_hydrolysis():
    rxn = MasterMix("""\
            Reagent                    Stock    Volume
            ========================  ======  ========
            water                             to 10 µL  
            IVTT reaction                       2.5 µL
            EDTA                      500 mM   1.25 µL
            sodium carbonate, pH 9.2  100 mM      5 µL
    """)

    return f"""\
Setup an NaHCO₃ hydrolysis reaction:

{rxn}

- Add EDTA before the carbonate buffer, to avoid 
  precipitating Mg²⁺.
- Incubate at 95°C for 10 min.
"""
def rnase_digestion():
    rxn = MasterMix("""\
            Reagent         Stock    Volume
            =============  ======  ========
            water                  to 10 µL  
            IVTT reaction            2.5 µL
            EDTA           500 mM   1.25 µL
            RNase A                  1.0 µL
    """)

    return f"""\
Setup an RNase A digestion:

{rxn}

- Incubate at 60°C for 30 min.
"""

p = Protocol()

p += load('nebex p27 -v 12 -d 75/5')
p += hydroxide_hydrolysis()
p += carbonate_hydrolysis()
p += rnase_digestion()
p += load('gel sdsmax −lysate,−p27,−treatment,NaOH,NaHCO₃,RnaseA -S')
p += load('laser blue')
p += load('gelgreen -r')

print(p)
