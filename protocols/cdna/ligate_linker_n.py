#!/usr/bin/env python3

"""\
Ligate linker-N to the mRNA.

Usage:
    ligate <n> [-v <µL>] [-x <percent>] [-i <time>] [-m <reagents>] [options]

Arguments:
    <n>
        The number of reactions to set up.

Options:
    -v --volume <µL>                [default: 40]
        The volume of each ligation reaction in µL.

    -x --extra <percent>            [default: 10]
        How much extra master mix to prepare.

    -i --incubate <time>            [default: 10 min]
        How long to incubate the reaction at the temperature indicated by the 
        `-I` flag.  Include a unit.

    -I --incubate-temp <temp>       [default: 25°C]
        What temperature to incubate the reaction at.  Include a unit.

    -m --master-mix <reagents>      [default: peg,lig]
        Include the indicated reagents in the master mix.  The following 
        reagents are understood:

        pnk: T4 PNK
        lig: T4 RNA ligase
        rna: annealed mRNA/linker
        peg: PEG-6000

    -M --no-master-mix
        Exclude all optional reagents from the master mix, i.e. `-m ''`.

    -k --pnk
        Add T4 PNK to the reaction, e.g. if using non-phosphorylated primers.

    -p --peg
        Include PEG-6000 in the reaction.  Many T4 RNA ligase protocols 
        recommend this, but in my hands it does not improve yield.

    -Q --no-quench
        Leave out the 65°C incubation to quench the reaction.  This is useful 
        if the reaction will be quenched by a downstream step anyways.
"""

import stepwise, docopt
from inform import plural

# I can't remember where exactly this protocol came from.  Some thoughts:
# 
# - 10x diluted PBS: my own optimizations of the annealing reaction.
# 
# - 20 U T4 RNA ligase: Probably [Naimudden2016]_, even though I'm using 10x 
#   less mRNA/linker (I probably assumed that more enzyme wouldn't hurt, and 
#   this is already about as little as I can pipet).
# 
# - BSA: Probably because Takara included BSA with the enzyme, although I'm 
#   using a different concentration that their online protocol suggests.
# 
# - 10 min incubation time at 25°C: [Naimudden2016]_.  I don't know where the 
#   65°C incubation/heat denaturation came from, though.

args = docopt.docopt(__doc__)
master_mix = '' if args['--no-master-mix'] else args['--master-mix']

ligate = stepwise.MasterMix.from_text("""\
Reagent                  Stock       Volume  MM?
=====================  =======  ===========  ===
nuclease-free water             to 40.00 µL  yes
T4 DNA ligase buffer       10x       4.0 µL  yes
BSA                       0.1%       4.0 µL  yes
PEG 6000                   50%      20.0 µL  yes
T4 PNK                 10 U/µL      0.33 µL
T4 RNA ligase          40 U/µL       0.5 µL
annealed mRNA/linker   1.25 µM       4.0 µL
""")

ligate.num_reactions = n = eval(args['<n>'])
ligate.extra_percent = eval(args['--extra'])
ligate.hold_ratios.volume = eval(args['--volume']), 'µL'
ligate['T4 PNK'].master_mix = 'pnk' in master_mix
ligate['T4 RNA ligase'].master_mix = 'lig' in master_mix
ligate['annealed mRNA/linker'].master_mix = 'rna' in master_mix
ligate['PEG 6000'].master_mix = 'peg' in master_mix

if not args['--pnk']:
    del ligate['T4 PNK']
if not args['--peg']:
    del ligate['PEG 6000']

protocol = stepwise.Protocol()

protocol += f"""\
Setup {plural(n):# ligation reaction/s}:

{ligate}
"""

protocol += f"""\
Incubate the {plural(n):ligation reaction/s} as follows:

- {args['--incubate-temp']} for {args['--incubate']}.
{'' if args['--no-quench'] else '- 65°C for 10 min.'}
"""

print(protocol)
