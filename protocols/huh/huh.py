#!/usr/bin/env python3

"""\
Attach a DNA sequence to a protein fused to a Rep domain.

Usage:
    huh <protein> <dna> [options]

Arguments:
    <protein>
        The names of 1 or more proteins to attach to the below DNA sequence(s).
    <dna>
        The names of 1 or more DNA sequences to attach to the above protein(s).  These names 
        should be present in the PO₄ database and have associated stock 
        concentrations and molecular weights.  Note that PO₄ can't 
        automatically calculate molecular weights for non-standard nucleotides 
        like iSP9, so you may need to provide a molecular weight manually for 
        such constructs.

Options:
    -n --num-reactions <int>
        The number of reactions to set up.  By default, this is the number of 
        DNA sequences given.

    -b --buffer <name>              [default: Vega-Rocha Rep buffer]
        Which buffer to use.  The default is the buffer from [VegaRocha2007].  
        The buffer is assumed to be 10x.

    -B --buffer-divalent-conc <mM>  [default: 2.5]
        The concentration of divalent cation in the buffer, in mM.  This is 
        used to calculate the appropriate concentration of EDTA to add to the 
        negative control.

    -c --cas9-stock-conc <µM>       [default: 22]
        The stock concentration of dCas9-PCV2, in µM.

    -S --no-sgrna
        Exclude sgRNA from the reaction.

    -m --master-mix <reagents>      [default: cas9,sgrna,buffer]
        A comma-separated list of reagents to include in the master mix.  The 
        following reagents are understood: cas9, sgrna, dna, edta, buffer
"""

import docopt
import stepwise
import po4

# This script should be broken in two: one for proteins in general and another 
# for Cas9.  The only real difference for Cas9 is that the reaction includes 
# sgRNA and I know some concentrations up front.

args = docopt.docopt(__doc__)
p = stepwise.Protocol()

huh = stepwise.MasterMix.from_text("""\
Reagent      Stock    Volume  MM?
==========  ======  ========  ===
water               to 10 µL  yes
buffer         10x      1 µL  yes
dCas9-PCV2    1 µM      1 µL  yes
sgRNA         1 µM      1 µL  yes
EDTA        500 mM      1 µL
DNA         200 nM      5 µL
""")

db = po4.load_db()

dna = args['<dna>']
dna_nM = min([round(db[x].conc_nM) for x in dna])

huh.num_reactions = int(args['--num-reactions'] or len(dna))
huh.extra_min_volume = '0.5 µL'

huh['buffer'].name = args['--buffer']
huh['dCas9-PCV2'].hold_conc.stock_conc = args['--cas9-stock-conc'], 'µM'
huh['sgRNA'].hold_conc.stock_conc = 10, 'µM'
huh['EDTA'].hold_stock_conc.conc = 30 / 2.5 * float(args['--buffer-divalent-conc']), 'mM'
huh['DNA'].hold_conc.stock_conc = dna_nM, 'nM'
huh['DNA'].name = ','.join(dna)

mm = args['--master-mix'].split(',')
huh['dCas9-PCV2'].master_mix = 'cas9' in mm
huh['sgRNA'].master_mix = 'sgrna' in mm
huh['DNA'].master_mix = 'dna' in mm
huh['EDTA'].master_mix = 'edta' in mm
huh['buffer'].master_mix = 'buffer' in mm

if args['--no-sgrna']:
    del huh['sgRNA']
else:
    p += stepwise.load('zap')
    huh['water'].name = 'nuclease-free water'

p += f"""\
Attach DNA to HUH-tagged dCas9 [1-3]:

{huh}

- Add each reagent in order.
- Mix after adding EDTA (and before adding DNA).
- Incubate at 37°C for 15 min.
"""

# 500 ng/band
p.footnotes[1] = """\
Invitrogen recommends loading no more than 250 
ng/band on Bolt SDS PAGE gels, and the detection 
limit for Coomassie (as an IR dye [Butt2013]) is 
at least 10 ng/band.  That corresponds to a range 
of ≈1.4-0.1 pmol Cas9-PCV2/band (MW: 176.8 kDa).  
I probably want to be on the high end of that.
"""

p.footnotes[2] = """\
The EDTA reaction is a negative control 
established in [VegaRocha2007].  They used 2.5 mM 
divalent metal and 30 mM EDTA to prevent coupling.  
"""

p.footnotes[3] = """\
f16 and f12 are 414 bp.  At that length, 50 ng/µL 
(a typical PCR yield) corresponds to ≈150 nM.  
[Lovendahl2017] used a 10:1 DNA:protein ratio to
maximize the amount of coupled protein.  I'm going 
to use a 1:1 ratio instead, both because I don't 
want a lot of unbound DNA in my qPCR reactions and 
because a 10:1 ratio would use a lot of material.
"""

print(p)
