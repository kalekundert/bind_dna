#!/usr/bin/env python3

"""\
Usage:
    prepare_dna_template_via_digest.py <templates>... [options]

Arguments:
    <templates>
        The names of the templates to prepare.

Options:
    -d --dna-conc NG_PER_UL  [default: 100]
        The concentration of the DNA to prepare (in ng/µL).

    -D --dna-ng NANOGRAMS  [default: 5000]
        How much DNA to prepare (in ng).
"""

import docopt
import dirty_water

# This script is a poster-child for why my system for writing protocols sucks.  
# I'll outline some of the problems below.

args = docopt.docopt(__doc__)
templates = args['<templates>']
dna_ng = float(args['--dna-ng'])
dna_ng_uL = float(args['--dna-conc'])

# Work out how much to add of each reagent.

# This is way too much code.  My library should handle these basic kinds of 
# volume manipulations.  The actual volume of the reaction is kinda hidden, and 
# the volume used up by the non-template reagents is calculated twice.  Some of 
# this is because I want water to be the first reagent.  I should be able to 
# specify the order reagents are added, independently of how they appear in the 
# source.

xmni_U = 1 * (dna_ng / 1000)            # 1U per 1 µg DNA, see ./xmni_sites_in_lambda.py
xmni_U_uL = 20                          # NEB R0194{S,L}: 20 U/µL
xmni_uL = max(xmni_U / xmni_U_uL, 0.5)  # Don't use less than 0.5 µL.

max_dna_uL = 50 - 5 - xmni_uL
dna_uL = min(dna_ng / dna_ng_uL, max_dna_uL)
dna_ng = dna_uL * dna_ng_uL

digest = dirty_water.Reaction()
digest.num_reactions = len(templates)

if dna_uL < max_dna_uL:
    digest['Water'].std_volume = max_dna_uL - dna_uL, 'µL'
    digest['Water'].master_mix = True

template_names = ','.join(templates)
digest[template_names].std_volume = dna_uL, 'µL'
digest[template_names].std_stock_conc = dna_ng_uL, 'ng/µL'

digest['CutSmart buffer'].std_volume = 5, 'µL'
digest['CutSmart buffer'].std_stock_conc = '10x'
digest['CutSmart buffer'].master_mix = True

# There should be a first-class way to include product numbers in these 
# protocols.
digest['XmnI (NEB R0194)'].std_volume = xmni_uL, 'µL'
digest['XmnI (NEB R0194)'].std_stock_conc = 20, 'U/µL'
digest['XmnI (NEB R0194)'].master_mix = True

# Print the protocol.

protocol = dirty_water.Protocol()

protocol += """\
Digest the template plasmid(s) with XmnI:

{digest}

- Incubate 37°C for 1h, then at 65°C for 20 min."""

# This is copied-and-pasted from another protocol, with only the bulleting 
# changed.  How can I avoid this duplication?  Part of the problem is that the 
# two protocols live in different places.  Also, it'd be nice to have an easy 
# way to hide the notes.
protocol += """\
Do a phenol-chloroform extraction to remove any  
RNase leftover from the plasmid prep.

- Dilute the reaction to 500 µL with water.
   
- Add 500 μL phenol:chloroform:isoamyl alcohol, pH 
  8.0 (i.e. the correct pH for purifying DNA).

- Vortex vigorously to mix the phases.

- Spin in a microfuge at top speed for 1-2 min to 
  separate the phases.

- Transfer the aqueous (upper) phase to a new 
  tube, being careful not to transfer any of the 
  protein at the phase interface.

- Repeat the phenol:chloroform:isoamyl alcohol 
  extraction.

- Extract the sample with an equal volume of 
  chloroform:isoamyl alcohol to remove phenol.
"""

# My existing EtOH precipitation protocol is much more verbose that this, but 
# I'd again like to be more DRY.
protocol += f"""\
Precipitate and resuspend the DNA.

- Add to the purified DNA:

  - 2 volumes 100% EtOH
  - 0.1 volumes 3M NaOAc, pH 5.2
  - 1 µL glycogen (optional)

- Incubate at -80°C for 20-60 min.

- Spin max speed, 15 min, 4°C.

- Remove supernatant and rinse pellet with 500 µL 
  70% ethanol.
  
- Spin max speed, 15 min, 4°C.

- Remove supernatant and air dry the pellet for 15 
  min at room temperature.

- Resuspend the pellet in {dna_ng/500} µL nuclease-free
  water (for a concentration of 500 ng/µL)."""

# References:
# https://www.neb.com/protocols/0001/01/01/dna-template-preparation-e2040

# I'd like a way to make this information 
# available in the actual protocol, like with a 
# verbose option or something.

print(protocol)
