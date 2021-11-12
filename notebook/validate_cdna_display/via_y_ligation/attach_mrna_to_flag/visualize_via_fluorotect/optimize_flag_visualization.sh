#!/usr/bin/env bash
set -euo pipefail

# Quick experiment to make sure I can visualize FLAG expression.
# - I'm going to wait to do this experiment until I can make the correct 
#   loading dye.

# Would I rather test load volume of loading buffer?
#
# - Loading less material will reduce signal and background.  The background is 
#   high in #99, but the signal is also low.  Note that I also had to use 
#   maximum sensitivity in the LysGreen channel for that experiment.  I'm 
#   skeptical that reducing the loading volume will make any difference.
#
# - The loading buffer formulations should be identical other than the dye, so 
#   I don't really have any reason to expect a difference.  
#
# - Gel type
#

# Load volumes:
# - Want to try two volumes:
#   - 2.5 µL: What I normally do for these reactions.
#   - 1.0 µL: What the PUREfrex manual recommends with FluoroTect.
#
# Conditions:
# - Load volume:
#   - I normally load 2.5 µL of each reaction.
#   - The PUREfrex manual recommends loading 1 µL with FluoroTect.
#   - But, my protein is very small.
#   - A side-by-side comparison won't be informative.  I know that adding more 
#     will increase signal and background.  I can tweak the sample volume as I 
#     do other experiments.
#   - In general, I think I would prefer 2.5 µL
#     - 1 µL seemed low if anything in #99
#       - The unreacted band is bright, but the protein band is dim.  If I 
#         eliminate the unreacted band with RNase A, I'll probably want the 
#         peptide band to be brighter if anything.
#   - I might use 1 µL for this experiment though, just to stretch my sample 
#     and maybe to get a sense for visualizing small qauntities.
#
# - Gel type:
#   - The FLAG peptide is very small (2.2 kDa).
#   - Proteins that small basically run at the dye front in Bolt gels.
#   - Tris-tricine gels are optimized for peptides.
#   - Two options:
#     - 16%: Best resolution of any gel around 2 kDa.
#     - 10-20%: Not as good resolution as 16% around 2 kDa, but overall range 
#       comparable to Bolt 4-16%.  Unsure if this would differ from Bolt 
#       substantially.
#   - I should try both.  I don't know what they'll look like, and I might want 
#     to use both going forward.
#   - Include Bolt as a control, too.
#   - Crystal violet loading dye for all.
sw cond optimize_flag_visualization_cond.xlsx \
  -k 'A: Bolt Bis-Tris-MES 4-16%' \
  -k 'B: Novex Tris-Tricine 10-20%' \
  -k 'C: Novex Tris-Tricine 16%' |

# Template:
# - I definitely want to use f110/f111 for this experiment, because the FLAG 
#   peptide is especially tough to see given how small it is.
#
# Template concentration:
# - Based on #99, 2 nM final concentration for the DNA template seems too low.
# - I could use 800 nM of the RNA template instead.  That seems to be the best  
#   based on what I know so far.
# - I could also use PURExpress...  I have more of it, but I don't have an 
#   optimized concentration... I'm on the fence.
sw ivtt f111 -p frex1/lys -n 2 -v 10 -c 800 -C 10000 -r |
#sw ivtt f110 -p purex/lys -n 2 -v 10 -c 2 -C 75 |

# RNase A volume:
# - The PUREfrex 2.0 manual calls for 1 µL 1 mg/mL RNase A per 10 µL reaction.  
# - I don't have any RNase A on hand, but I do have RNase cocktail (Invitrogen 
#   AM2286).  This cocktail contains:
#   - 500 U/mL RNase A
#   - 20,000 U/mL RNase T1
# - The manual for the cocktail says to replace RNase A at equivalent RNase A 
#   concentrations.
# - I can't easily find a way to convert between mass and units, so I'm just 
#   going to use the same proportion as in the PUREfrex manual and hope that 
#   works.
sw step "Split each reaction into 2x 4 µL aliquots." |
sw step "Add 0.4 µL RNase cocktail (Invitrogen AM2286) to the +RNase A aliquots." |
sw step "Incubate at 37°C for 15 min." |

sw gel bolt/ivtt/flag 4 -v 1 -S |
sw gel tricine/ivtt/flag 4 -v 1 -p '10-20' -S |
sw gel tricine/ivtt/flag 4 -v 1 -p '16'


