#!/usr/bin/env bash
set -euo pipefail

# Binding buffer:
# - Storage buffer (expt 121) + BSA + ssDNA
#   
#   - 10 mM tris
#   - 50 mM NaCl
#   - 1 mM MgCl₂
#   - 90 µM ZnCl₂
#   - 5 mM DTT
#   - 200 µg/mL BSA
#   - 80 µg/mL ssDNA
#   - pH=7.5
#
# - Why not tween?
#
#   - My previous EMSA reactions (buffer from [Lam2011]) used tween.
#   - The [Zykovich2009] buffer doesn't have detergent, but many others do.
#   - Leaving it out just to follow [Zykovich2009], no really good reason.

sw binding_buffer.txt |
sw step "Prepare 80 µL 1x buffer:
  ~64 µL water
  ~18 µL 5x Zif268 binding buffer" |

# Steps:
# - Want to fit everything on a single gel, i.e. 15 lanes.
# - That means 7 lanes per target.
sw cond confirm_zif268_activity_cond.xlsx |

# DNA concentration:
# - 2019/08/28:
#   - 625 pM in binding reaction: 2 µL × 7.5 nM / 24 µL
#   - Loaded 10 µL of 32 µL sample: 4.68 pmol target DNA
#   - Got good signal
# - [Hellman2007]
#   - 18.75 fmol target
#
# - The standard native PAGE protocol calls for 1.25 µL sample brought to a 
#   total volume of 5 µL before loading.
#   -  4.68 fmol / 1.5 µL = 3.74 nM in binding reaction
#   - 18.75 fmol / 1.5 µL = 15.0 nM in binding reaction
#
#   - This seems like a reasonable range.
#   - I'm inclined to pick a round number in the middle of the range, so I'm 
#     going to go with 10 nM.
#
# - 2021/08/17: The signal with 10 nM DNA has been too faint, so I'm going to 
#   double the volume of DNA being added to the reaction to reach 20 nM.  
#   That's probably still not enough to get good signal, but it should be a 
#   step in the right direction.  I'm not going to increase the protein 
#   concentration accordingly; see below.
#
# Protein concentration:
# - [Hellman2007] recommends titrating against constant DNA.  Seems reasonable.
# - [Hellman2007] concentrations: 0.75x to 21x in 9 (unevenly spaced) steps
# - I'm doing 0.5x to 16x in 6 (evenly spaced) steps.  Seems comparable.
#
# - 2021/08/17: I want to have more lanes with <1x protein, and I also want to 
#   increase the quantity of DNA in the reaction.  So I'm going to double the 
#   quantity of DNA while holding the quantity of protein constant, meaning 
#   that my titration will become 8x to 0.25x.
#
# Protein diluent:
# - The protein is already in storage buffer, so it makes sense to use that for 
#   the diluent.
# - I wouldn't want to use binding buffer, because then different titrations 
#   would get slightly different amounts of BSA, ssDNA, etc.

sw serial 10µL 50nM x 2 6 -0 -m "PCV2-Zif268" -d "1x Zif268 storage buffer" |
sw reaction \
  -s binding \
  "water;                ; to 10 µL; +" \
  "binding buffer;     5x;     2 µL; +" \
  "f128,f129;      100 nM;     2 µL; +" \
  "PCV2-Zif268;       10x;     1 µL; -" \
  -n 7 \
  --extra-reactions 1 \
  -i "Make a separate master mix for each target." |

# Incubation time:
# - Needs to be long enough to achieve equilibrium [Hellman2007], or results 
#   will not be reproducible.
#   
#   - 30 min is usually enough.
#   - A simple control: Setup 30- and 60-min reactions, and check that they 
#     appear the same.
#
# - I did 1h on 2019/8/28.  I might just do that again.

sw step "Incubate at room temperature for 1h." |

sw gel emsa/zif 14




  
