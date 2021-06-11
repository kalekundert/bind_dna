#!/usr/bin/env bash
set -euo pipefail

sw cond ligate_unreacted_mrna_cond.xlsx |

sw zap |

sw dilute o237 f11 'f11 (unreacted):2µM' -c '1µM' -v 1 -d 'nuclease-free water' |

sw step "Keep these controls on ice until ready to run the gel." |

# mRNA gene:
# - Initially I was planning to use both FLAG (f111) and Zif268 (f11).
#
# - But I think FLAG will be more informative.
#
#   - My unreacted f111 should be very clean, because the ligated/unligated 
#     bands were very far apart.
#   - In contrast, my unreacted f11 should be relatively dirty, because the 
#     bands were very close.
#   - My unreacted f111 is more concentrated (2 µM) than my unreacted f11 (1.5 
#     µM).
#   - The FLAG/Zif268 ligation reactions have always had comparable yield, so 
#     there's no reason to think that I'll miss anything important if I only 
#     use FLAG.
#
# - Counter-point:
#   
#   - FLAG and Zif268 have different annealing sequences.  That may be 
#     important.
#
# mRNA concentration:
# - [Reyes2021] has final concentrations:
#   - mRNA: 4 µM
#   - linker: 6 µM
#
# - Normally my stocks are 10 µM, so I use final concentrations:
#   - mRNA: 2 µM
#   - linker: 3 µM
#
# - This time my f111 stock is 2 µM, so my final concentrations will be lower.  
#   I'll try to compensate by increasing volume and not using excess linker.
#
# mRNA volume:
# - I'm just going to use a whole 5 µL aliquot.
#
#   - This is really the only experiment I need this RNA for, so there's no 
#     reason to try to conserve too much.
#
#   - Note that 1 µL was already used for the controls above, so only 4 µL 
#     remain for this reaction.
#
#   - I can't add much more than 4 µL mRNA to a 10 µL reaction anyways.  (If I 
#     didn't add any water at all, I could add 6.5 µL 2 µM mRNA and 1.3 µL 10 
#     µM linker.)
#
# Ligase concentration
# - Should I scale the ligase concentration with the mRNA concentration?
# - No:
#   - The ligase is already hard to pipet, decreasing the volume further would 
#     make that worse.
#   - For the hypothesis I'm testing (something about the mRNA inhibits 
#     ligation), there's no harm in having an excess of enzyme.
#
# Reaction volume:
# - 10 µL is more than I need:
#   - I only need 1 µL 1 µM mRNA to run a gel.
#
# - A 10 µL reaction has 0.2 µL ligase.  I'm already worried that I can't pipet 
#   such a small volume accurately.  Decreasing the reaction volume further 
#   would exacerbate the problem.  Some ways around this:
#
#   - Master mix: This won't help a lot, because so much of the volume of the 
#     reaction will be the mRNA.  But it's better than nothing.
#
#   - Echo.  I'd need to buy plates...
#
# Reagent order:
# - Previously I followed [Reyes2021] exactly and added everything but the 
#   ligase to the annealing reaction.
#
#   - [Reyes2021] is the best ligation protocol I have so far, although the 
#     differences between the protocols are small.
#
# - This time I decided to follow my instincts and add as little as possible to 
#   the annealing reaction:
#
#   - Increased RNA/DNA concentration should favor annealing.
#   - Larger volume of ligase master mix should improve pipetting accuracy.
#   - Add ligase buffer to annealing reaction because salt helps annealing.  
#     There isn't actually much salt in this buffer (50 mM Tris, 10 mM MgCl₂), 
#     but anything is better than nothing.  Only add buffer to 1x, so that the 
#     ligase master mix is also in 1x buffer (to keep the ligase happy).
#
# - Probably it'd be better to do a head-to-head comparison of this strategy vs 
#   [Reyes2021], but I'm just gonna cowboy this one.
sw reaction \
  -s "annealing reaction/s" \
  'T4 RNA ligase buffer ;     10x ;   0.67 µL ; +' \
  'f111                 ;    2 µM ;      5 µL ; -' \
  'o237                 ;   10 µM ;      1 µL ; +' \
  -v 5.33 \
  -n 2 |

sw step "Incubate as follows:~90°C for 30s~Cool to 25°C at 1°C/s" |

# Extra master mix:
# - All the master mix components are cheap, so I can make lots of extra.
# - I think I can pipet 0.5 µL ligase accurately.
sw reaction \
  -s "ligation reaction/s" \
  'nuclease-free water  ;         ;  to 10 µL ; +' \
  'T4 RNA ligase buffer ;     10x ;   0.33 µL ; +' \
  'ATP                  ;   10 mM ;      1 µL ; +' \
  'T4 RNA ligase        ; 10 U/µL ;    0.2 µL ; +' \
  'annealed mRNA        ;  1.5 µM ;   6.67 µL ; -' \
  -v 3 \
  -X 0.5 \
  -n 2 |

sw reaction \
  -s "-ligase control/s" \
  'nuclease-free water  ;         ;  to 10 µL ; +' \
  'T4 RNA ligase buffer ;     10x ;   0.33 µL ; +' \
  'ATP                  ;   10 mM ;      1 µL ; +' \
  'annealed mRNA        ;  1.5 µM ;   6.67 µL ; -' \
  -v 3 \
  -X 0.5 \
  -n 2 |

# Incubation time:
# - From [Reyes2021].  In my experience (expt #3), the reaction is complete 
#   within 10 min.
sw step "Incubate at 25°C for 30 min." |
sw step "Label the products: f113" |

sw gel urea/mrna/flag 7 -s gelgreen_cy5

