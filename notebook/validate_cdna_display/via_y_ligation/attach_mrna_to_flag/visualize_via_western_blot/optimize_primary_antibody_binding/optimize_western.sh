#!/usr/bin/env bash
set -euo pipefail

ABCAM_PROTOCOL=https://tinyurl.com/t2n6tr63

# FLAG quantity:
#
# - Want to begin with an amount that will be easily detectable.
#
# - Detection limit depends on primary antibody, but manufacturer doesn't make 
#   any suggestions.
#
# - Abcam notes that one of their polyclonal anti-FLAG antibodies detects 25 ng 
#   and 100 ng of positope, which is 53 kDa protein containing numerous epitope 
#   tags.  For comparison, the FLAG peptide is 1 kDa.  
#
#   My antibody is obviously different, but this should be a reasonable 
#   starting point.  50x less than 100 ng is 2 ng.
#
#   https://www.abcam.com/ddddk-tag-binds-to-flag-tag-sequence-antibody-ab1162.html

# FLAG concentration:
# - If I want to load 2 ng total, 1 ng/µL is a reasonable concentration.
# - I'm ordering 500 µg of peptide.
# - That works out to: 250,000 µL
# - Ok, I'll need a more conc stock

sw gel bolt/mes "FLAG peptide" -n 6 -c 2 -v 1 -S -L "5 µL Chameleon 800" |

sw transfer_iblot.txt |

sw tbst |

sw step "Visualize the protein with Ponceau staining:
~ Soak membrane in Ponceau stain for 15 min.
~ Rinse membrane in water until bands are visible.
~ Cut each band into its own strip.
~ Label the side of the membrane with protein.
~ Repeat until stain is no longer visible (4-5x):
~~ Rinse for 5-10 min in 1x TBST." |

# Blocking buffer:
# - Abcam's general western blotting protocol calls for Tween-20 in the 
#   blocking buffer.
#   https://www.abcam.com/protocols/general-western-blot-protocol#Solutions%20and%20reagents%20running,%20transfer%20and%20blocking 
#
# - Abcam's fluorescent western blotting protocol douesn't really say what the 
#   blocking buffer should be:
#   https://www.abcam.com/secondary-antibodies/fluorescent-western-blot-protocol--irdye-secondary-antibodies
#
# - LiCor's fluorescent western guide recommends not adding detergent to 
#   the blocking buffer:
#   https://www.licor.com/documents/dq6jb8sgnkiwlas0g99b5hq0tsuvcyyz
#
# I'm inclined to trust LiCor, since it's the most clear and specific 
# reference.  I'll just block with 5% milk.

sw blocking_buffer.txt |
sw step "Incubate each membrane strip in blocking buffer for 1h at room 
temperature with gentle shaking." |
sw note $ABCAM_PROTOCOL |

# Primary antibody concentrations:
# - FUJIFILM recommends 1,000x-10,000x dilution.
#
# - [Reyes2021]_ used 12,000x dilution.
#
# - 3-step serial dilutions:
#   - 1000x, 3464x, 12000x
#   - 1000x, 3162x, 10000x
#
# - I think it makes sense to include the same conditions as [Reyes2021]_.  And 
#   it probably doesn't really matter.
#
# - I'll probably round the middle dilution to 4000x and do direct dilutions.

sw step "Prepare 3 dilutions of the primary antibody:
~12 mL TBST
~12 mL blocking buffer
~24 µL, 6 µL, 2 µL monoclonal mouse anti-DYKDDDDK tag" |

sw step "Place each membrane strip in 10 mL primary antibody dilution." |
sw note $ABCAM_PROTOCOL |

sw step "Incubate 3 of the strips (one for each dilution) at room temperature for 1h, with gentle shaking.
~Continue visualizing these strips immediately after the blocking is complete; don't wait for the overnight incubation to finish." |

sw step "Incubate the remaining strips at 4°C overnight, with gentle shaking." |
sw note $ABCAM_PROTOCOL |

sw step "Discard the primary antibody solution." |
sw note $ABCAM_PROTOCOL |

sw step "Wash the membrane as follows:
~Rinse 2x with TBST
~Wash 1x for 15 min in TBST (with shaking)
~Wash 3x for 10 min in TBST (with shaking)" |
sw note $ABCAM_PROTOCOL |
sw note "The specifics of the washes don't seem to be important.  I've seen other protocols that call for 3x rinses and 5x 5-10 min washes." |

sw step "Prepare a 10,000x dilution of the secondary antibody.  Keep the antibody in the dark as much as possible:
~17 mL TBST
~17 mL blocking buffer
~3.4 µL goat anti-mouse IRDye 800CW" |
sw note $ABCAM_PROTOCOL |

sw step "Incubate each of the strips in 10 mL secondary antibody dilution at room temperature for 1h." |
sw note $ABCAM_PROTOCOL |

sw step "Repeat the above washing procedure." |
sw note $ABCAM_PROTOCOL |

sw laser nir -i "Place the membrane in the scanner protein-side down." |

sw step "Dry the membrane strips overnight between two sheets of filter paper, covered by aluminum foil." |

sw laser nir -i "Place the membrane in the scanner protein-side down."













