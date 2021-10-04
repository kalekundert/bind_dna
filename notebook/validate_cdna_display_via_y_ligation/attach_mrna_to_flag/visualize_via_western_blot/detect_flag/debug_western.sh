#!/usr/bin/env bash
set -euo pipefail

ABCAM_PROTOCOL=https://tinyurl.com/t2n6tr63

# FLAG concentration:
# - In :expt:`126`, I used 2 ng FLAG.  This is roughly the molar amount of a 
#   FLAG-labeled protein that Abcam claims can be detected by their primary 
#   (I'm using a different primary, but I couldn't find any such 
#   recommendations for my antibody.)
#
# - Here I want to use much more.  I'm going to go for 1 µg, which is the most 
#   I can do with 1 µL of my stock.  This should be detectable without any 
#   doubts.

sw tbst |

sw blocking_buffer.txt |

#sw gel bolt/mes "FLAG peptide (r5)" -c 1000 -v 1 -r 30 -S -L "5 µL Chameleon 800" |
sw step "Load 10 µL of FLAG-fusion positive control." |
sw gel bolt/mes -M -S -L "5 µL Chameleon 800" |

sw transfer_iblot.txt |

sw step "Trim off unused parts of the membrane and mark the side facing the gel." |

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
# reference.  It also makes sense: how will casein bind the membrane if it's 
# surrounded by detergent?  I'll just block with 5% milk.

sw step "Incubate the membrane in 20 mL blocking buffer for 1h at room 
temperature with gentle shaking." |
sw note $ABCAM_PROTOCOL |

# Primary antibody concentrations:
# - FUJIFILM recommends 1,000x-10,000x dilution.
# - [Reyes2021]_ used 12,000x dilution.
# - I'm going to use 1000x until I detect some signal.

sw step "Prepare the primary antibody:
~5 mL TBST
~5 mL blocking buffer
~10 µL monoclonal mouse anti-DYKDDDDK tag (r4)" |

sw step "Incubate the membrane in the primary antibody solution overnight at 
4°C, in the dark and with gentle shaking." |
sw note $ABCAM_PROTOCOL |

sw step "Store unused blocking buffer at 4°C overnight." |
sw note $ABCAM_PROTOCOL |

# Wash steps:
# - I want to do less washing, because lack of signal could be caused by 
#   overwashing.
#
# - Fitzy says that she basically just rinses a few times, maybe with a 10 min 
#   wash.
#
# - This thread has several recommendations in the 2-5x 5-10 min range:
#   http://www.protocol-online.org/biology-forums-2/posts/21073.html#
#
#   Most of the recommendations call for 20-30 min of total washing.
#
# - I'm going to try 3x 5 min washes.  That's on the low end of what I've seen 
#   recommended, but it's not unreasonable.
sw step "Wash the membrane as follows:
~Wash 3x for 5 min in TBST (with shaking)" |
sw note $ABCAM_PROTOCOL |

sw step "Prepare a 1000x dilution of the secondary antibody.  Keep the antibody 
in the dark as much as possible:
~5 mL TBST
~5 mL blocking buffer
~10 µL goat anti-mouse IRDye 800CW (r3)" |
sw note $ABCAM_PROTOCOL |

sw step "Incubate the membrane in the secondary antibody solution at room 
temperature for 1h, with shaking." |
sw note $ABCAM_PROTOCOL |

sw step "Repeat the above washing procedure." |
sw note $ABCAM_PROTOCOL |

sw laser nir -i "Place the membrane in the scanner protein-side down."













