#!/usr/bin/env bash
set -euo pipefail


## Strains:
## - p224 TGG/AAA: single-plasmid construct that works the best.
## - Can't think of reasons to use any other constructs.
##
## Volume:
## - Need enough for 5 purifications.
## - Direct-zol and Quick-RNA kits call for 1e8 bacterial cells.
## - TRIzol calls for 1e7 cells.
## - Rule of thumb is that "each 0.1 OD unit is roughly equivalent to 1e8 
##   cells/ml" [Elbing2018, DOI: 10.1002/cpmb.81]
## - Overnight culture are roughly OD 1.0, i.e. 1e9 cells/mL [Elbing2018].
## - I don't want saturated culture, so I'll aim for OD 0.1 and use 1 mL cells.
##
## Flask vs tubes?
## - For cultures larger than 5 mL, I usually prefer to use flasks.
## - But 6 mL is close to 5 mL, and I want to maintain my usual culture 
##   conditions as much as possible.  So I'll stick with tubes.
#sw grow sz224 sz228 -v 6 -d 0.1 -t 4 |

## Extract with:
## - TRIzol
## - Direct-zol - DNase
## - Direct-zol + DNase
## - Quick-RNA miniprep - DNase
## - Quick-RNA miniprep + DNase
##
## Was on the fence about whether to buy Direct-zol and Quick-RNA, since I'll 
## probably only end up using one or the other, but they're cheap and I might as 
## well find out if one works better than the other.
#sw step "Purify each sample using TRIzol, Direct-zol (+/−DNase), and Quick-RNA miniprep (+/−DNase), according to manufacturers' instructions." |
#sw trizol |

#sw reverse_transcribe |
#sw sub "reverse transcription reactions" "reverse transcription reactions (including −RT controls)" |

# qPCR
# - no template
# - DNA template
# - For each sample:
#   - no RT 
#   - +RT
sw step "Prepare −template and +template controls." |

# Number of reactions:
# - 80 total reaction (ignoring +/− template controls)
#   - 2 primers
#   - 5 extractions
#   - +/− RT
#   - 2 replicates
#   - 2 templates
#
# - Separate master mixes for each primer.
# - Duplicates option accounts for replicates.
# - That means N=5*2*2=20 (extraction, RT, template) for each master mix.
sw pcr -p ssoadv \
  sz224,sr1,sr2 \
  sz228,sr1,sr2 \
  -n 2*2*5+2 \
  -d 2 \
  --skip-thermocycler |

sw sub 228 '228 (5 extractions, ±RT)' |

sw pcr -p ssoadv \
  sz224,sr3,sr5 \
  sz228,sr3,sr5 \
  -n 2*2*5+2 \
  -d 2 |

sw sub 228 '228 (5 extractions, ±RT)'

