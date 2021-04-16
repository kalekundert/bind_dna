#!/usr/bin/env bash
set -euo pipefail

sw zap |

# Volume:
# - f15 aliquots are 5 µL, so I tried to increase the volume as much 
#   as possible (to make pipetting more accurate) without using more than one 
#   aliquot.
#
# - I haven't made f15 yet, so I don't know it's concentrations, but I'm going 
#   to aim for 10 µM.  This has been hard to reach in the past, but I'm using 
#   bigger spin columns now.  I might also scale up the reaction a bit.
sw serial 2.6 '10000 nM' to 100 7 -0 -m f15 -d 'nuclease-free water' |

sw ivtt f15 -p purefrex1 -v 2.5 -n 8 -c 2000 -C 10000  |

sw gel bolt/ivtt 8 

