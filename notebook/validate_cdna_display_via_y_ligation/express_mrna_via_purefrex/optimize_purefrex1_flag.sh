#!/usr/bin/env bash
set -euo pipefail

sw zap |
sw cond optimize_purefrex1_flag_cond.xlsx |

# Volume:
# - f111 aliquots are 5 µL, so I tried to increase the volume as much 
#   as possible (to make pipetting more accurate) without using more than one 
#   aliquot.
sw serial 2.6 '10000 nM' to 100 6 -0 -m f111 -d 'nuclease-free water' |

# mRNA concentration:
# - 1 µM was the highest concentration I tested for PUREfrex 2.0 with mWasabi.   
#   That didn't saturate expression, though, so I want to try going higher to 
#   see if I can find the peak.
# - I'm going to make 20 nM the lowest concentration (100 nM template stock), 
#   because I saw very little expression below that level with PUREfrex 2.0.  

sw step "Dilute f110 to 10 nM so that it will have a final concentration of 2 nM in the following reaction." |

sw ivtt f111 -p frex1/lys -v 2.5 -n 8 -c 2000 -C 10000  |

sw step "Add 0.5 µL 1 mg/mL RNase A to each reaction." |
sw step "Incubate at 37°C for 15 min." |

# Run time:
# - When I've run gels for 42m in previous experiments, the 3 kb band has been 
#   right at or just past the bottom.
# - My FLAG peptide is 2.2 kDa (assuming it is translated all the way to the 
#   end of the transcript; there is no stop codon).
# - 42m is probably too long; I'll try 30m to start.
sw gel bolt/ivtt 8 -v 1 -r 30

