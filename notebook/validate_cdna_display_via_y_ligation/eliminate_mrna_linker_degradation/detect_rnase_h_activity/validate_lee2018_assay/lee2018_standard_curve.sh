#!/usr/bin/env bash
set -euo pipefail

sw cond lee2018_standard_curve_cond.xlsx |
sw note "https://doi.org/10.1007/s00604-017-2610-8" |

# It would make more sense to me if the ligation reactions were done in a lower 
# volume.  This would maximize the concentrations of the oligos for the 
# annealing step, minimize the amount ligase buffer that ends up in the final 
# reaction, and allow the ligase to be pipeted more accurately.  But the 
# [Lee2015] protocol is very detailed, so I'll follow it exactly.
sw reaction \
  -t "template master mix/es" \
  "nuclease-free water                   ;           ; to 17 µL ; +" \
  "E. coli DNA ligase buffer (NEB B0205) ;       10x ;     2 µL ; +" \
  "o226                                  ;      1 µM ;     2 µL ; +" \
  -n 10 \
  -V 1.1 \
  |

sw reaction \
  -t "control ligation/s" \
  "template master mix                   ;           ;    17 µL ; +" \
  "o227,o228                             ;      2 µM ;     2 µL ; -" \
  -n 2 \
  |

sw reaction \
  -t "standard curve ligation/s" \
  "template master mix                   ;           ;    17 µL ; +" \
  "o229                                  ;      2 µM ;     2 µL ; +" \
  -n 8 \
  |

sw step "Run the following thermocycler protocol:~80°C for 5 min~Cool to 16°C at 0.1°C/sec.~16°C for 15 min." |

sw step "Add 1 µL E. coli DNA ligase (NEB M0205, 10 U/µL) to each reaction." |

sw step "Incubate at 16°C for 1h." |

sw note "Store at 4°C until ready for use." |

sw step "Repeat the following dilution twice to make a 10 U/mL RNase H stock solution:~21.36 µL 1x RNase H buffer (NEB B0297)~1 µL RNase H (NEB M0297)" |

# Much more than I need; to help pipet accurately and to facilitate 
# multichannel pipetting.
sw serial 10 10U/mL to 0.1 7 -0 -m "RNase H" -d "1x RNase H buffer"|

# Volumes adjusted to match stock concentrations I have on hand:
# - dNTPs
# - BSA
sw reaction \
  -t "RCA" \
  "nuclease-free water                   ;           ; to 40 µL ; +" \
  "dNTP mix (Takara 3261)                ;     10 mM ;     1 µL ; +" \
  "SYBR Green II (Invitrogen S7564)      ;       20x ;     2 µL ; +" \
  "BSA (NEB B9000)                       ;  20 mg/mL ;   0.2 µL ; +" \
  "Φ29 buffer (NEB B0269)                ;       10x ;     4 µL ; +" \
  "Φ29 DNA polymerase (NEB M0269S)       ;   10 U/µL ;     1 µL ; +" \
  "ligation reaction                     ;           ;    20 µL ; -" \
  "RNase H dilution                      ; 0-10 U/mL ;     4 µL ; -" \
  -n 10 \
  -i "Use white tubes." \
  |

sw step "Incubate at 30°C for 20 min, measuring fluorescence (ex:470, em:520) every minute."

