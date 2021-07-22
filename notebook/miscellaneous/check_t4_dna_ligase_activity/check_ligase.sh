#!/usr/bin/env bash
set -euo pipefail

sw reaction \
    -s "T4 DNA ligase control" \
    "water                          ;           ; to 20 µL ; +" \
    "T4 DNA ligase buffer           ;       10x ;     2 µL ; +" \
    "Lambda DNA/HindIII DNA marker  ; 500 ng/µL ;     1 µL ; +" \
    "T4 DNA ligase                  ;    1 U/µL ;     1 µL ; -" \
    -n 2 \
    |

sw note "https://tinyurl.com/24e725tj" |

sw step "Setup another copy of the above reaction with fresh T4 DNA ligase buffer." |

sw step "Incubate at room temperature for 10 min." |

sw e_gel
