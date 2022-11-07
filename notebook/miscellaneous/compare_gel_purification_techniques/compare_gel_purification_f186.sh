#!/usr/bin/env bash
set -euo pipefail

# DNA quantity:
# - 900 ng at 20 ng/µL works out to 45 µL.
# - Given a ≈50 µL miniprep, that leaves 5 µL for Nanodropping/running a gel.
#
# Wait, no reason to use freshly miniprepped DNA.  But I will run gel of 
# miniprepped DNA, and if I use the same, I will be able to calculate true 
# yield.  
sw digest p236 NotI-HF,HindIII-HF -d 0.9 |

sw step "Label the product: f186" |

sw gel agarose/purify/2kb f186 -n3 -v 15 -l 18 |

sw freeze_and_squeeze 2.4kb |

sw electroelute |

sw step "Purify the desired band using the Qiagen QIAEX II gel extraction kit."

sw step "Measure yields by Nanodrop and by 1% E-Gel EX."
