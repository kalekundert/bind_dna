#!/usr/bin/env bash
set -euo pipefail

sw zap |

sw cond check_salt_degradation_cond.xlsx |

# mRNA concentration:
# - In IVTT reaction: 3.3 µL × 1000 nM / 10 µL = 330 nM
# - I have aliquots that are about that.

sw cdna/couple f89 -C 330nM |

sw gel bolt/gfp-mrna 4
