#!/usr/bin/env bash
set -euo pipefail


sw zap |

# Controls:
# - −linker: to see the difference between coupled/uncoupled.
# - −mRNA: to see which fluorescent bands are always there.
# 
# I did salt titrations when trying to reproduce [Reyes2021], but my goal this 
# time is just to reproduce [Doshi2014], so I'll worry about optimization 
# later.
sw cond validate_puro.xlsx |

# IVTT volume:
# - 12 µL is max that can be loaded on desalting column.
# - I'll go with 10 µL: round number close to max, should se something.
sw ivtt f146 f145 -p purex/doshi2014 |

sw step "Incubate on ice for 10 min." |

# Salt concentration
# - [Doshi2014] calls for:
#   - 2.5 µL 1M MgCl₂
#   - 12.5 µL 2M KCl
# - Total volume: 40 µL
# - Final concentrations:
#   - MgCl₂: 62.5 mM
#   - KCl: 625 mM
#   - These are comparable with [Reyes2021].
#
# - My stocks:
#   - MgCl₂: 2M
#   - KCl: 3M
# - My volumes:
#   - total: 13.15
#   - 0.41 µL 2M MgCl₂
#   - 2.74 µL 3M KCl
sw step "Add 3.15 µL of the following salt solution to each reaction:
~3 µL 2M MgCl₂
~20 µL 3M KCl" |

sw step "Incubate as follows:~RT for 10 min~−20°C overnight" |

# Desalting
# - Can't run high-K solution with SDS PAGE.
# - [Doshi2014] calls for bead purification, but that will lose the −linker 
#   control.
# - Zeba desalting columns:
#   - 7 kDa MWCO
#   - Nb-GFP is 14.8 kDa
#   - Should be ok
sw zeba |

sw gel bolt/ivtt/nb-gfp 3
