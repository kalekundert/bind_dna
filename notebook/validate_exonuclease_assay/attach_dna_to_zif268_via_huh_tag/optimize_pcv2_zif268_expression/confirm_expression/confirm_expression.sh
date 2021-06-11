#!/usr/bin/env bash
set -euo pipefail

sw step "Transform p106 and p107 into Lemo21(DE3) cells." |
sw step "Label the products: s8, s9" |

sw cond confirm_expression_cond.xlsx \
  -k "NI: non-induced" \
  -k "I: induced" \
  -k "CL: cleared lysate" \
  -k "FT: flow-through" \
  -k "W1: wash #1" \
  -k "W1: wash #2" \
  -k "E1: eluate #1" \
  -k "E2: eluate #2" \
  |

sw ni_nta/denaturing_buffers |
sw ni_nta/grow_miniprep |
sw ni_nta/miniprep_denaturing |

sw gel bolt/max 16
