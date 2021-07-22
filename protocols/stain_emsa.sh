#!/usr/bin/env bash
set -euo pipefail

# Simultaneously visualize protein and Cy5-labeled DNA in native PAGE gels.

sw warn "the 'stain_emsa' protocol hasn't been tried before" |

# Soak time:
# - The SYPRO Orange manual <https://tinyurl.com/r7p6juer> calls for soaking in 
#   0.05% SDS for 30 min as part of the process for removing Triton X-100 from 
#   gels.
# - This handbook <https://tinyurl.com/2f59nzf8> calls for soaking in 0.05% SDS 
#   in order to stain native gels with SYPRO Orange/Red, but doesn't specify a 
#   time.  I think it's reasonable to use the time from the manual, though.
sw step "Soak gel in 0.05% SDS for 30 min." |

sw stain sypro-orange -I |

sw laser blue red
