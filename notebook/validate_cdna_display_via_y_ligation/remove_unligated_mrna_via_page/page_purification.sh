#!/usr/bin/env bash
set -euo pipefail

sw make f111 |
sw zymo-100 |

# Just the anneal/ligate steps
sw make f113 |

# - Start with pre-cast gel.
# - Ask: How much can I load in a single lane?
#   - Max: 12 lanes → 1/12 reaction (≈10 µg)
#   - serial dilution by 2 → 2**12 ≈ 1000 → 10 ng
#
sw gel page -I

