#!/usr/bin/env bash
set -euo pipefail

# Template:
# - Want a control that I know works.
# - I'll start with purified plasmid, then repeat with colonies.
# - I would just start with colonies, but I don't have any positive controls.

sw pcr \
  p181,o246,o247 \
  p181,o246,o250 \
  p181,o246,o251 \
  p181,o246,o252 \
  p181,o246,o253 |

sw e_gel
