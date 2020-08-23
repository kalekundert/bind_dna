#!/usr/bin/env bash
set -euo pipefail

stepwise custom A | stepwise <<EOF | stepwise custom D
- B
- C
EOF

