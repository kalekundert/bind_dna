#!/usr/bin/env bash
set -euo pipefail

sw gel bolt/gfp-mrna f85,f85 -S |
sw laser red |
sw step "Cut the gel to separate the two identical lanes, and keep one half as a −transfer control." |
sw northern/transfer_semidry |
sw step "Image the gel (as above) every 30 min, to check how much of the RNA has been transferred."
