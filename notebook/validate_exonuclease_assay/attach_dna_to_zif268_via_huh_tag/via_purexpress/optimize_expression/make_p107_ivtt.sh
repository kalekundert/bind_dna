#!/usr/bin/env bash
set -euo pipefail

# Prepare p107 for IVTT expression.  Specifically, this means:
# - Phenol-chloroform extraction to remove RNase activity.
# - Lyophilize to increase concentration.

sw step "Transform 0.5 µL p107 (124 ng/µL) with 10 µL MACH1 cells." |
sw step "Miniprep directly from the resulting lawn.  Elute in 200 µL water." |
sw phenol |
sw lyo -c '100 nM (408.9 ng/µL)'
