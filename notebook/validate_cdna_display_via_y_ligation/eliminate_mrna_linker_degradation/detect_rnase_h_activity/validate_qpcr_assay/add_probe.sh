#!/usr/bin/env bash
set -euo pipefail

sw step "Add 0.4 µL annealed probe to each RNase H sample." |
sw note "1 pmol probe, ≈100 nM final probe conc" "sample" |
sw step "Incubate at 37°C for 20 min." |
sw note "No heat inactiavtion; see expt #87"
