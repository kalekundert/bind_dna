#!/usr/bin/env bash
set -euo pipefail

sw purex f62 f63 f65 f66 -v5 -I |
sw step "Add 0.5 µL RNase A/RNase T1 cocktail to each reaction." |
sw step "Divide each reaction into 2x 2.5 µL aliquots." |
sw step "Incubate the aliquots at 37°C and 65°C respectively for 15 min." |
sw gel natbr 8
