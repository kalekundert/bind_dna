#!/usr/bin/env bash
set -euo pipefail

sw reaction \
  "nuclease-free water ;          ; to 6.00 µL ; +" \
  "RNase H buffer      ;      10x ;    0.60 µL ; +" \
  "protected mRNA      ;     1 µM ;    1.00 µL ; -" \
  "RNase H             ; 5000U/mL ;    0.12 µL ; +" \
