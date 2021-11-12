#!/usr/bin/env bash
set -euo pipefail

sw step "Resuspend each colony in 20 µL EB." |

sw pcr \
  p184,o2,o262 \
  p184,o2,o185 \
  p184,o2,o188 \
  p185,o2,o262 \
  p185,o2,o185 \
  p185,o2,o188 \
  p186,o2,o262 \
  p186,o2,o185 \
  p186,o2,o188 \
|

sw e_gel

