#!/usr/bin/env bash
set -euo pipefail

sw dilute \
  s4:0.164OD \
  s5:0.123OD \
  s36:0.243OD \
  s37:0.152OD \
  s38:0.128OD \
  s16:0.252OD \
  s22:0.222OD \
  sz213:0.256OD \
  sz214:0.238OD \
  sz215:0.210OD \
  -c 0.10OD \
  -V 146.25 \
  -d '1x NM'
