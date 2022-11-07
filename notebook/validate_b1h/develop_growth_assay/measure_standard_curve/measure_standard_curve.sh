#!/usr/bin/env bash
set -euo pipefail

# Compare 2-plasmid vs 1-plasmid for standard curve: which is more linear?

# Target  Strain  Plasmids
# ------  ------  ----------
# 2TGG    s4      p167, p166
# 2AAA    s5      p168, p166
# 2TAG    s36     p210, p166
# 2TTG    s37     p211, p166
# 2CAG    s38     p212, p166
# ------  ------  ----------
# 2TGG    s16     p189
# 2AAA    s22     p202
# 2TAG    sz213   p213
# 2TTG    sz214   p214
# 2CAG    sz215   p215

sw b1h/od_kinetic s4 s5 s36 s37 s38 s16 s22 sz213 sz214 sz215
