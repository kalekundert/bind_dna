#!/usr/bin/env bash
set -euo pipefail

klab_reverse_translate -S GSR -D GGGTCTCGC -R BsaI
klab_reverse_translate -S ERFP -D GAACGTTTTCCA -R XmnI
