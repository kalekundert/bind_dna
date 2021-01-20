#!/usr/bin/env bash
set -euo pipefail

# Steps:
# - titrate RNase H
# - add probe
# - qPCR

sw zap |
sw cond validate_qpcr_assay_conditions.xlsx |
sw anneal 'probe RNA (f11)' 'optimal target oligo' -c 2 -C 10 -v 6 |
sw serial 10 100U/mL to 0.1 7 -0 -m 'RNase H' -d '1x RNase H reaction buffer' |
sw add_probe.sh |
sw smart_mmlv.py 'crude samples + probes' o214 -n 8 -m primer -t 0.5 -T 100 |
sw qpcr 'crude RT products' o214 o215 -a 'optimal Ta' -l 88 -n 10 -v 66 -m primers

