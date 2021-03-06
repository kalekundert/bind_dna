#!/usr/bin/env bash
set -euo pipefail

# 1 µM = 1 pmol/µL
#sw zap |
#sw smart_mmlv.py f11 o214 -m rna,primer |
#sw serial 20µL 50nM / 10 6 -m "RT product" -d 'nuclease-free water' -0 |
#sw qpcr 'RT dilution,NTC,NRT' o214 o215 -a 60 -l 88 -n 24 -m primers -V 5

sw zap |
sw cond validate_rt_qpcr_conditions.xlsx |
sw serial 10µL 100nM / 10 6 -0 -m f11 -d "nuclease-free water" |
sw smart_mmlv.py f11 o214 -n 7 -m primer -t 0.1 -T 100 |
sw qpcr 'RT products' o214 o215 -a 'optimal Ta' -l 88 -n 8 -v 66 -m primers -T '5 nM'
